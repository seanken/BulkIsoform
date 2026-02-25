version 1.0

workflow BulkPipeline {
    input {
        Array[File] fq1
        Array[File] fq2
        String outdir
        String ref_comb
        File gtf = ref_comb + "/genes/genes_new.gtf"
        File refflat = ref_comb + "/RefFlat/refflat.txt"
        File salmon_ref = ref_comb + "/Salmon/ref.idx"
        File rsem_ref = ref_comb + "/RSEM/ref"
        File star_ref = ref_comb + "/STAR_ref"
        String isoMethod = "Salmon"
        Int s = 2
        String gene_name = "gene_name"
        Int numGibbsSamples = 0
    }

    call MapReads {
        input:
            fq1 = fq1,
            fq2 = fq2,
            star_ref = star_ref
    }

    call CountReads {
        input:
            gtf = gtf,
            mapped_bam = MapReads.mapped_bam,
            mapped_bam_bai = MapReads.mapped_bam_bai,
            strand = s,
            gene_name = gene_name
    }

    call PicardQC {
        input:
            mapped_bam = MapReads.mapped_bam,
            mapped_bam_bai = MapReads.mapped_bam_bai,
            refflat = refflat,
            strand = s
    }

    if (isoMethod == "RSEM") {
        call RunRSEM {
            input:
                rsem_ref = rsem_ref,
                fq1 = fq1,
                fq2 = fq2,
                strand = s
        }
    }

    if (isoMethod == "Salmon") {
        call RunSalmon {
            input:
                salmon_ref = salmon_ref,
                fq1 = fq1,
                fq2 = fq2,
                gtf = gtf,
                numGibbsSamples = numGibbsSamples
        }
    }

    output {
        File mapped_bam_out = MapReads.mapped_bam
        File mapper_qc = MapReads.mapper_qc
        File exon_counts = CountReads.exon_counts
        File exon_counts_sum = CountReads.exon_counts_sum
        File with_intron_counts = CountReads.with_intron_counts
        File with_intron_counts_sum = CountReads.with_intron_counts_sum
        
        File picard_qc = PicardQC.picard_qc
        File picard_dup = PicardQC.picard_dup
        File? rsem_out = RunRSEM.rsem_out
        File? salmon_out = RunSalmon.salmon_out
    }
}

task MapReads {
    input {
        Array[File] fq1
        Array[File] fq2
        File star_ref
    }

    command <<<
        STAR --genomeDir ~{star_ref} \
             --readFilesIn ~{sep=',' fq1} ~{sep=',' fq2} \
             --outSAMattributes NH HI AS nM \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix results \
             --readFilesCommand zcat \
             --outStd BAM_SortedByCoordinate \
             --outSAMunmapped Within > mapped.bam

        samtools index mapped.bam
    >>>

    output {
        File mapped_bam = "mapped.bam"
        File mapped_bam_bai = "mapped.bam.bai"
        File mapper_qc = "resultsLog.final.out"
    }

    runtime {
        docker: "quay.io/biocontainers/star:2.7.10a--h9ee0642_0"
    }
}

task CountReads {
    input {
        File gtf
        File mapped_bam
        File mapped_bam_bai
        Int strand
        String gene_name
    }

    command <<<
        featureCounts -p -s ~{strand} -t exon -g ~{gene_name} \
                      -a ~{gtf} -o counts.exons.txt ~{mapped_bam}
        featureCounts -p -s ~{strand} -t gene -g ~{gene_name} \
                      -a ~{gtf} -o counts.introns.txt ~{mapped_bam}
    >>>

    output {
        File exon_counts = "counts.exons.txt"
        File exon_counts_sum = "counts.exons.txt.summary"
        File with_intron_counts = "counts.introns.txt"
        File with_intron_counts_sum = "counts.introns.txt.summary"
    }

    runtime {
        docker: "quay.io/biocontainers/subread:2.0.3--h7132678_0"
    }
}

task PicardQC {
    input {
        File mapped_bam
        File mapped_bam_bai
        File refflat
        Int strand
    }

    String picard_strand = if strand == 2 then "SECOND_READ_TRANSCRIPTION_STRAND" else if strand == 1 then "FIRST_READ_TRANSCRIPTION_STRAND" else "NONE"

    command <<<
        echo ~{picard_strand}
        mkdir tmp
        CollectRnaSeqMetrics \
             I=~{mapped_bam} O=output.QC.txt \
             STRAND=~{picard_strand} REF_FLAT=~{refflat} TMP_DIR=tmp
        MarkDuplicates \
             I=~{mapped_bam} O=dup.bam M=metrics.dup.txt TMP_DIR=tmp
        rm dup.bam
    >>>

    output {
        File picard_qc = "output.QC.txt"
        File picard_dup = "metrics.dup.txt"
    }

    runtime {
        docker: "broadinstitute/picard:3.4.0"
    }
}


task RunRSEM {
    input {
        File rsem_ref
        Array[File] fq1
        Array[File] fq2
        Int strand
    }

    String rsem_strand = if strand == 2 then "reverse" else if strand == 1 then "forward" else "none"

    command <<<
        echo Hi
        echo ~{rsem_strand}
        mkdir RSEM_out
        rsem-calculate-expression --star --star-gzipped-read-file \
             --no-bam-output --strandedness ~{rsem_strand} --paired-end \
             --estimate-rspd ~{sep=',' fq1} ~{sep=',' fq2} \
             ~{rsem_ref} RSEM_out/Samp
    >>>

    output {
        File rsem_out = "RSEM_out"
    }

    runtime {
        docker: "quay.io/biocontainers/rsem:1.3.3--pl5262hdcf5f25_4"
    }
}

task RunSalmon {
    input {
        File salmon_ref
        Array[File] fq1
        Array[File] fq2
        File gtf
        Int numGibbsSamples = 0
    }

    command <<<
        mkdir Salmon_out
        
        # Build salmon command with optional Gibbs sampling
        SALMON_CMD="salmon quant -i ~{salmon_ref} -l A --posBias --seqBias --gcBias --validateMapping -1 ~{sep=' ' fq1} -2 ~{sep=' ' fq2} -o Salmon_out"
        
        if [[ ~{numGibbsSamples} -gt 0 ]]; then
            SALMON_CMD="$SALMON_CMD --numGibbsSamples ~{numGibbsSamples}"
        fi
        
        eval $SALMON_CMD
        
        echo Make tx to gene file
        
        # Generate trans2gene.txt
       
        grep -v \# ~{gtf} | awk '{if($3=="transcript"){print $0}}' | awk -F 'gene_name ' '{print $2}' | awk '{print $1}' | sed 's/\;//g' | sed 's/\"//g' > gene.txt
        grep -v \# ~{gtf} | awk '{if($3=="transcript"){print $0}}' | awk -F 'transcript_id ' '{print $2}' | awk '{print $1}' | sed 's/\;//g' | sed 's/\"//g' > trans.txt
        paste trans.txt gene.txt | sort | uniq > Salmon_out/trans2gene.txt
        rm gene.txt trans.txt
        
        # Generate trans2gene.ens.txt
        grep -v \# ~{gtf} | awk '{if($3=="transcript"){print $0}}' | awk -F 'gene_id ' '{print $2}' | awk '{print $1}' | sed 's/\;//g' | sed 's/\"//g' > gene.txt
        grep -v \# ~{gtf} | awk '{if($3=="transcript"){print $0}}' | awk -F 'transcript_id ' '{print $2}' | awk '{print $1}' | sed 's/\;//g' | sed 's/\"//g' > trans.txt
        paste trans.txt gene.txt | sort | uniq > Salmon_out/trans2gene.ens.txt
        rm gene.txt trans.txt
    >>>

    output {
        File salmon_out = "Salmon_out"
    }

    runtime {
        docker: "combinelab/salmon:1.10.3"
    }
}
