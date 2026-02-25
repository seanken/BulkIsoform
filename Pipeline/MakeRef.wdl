version 1.0

task BuildSTAR {
    input {
        File fa
        File gtf
    }
    
    command <<<
        STAR --runMode genomeGenerate --genomeDir STAR_ref --genomeFastaFiles ~{fa} --sjdbGTFfile ~{gtf}
        mkdir genes
        cp ~{gtf} genes/genes_new.gtf
    >>>
    
    output {
        File star_ref = "STAR_ref"
        File genes_dir = "genes"
    }
    
    runtime {
        docker: "quay.io/biocontainers/star:2.7.10a--h9ee0642_0"
    }
}

task BuildRSEM {
    input {
        File fa
        File gtf
    }
    
    command <<<
        mkdir RSEM
        rsem-prepare-reference --gtf ~{gtf} --star ~{fa} RSEM/ref
    >>>
    
    output {
        File rsem_ref = "RSEM"
    }
    
    runtime {
        docker: "quay.io/biocontainers/rsem:1.3.3--pl5.22.0_2"
    }
}

task BuildSalmon {
    input {
        File fa
        File rsem_ref
    }
    
    command <<<
        mkdir Salmon
        cat ~{rsem_ref}/ref.idx.fa ~{fa} > comb.fa
        grep "^>" ~{fa} | cut -d " " -f 1 | sed -e 's/>//g' > decoys.txt
        salmon index -t comb.fa -i Salmon/ref.idx --decoys decoys.txt
        rm comb.fa
    >>>
    
    output {
        File salmon_ref = "Salmon"
    }
    
    runtime {
        docker: "combinelab/salmon:1.10.3"
    }
}

task BuildRefFlat {
    input {
        File gtf
    }
    
    command <<<
        gtfToGenePred -genePredExt -geneNameAsName2 -ignoreGroupsWithoutExons ~{gtf} temp.txt
        awk 'BEGIN { OFS="\t"} {print $12, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}'  temp.txt > refflat.txt
    >>>
    
    output {
        File refflat = "refflat.txt"
    }
    
    runtime {
        docker: "quay.io/biocontainers/ucsc-gtftogenepred:469--h664eb37_1"
    }
}

task CopyToOutdir {
    input {
        File star_ref
        File genes_dir
        File rsem_ref
        File salmon_ref
        File refflat
        String outdir
    }
    
    command <<<
        mkdir -p ~{outdir}
        cp -r ~{star_ref} ~{outdir}/STAR_ref
        cp -r ~{genes_dir} ~{outdir}/genes
        cp -r ~{rsem_ref} ~{outdir}/RSEM
        cp -r ~{salmon_ref} ~{outdir}/Salmon
        mkdir -p ~{outdir}/RefFlat
        cp ~{refflat} ~{outdir}/RefFlat/refflat.txt
    >>>
    
    output {
        File final_outdir = outdir
    }
}

workflow MakeRef {
    input {
        File fa
        File gtf
        String outdir="ref"
    }
    
    call BuildSTAR {
        input:
            fa = fa,
            gtf = gtf
    }
    
    call BuildRSEM {
        input:
            fa = fa,
            gtf = gtf
    }
    
    call BuildSalmon {
        input:
            fa = fa,
            rsem_ref = BuildRSEM.rsem_ref
    }
    
    call BuildRefFlat {
        input:
            gtf = gtf
    }
    
    call CopyToOutdir {
        input:
            star_ref = BuildSTAR.star_ref,
            genes_dir = BuildSTAR.genes_dir,
            rsem_ref = BuildRSEM.rsem_ref,
            salmon_ref = BuildSalmon.salmon_ref,
            refflat = BuildRefFlat.refflat,
            outdir = outdir
    }
    
    output {
        File output_directory = CopyToOutdir.final_outdir
    }
}
