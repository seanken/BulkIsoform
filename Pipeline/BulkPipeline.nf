//input parameters
params.fq1 //The first fastq, can be a list, should all be gzipped 
params.fq2 //The second fastq, can be a list, should all be gzipped
params.outdir //the out directory to publish to
params.ref_comb //Combined referencei
params.gtf="${params.ref_comb}/genes/genes_new.gtf" //GTF reference corresponding to STAR reference
params.refflat="${params.ref_comb}/RefFlat/refflat.txt" //ref flat file, later extend pipeline to make on its own
params.salmon_ref="${params.ref_comb}/Salmon/ref.idx" //Salmon reference
params.rsem_ref="${params.ref_comb}/RSEM/ref" //RSEM reference
params.star_ref="${params.ref_comb}/STAR_ref" //STAR reference
params.isoMethod="Salmon" //Either RSEM or Salmon or Both

params.picard="$projectDir/jars/picard.jar" //PICARD jar file
params.makeGene="$projectDir/scripts/makeGene.sh" //bash script to make gtf into genes bed file
params.geneToTrans="$projectDir/scripts/makeGeneToTrans.sh"


//Maps reads to a reference with STAR
process MapReads
{
publishDir "${params.outdir}/Bam", mode: 'rellink'

input:
env read1 from params.fq1
env read2 from params.fq2
path star_ref, stageAs:"ref" from params.star_ref

output:
path "mapped.bam" into mapped_bam
path "mapped.bam.bai" into mapped_bam_bai


'''
STAR --genomeDir ref --readFilesIn $read1 $read2 --outSAMattributes NH HI AS nM  --outSAMtype BAM SortedByCoordinate --outFileNamePrefix results --readFilesCommand zcat --outStd BAM_SortedByCoordinate > mapped.bam

samtools index mapped.bam
'''

}


//Count reads, both with and without introns counted, assumes antisense reads
process CountReads
{
publishDir "${params.outdir}/Counts", mode: 'rellink'
input:
path gtf, stageAs:"genes.gtf" from params.gtf
path mapped_bam, stageAs:"mapped.bam" from mapped_bam

output:
path "counts.exons.txt" into exon_counts
path "counts.introns.txt" into with_intron_counts


'''
featureCounts -p -s 2 -t exon -g gene_name -a genes.gtf -o counts.exons.txt mapped.bam 
featureCounts -p -s 2 -t gene -g gene_name -a genes.gtf -o counts.introns.txt mapped.bam 
'''

}

//Runs Picard too collect RNASeq Metrics
process PicardQC
{

publishDir "${params.outdir}/QC", mode: 'rellink'
input:
path pic_jar, stageAs:"picard.jar" from params.picard
path mapped_bam, stageAs:"mapped.bam" from mapped_bam
path mapped_bam, stageAs:"mapped.bam.bai" from mapped_bam_bai
path refflat, stageAs:"refFlat.txt" from params.refflat

output:
path "output.QC.txt" into Picard_QC


'''
java -jar picard.jar CollectRnaSeqMetrics I=mapped.bam O=output.QC.txt STRAND=SECOND_READ_TRANSCRIPTION_STRAND REF_FLAT=refFlat.txt
'''

}


//Uses regtools to extract junction level information from the bam files
process ExtractJunctions
{

publishDir "${params.outdir}/Junctions", mode: 'rellink'
input:
path mapped_bam, stageAs:"mapped.bam" from mapped_bam
path mapped_bam, stageAs:"mapped.bam.bai" from mapped_bam_bai

output:
path "regtools.junc" into regtools_junc

'''
regtools junctions extract -a 8 -m 50 -M 500000 -s 2 -o regtools.junc mapped.bam
'''

}

//annotates the junc files with gene of origin information
process AnnotateJunctions
{

publishDir "${params.outdir}/AnnJunctions", mode: 'rellink'

input:
path reg_jun, stageAs:"regtools.junc" from regtools_junc
path gtf, stageAs:"genes.gtf" from params.gtf
path makeGene, stageAs:"makeGene.sh" from params.makeGene

output:
path "juncfile.genes.bed" into ann_juncs

'''
source makeGene.sh genes.gtf genes.bed
bedtools intersect -S -a regtools.junc -b genes.bed -wa -wb > juncfile.genes.bed
'''
}



process RunRSEM
{
publishDir "${params.outdir}", mode: 'rellink'

input:
env RSEM from params.rsem_ref
env read1 from params.fq1
env read2 from params.fq2

output:
path "RSEM_out" into RSEM_count

when:
params.isoMethod!="Salmon"

'''
mkdir RSEM_out
rsem-calculate-expression --star --star-gzipped-read-file --no-bam-output --strandedness reverse --paired-end --estimate-rspd $read1 $read2 $RSEM RSEM_out/Samp
'''
}


process RunSalmon
{
publishDir "${params.outdir}", mode: 'rellink'

input:
path SALMON, stageAs:"ref" from params.salmon_ref
env read1 from params.fq1
env read2 from params.fq2
path gtf, stageAs:"genes.gtf" from params.gtf
path geneToTrans, stageAs:"makeGeneToTrans.sh" from params.geneToTrans

output:
path "Salmon_out" into Salmon_count

when:
params.isoMethod!="RSEM"

'''
mkdir Salmon_out
read1_new=$(echo $read1 | sed 's/,/ /g')
read2_new=$(echo $read2 | sed 's/,/ /g')
salmon quant -i ref -l A --posBias --seqBias --gcBias --validateMapping -1 $read1_new -2 $read2_new -o Salmon_out
echo Make tx to gene file
source makeGeneToTrans.sh genes.gtf Salmon_out/trans2gene.txt
'''
}

