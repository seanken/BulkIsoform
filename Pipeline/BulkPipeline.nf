//input parameters
params.fq1 //The first fastq, can be a list, should all be gzipped 
params.fq2 //The second fastq, can be a list, should all be gzipped
params.star_ref //STAR reference
params.gtf //GTF reference corresponding to STAR reference
params.refflat //ref flat file, later extend pipeline to make on its own

params.picard="$projectDir/jars/picard.jar" //PICARD jar file
params.makeGene="$projectDir/scripts/makeGene.sh" //bash script to make gtf into genes bed file



//Maps reads to a reference with STAR
process MapReads
{
input:
env read1 from params.fq1
env read2 from params.fq2
path star_ref, stageAs:"ref" from params.star_ref

output:
path "mapped.bam" into mapped_bam
path "mapped.bam.bai" into mapped_bam_bai


'''
STAR --genomeDir star_ref --readFilesIn $read1 $read2 --outSAMattributes NH HI AS nM  --outSAMtype BAM SortedByCoordinate --outFileNamePrefix results --readFilesCommand zcat --outStd BAM_SortedByCoordinate > mapped.bam

samtools index mapped.bam
'''

}


//Count reads, both with and without introns counted, assumes antisense reads
process CountReads
{
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
input:
path pic_jar, stageAs:"picard.jar" from params.picard
path mapped_bam, stageAs:"mapped.bam" from mapped_bam
path mapped_bam, stageAs:"mapped.bam.bai" from mapped_bam_bai
path refflat, stageAs:"refFlat.txt" from params.refflat

output:
path "output.QC.txt" into Picard_QC


'''
java -jar picard.jar CollectRnaSeqMetrics I=mapped.bam O=output.QC.txt STRAND=FIRST_READ_TRANSCRIPTION_STRAND REF_FLAT=refFlat.txt
'''

}


//Uses regtools to extract junction level information from the bam files
process ExtractJunctions
{

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
input:
path reg_jun, stageAs:"regtools.junc" from regtools_junc
path gtf, stageAs:"genes.gtf" from params.gtf
path makeGene, stageAs:"makeGene.sh" from params.makeGene

output:
path "juncfile.genes.bed" into ann_juncs

'''
makeGene.sh genes.gtf genes.bed
bedtools intersect -S -a regtools.junc -b gene.bed -wa -wb > juncfile.genes.bed
'''
}




