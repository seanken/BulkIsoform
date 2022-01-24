params.fa
params.gtf
params.outdir
params.refflat="$projectDir/scripts/MakeRefFlat.sh"

process BuildSTAR //Builds reference for STAR and copies over GTF into reference structure being created
{
publishDir "${params.outdir}", mode: 'move'
input:
path fa, stageAs:"genome.fa" from params.fa
path gtf, stageAs:"genes.gtf" from params.gtf

output:
path "STAR_ref" into STAR
path "genes" into GTF

'''
STAR --runMode genomeGenerate --genomeDir STAR_ref --genomeFastaFiles genome.fa --sjdbGTFfile genes.gtf
mkdir genes
cp genes.gtf genes/genes_new.gtf
'''

}


process BuildRSEM //Builds RSEM reference with STAR
{
publishDir "${params.outdir}", mode: 'copy'
input:
path fa, stageAs:"genome.fa" from params.fa
path gtf, stageAs:"genes.gtf" from params.gtf

output:
path "RSEM" into RSEM

'''
mkdir RSEM
rsem-prepare-reference --gtf genes.gtf --star genome.fa RSEM/ref
'''

}



process BuildSalmon //Uses RSEM reference to build Salmon reference with full genome decoy
{
publishDir "${params.outdir}", mode: 'move'

input:
path fa, stageAs:"genome.fa" from params.fa
path RSEM, stageAs:"RSEM" from RSEM

output:
path "Salmon" into SalmonRef

'''
mkdir Salmon
cat RSEM/ref.idx.fa genome.fa > comb.fa
grep "^>" genome.fa | cut -d " " -f 1 | sed -e 's/>//g' > decoys.txt
salmon index -t comb.fa -i Salmon/ref.idx --decoys decoys.txt
rm comb.fa
'''

}


process BuildRefFlat //Creates RefFlat file for use with PICARD Tools
{
publishDir "${params.outdir}/RefFlat", mode: 'move'
input:
path params.gtf, stageAs:"genes.gtf" from params.gtf
path refflat, stageAs:"MakeRefFlat.sh" from params.refflat

output:
path "refflat.txt" into RefFlatFile

'''
source MakeRefFlat.sh genes.gtf refflat.txt
'''

}
