# Bulk Isoform Nextflow Pipeline

This pipeline is built to take in raw fastq files from paired end bulk RNA-seq data and perform numerous basic analysis on it, with the goal of producing information used for Isoform analysis. It requires Java, nextflow, and (optionally) Conda (without Conda several other tools will need to be installed seperately). 

## Set Up

You need to have Java installed (v8 or later) and some version of Conda, as well as nextflow (https://www.nextflow.io/). To install this, download it from git hub with:
```
git clone https://github.com/seanken/BulkIsoform.git
```

You also need to install the requirements, see below--the easiest way is to start a conda environment and run 'source scripts/makeConda.sh'. Once you do that, however, you are then ready to run!


## Requirements

The pipeline is set up so that all you need to run it is Conda, nextflow, and java as mentioned above. It is, however, possible to run the pipeline without conda. In that case, however, the number of requirements greatly increases. In that case, in addition to Java and nextflow, you will need the following installed/on the PATH:
1) regtools 
2) samtools
3) bedtools
4) STAR
5) subread (for featureCounts)
6) salmon
7) RSEM
8) gtftogenepred (for making reference only)

If you have conda installed you simply need to create a new environment and use the script scripts/makeConda.sh to install everything except RSEM. We plan on adding alternative methods to get the proper environment later.

## Generating Reference

In order to run the pipeline, one needs a gtf, a STAR reference, a salmon reference, a refflat file, and a RSEM reference. It is possible to generate each of them seperately, but we also offer an option to generate a joint reference. To do this, one needs a gtf file and a fasta file from the genome/transcriptome of the species of interest. 

```
nextflow=/path/to/nextflow
pipeline=/path/to/
gtf=/path/to/gtf
fa=/path/to/fasta
output=/path/to/output/directory

$nextflow $pipeline --gtf $gtf --fa $fa --outdir $output
```

This will create a reference in the output directory that we can use for the pipeline.

## How To Run Pipeline

```
nextflow=/path/to/nextflow
pipeline=/path/to/BulkPipeline.nf

$nextflow $pipeline [options]
```

The pipeline has a mix of optional and required options. The options:

#### --fq1 
The fastq for read 1. Can be a comma seperated list. Must all by gzipped.
#### --fq2 
The fastq for read 2. Can be a comma seperated list. Must all by gzipped.
#### --ref_comb
A combined reference generated as above.
#### --star_ref 
The STAR reference to be used, not needed if combined reference is given.
#### --gtf 
The GTF to use, corresponding to the STAR reference, not needed if combined reference is given.
#### --outdir 
The name of the out directory to publish results to.
#### --refflat 
A ref flat file to use with PICARD. Should correspond to the GTF. Needed for QC step. Not needed if combined reference is given.
#### --salmon_ref
A Salmon reference, not needed if combined reference is given. Suggested to use decoys.
#### --rsem_ref
A RSEM reference, not needed if combined reference is given.

## To Be Added

This is an initial version of the pipeline. There are many things to be added, in particular an updated pipeline to run multiple samples at once using a sample sheet and updates to run with different strandedness (currently has read 1 to be antisense and read 2 to be sense).
