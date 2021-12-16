# Bulk Isoform Nextflow Pipeline

This pipeline is built to take in raw fastq files from paired end bulk RNA-seq data and perform numerous basic analysis on it, with the goal of producing information used for Isoform analysis. It requires Java, nextflow, and (optionally) Conda (without Conda several other tools will need to be installed seperately). 

## Set Up

You need to have Java installed (v8 or later) and some version of Conda, as well as nextflow (https://www.nextflow.io/). To install this, download it from git hub with:
```
git clone https://github.com/seanken/BulkIsoform.git
```

You are then ready to run!


## Requirements

The pipeline is set up so that all you need to run it is Conda, nextflow, and java as mentioned above. It is, however, possible to run the pipeline without conda. In that case, however, the number of requirements greatly increases. In that case, in addition to Java and nextflow, you will need the following installed/on the PATH:
1) regtools 
2) samtools
3) bedtools
4) STAR
5) featureCounts

If you have conda installed you simply need to create a new environment and use the script scripts/makeConda.sh to install everything. We plan on adding alternative methods to get the proper environment later.

## Generating Reference

In order to run this pipeline you need a STAR reference and a matching gtf file. See STAR documentation for details. To run the QC step you also need a ref flat file corresponding to the gtf (note we plan to add an option to generate the file if it is not present).

## How To Run Pipeline

```
nextflow=/path/to/nextflow
pipeline=/path/to/BulkPipeline.nf

$nextflow $pipeline [options]
```

The pipeline has a mix of optional and required options.
