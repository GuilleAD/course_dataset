# Overview of this repo

> This repo houses all of the code and analysis conducted for the 2026 iteration of the DIYtranscriptomics course

## Abstract
Some stuff about your science that the code in this repo was intended to tackle.  If you have a published or draft manuscript, this can literally be your abstract.

## Code and files
> The locations of the core components of this repo are outlined in the file system map below. In short, there are the following main directories:

 - /analysis/code - contains the RScripts and Rmarkdown (.rmd) files that were used to generate the supplementary code file included in the manuscript the pipeline.  All code for manuscript figures based on RNAseq data are included in this document.  Also contains the study design files that describe each sample.
 - /readMapping/Kallisto_output_cDNA - includes Kallisto outputs and log files for mapping of raw fastq files from the course to the human reference transcriptome that includes all annotated cDNAs as of January 2026 (Release 115; Homo_sapiens.GRCh38.cdna.all.fa).  Raw fastq. files were omitted from this repo due to restrictions inherent to GitHub and can be found on the Sequence Read Archive (SRA) under accession XXXXX.
  - /readMapping/Kallisto_output_ncRNA- Kallisto outputs and log files for mapping of raw fastq files to the human reference file for non-coding RNAs as of January 2026 (Release 115; Homo_sapiens.GRCh38.ncrna.fa)
 - /QC_results/fastqc - includes outputs from running [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) on all raw fastq files used in the course.
 - /QC_results/multiqc- Contains all outputs from running [multiqc](https://multiqc.info/), including the html summary report that summarizes fastqc and kallisto results from mapping raw fastq files to the human reference transcriptome.

