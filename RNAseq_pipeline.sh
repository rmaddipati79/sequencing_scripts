#!/bin/bash

#####################################################################################
# bash script for RNA sequencing QC, Alignment, Data munging, and Raw Count generation on PennHPC
#
# by David Balli, Stanger laboratory at the University of Pennsylvania, Abramson Family Cancer Research Center
#
#
# tools used in this pipe:
# FastQC version 0.11.2 - http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# Fastx_toolkit version 0.0.6 - http://hannonlab.cshl.edu/fastx_toolkit/download.html
# STAR version 2.3.0e - https://code.google.com/p/rna-star/
# Samtools versions 0.1.18 and 0.1.19 - http://samtools.sourceforge.net
# HTSeq-count version 0.6.1 - http://www-huber.embl.de/users/anders/HTSeq/doc/index.html
#
# Differential gene expression anaysis is performed in DESeq2 (v2.X)
#
#####################################################################################
#
# useage:
# bsub -q max_mem64 -n 15 -e e.RNAseq -o o.RNASeq sh ~/scripts/RNAseq_pipeline.sh ARGVs
#
#####################################################################################
#
# notes:
#
# $1 # Forward fastq read - be sure to include .fq, can use .fq.gz files
# $2 # Reverse fastq read - be sure to include .fq, can use .fq.gz files
# $3 # output base name for files
#
#
#####################################################################################
# run Fastqc to generate QC metrics for each individual fq read
mkdir fastqc/
bsub -q plus -n 10 -e e.fastq -o o.fastq sh ~/scripts/fastqc.sh

# load module STAR v2.3.0e from PMACS modules
# Align to mm10 UCSC
if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh
fi

module load STAR-2.3.0e

STAR --genomeDir ~/mm10_STAR/ --readFilesIn $1 $2 --runThreadN 16 --outSAMattributes Standard 

# convert sam to bam file, sorting by name (using -n flag)
# load samtools v0.1.19 for faster/better sorting

module load samtools-0.1.19

samtools view -uS Aligned.out.sam -b | samtools sort -n -@ 8 -m 800000000 - $3.sorted

module load samtools-0.1.19

samtools view $3.sorted.bam > $3.sorted.sam


# for counts generation - need a SAM file because:
# the python environment for HTSeq count requires the module "pysam"
# currently (2014-9-12), the cluster is not installing pysam correctly

# run HTseq count in fresh submission as it takes several hours for 65.0x10^6 RNAseq reads
# very low memory usage so use normal job queue  

bsub -n 20 -e e.counts -o o.counts sh ~/scripts/htseqcount_nonempty.sh $3.sorted.sam $3

# should generate tab delimited text file of "GENEID   # reads"
# once ht-seq count is done convert $3.sorted.sam to BAM file for long term storage (saves a lot of disc space)
# once all count files have been generated - scp counts_files/ to laptop for DESeq2 data analysis




