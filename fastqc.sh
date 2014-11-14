#!/bin/bash

# usage: bsub -q plus -n 10 -e e.fastqc -o o.fastqc sh fastqc.sh 
# or
# being called from RNAseq_pipeline.sh

if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh
fi

module load java-sdk-1.7.0

fastqc *.fq.gz -f fastq -o fastqc/
fastqc *.fq -f fastq -o fastqc/

## QC assement of each fq files using Fastqc
## fastqc/ is sub-directory in working seq directory
## repeat for each X.fq.gz or X.fq file
