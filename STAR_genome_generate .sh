#!/bin/bash

# using STAR for RNAseq alignment
# this script generates files need for STAR alignment to mm10 reference genome from UCSC
# http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/

# usage: bsub -q max_mem30 -n 10 -e e.star-generator -o o.star-generator sh ~/scripts/STAR_genome_generate.sh

if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh
fi
module load STAR-2.3.0e

STAR --runMode  genomeGenerate --genomeDir /home/dballi/mm10_STAR --genomeFastaFiles /home/dballi/mm10_index/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa --runThreadN 1