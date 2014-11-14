#!/bin/bash


#####################################################################################
# bash script for RNA sequencing Alignment, munging, GATK variant calling pipeline,
# Mutect/Somatic-Sniper Tumor-Normal variant comparrison and discovery of
# significantly mutated genes (SMGs) using MuSIC on PennHPC
#
# based on Broads GATK best practices:
# http://www.broadinstitute.org/gatk/guide/best-practices?bpm=DNAseq#data-processing-ovw
#
# by David Balli, Stanger laboratory at the University of Pennsylvania, Abramson Family Cancer Research Center
#
# tools used in this pipeline:
# FastQC version 0.11.2 - http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# Fastx_toolkit version 0.0.6 - http://hannonlab.cshl.edu/fastx_toolkit/download.html
# GATK version  3.2-2 - https://www.broadinstitute.org/gatk/download
# STAR version 2.3.0e - https://code.google.com/p/rna-star/
#
#####################################################################################
#
# useage:
# bsub -q max_mem64 -n 20 -e e.RNAseqVariant -o o.RNAseqVariant sh ~/scripts/RNASeq/GATK_RNAseq_Variant.sh ARGVs
#
#####################################################################################
# notes:
#
# ARGVs:
# $1 # mate 1 fastq - can be gz file  
# $2 # mate 2 fastq - can be gz file 
# $3 # output base name for files
#
#####################################################################################

# load java version 1.7.0 amd STAR-2.3.0e

if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh
fi

module load java-sdk-1.7.0
module load STAR-2.3.0e

# Alignment 1 with STAR 
STAR --genomeDir ~/mm10_STAR/ --readFilesIn $1 $2 --runThreadN 16 --outSAMattributes Standard 

# convert to bam and sort by genomic location
java -Xmx16g -jar ~/bin/AddOrReplaceReadGroups.jar \
SO=coordinate \
I=Aligned.out.sam \
O=$3.sorted.rg.bam \
RGLB=RNAseq_NYU \
RGPL=Illumina RGID=$3 \
RGPU=NYU_$3 \
RGSM=$3 \
VALIDATION_STRINGENCY=LENIENT

# mark pcr duplicates
java -Xmx16g -jar ~/bin/MarkDuplicates.jar \
I=$3.sorted.rg.bam \
O=$3.sorted.rg.ddup.bam \
M=$3.metrics.txt \
CREATE_INDEX=true \

### skipping this for now - weird error of "contig '*' not in referencer"
# split'n'trim and reassign mapping qualities
java -Xmx16g -jar ~/bin/gatk.jar \
-T SplitNCigarReads \
-R ~/mm10_STAR/mm10.fa \
-I $3.sorted.rg.ddup.bam \
-o $3.sorted.rg.ddup.split.bam \
-rf ReassignOneMappingQuality \
-RMQF 255 \
-RMQT 60 \
-U ALLOW_N_CIGAR_READS


# realignment calculation A
java -Xmx16g -jar ~/bin/gatk.jar \
-T RealignerTargetCreator \
-R ~/mm10_STAR/mm10.fa \
-U ALLOW_N_CIGAR_READS \
-o $3.output.intervals \
-I $3.sorted.rg.ddup.split.bam

# realignment calculation B
java -Xmx16g -jar ~/bin/gatk.jar \
-T IndelRealigner \
-R ~/mm10_STAR/mm10.fa \
-I $3.sorted.rg.ddup.split.bam \
-targetIntervals $3.output.intervals \
-U ALLOW_N_CIGAR_READS \
-o $3.sorted.rg.ddup.split.realigned.bam 

# base recalibration 1 (BQSR)
# we are using mm10 snp137 for known variants
java -Xmx16g -jar ~/bin/gatk.jar \
-T BaseRecalibrator \
-R ~/mm10_STAR/mm10.fa \
-I $3.sorted.rg.ddup.split.realigned.bam \
-L chr19 \
-knownSites ~/bin/mm10_dbsnp137.vcf \
-o $3.recal_data.table

# second pass - recommended by GATK
java -Xmx16g -jar ~/bin/gatk.jar \
-T BaseRecalibrator \
-R ~/mm10_STAR/mm10.fa \
-I $3.sorted.rg.ddup.split.realigned.bam \
-L chr19 \
-knownSites ~/bin/mm10_dbsnp137.vcf \
-BQSR $3.recal_data.table \
-o $3.post_recal_data.table

# finish BQSR 
java -Xmx16g -jar ~/bin/gatk.jar \
-T PrintReads \
-R ~/mm10_STAR/mm10.fa \
-I $3.sorted.rg.ddup.split.realigned.bam \
-BQSR $3.recal_data.table \
-o $3.sorted.rg.ddup.split.realigned.recal.bam

# calling variants with HaplotypeCaller
java -Xmx16g -jar ~/bin/gatk.jar \
-T HaplotypeCaller \
-R ~/mm10_STAR/mm10.fa \
-I $3.sorted.rg.ddup.split.realigned.recal.bam \
-recoverDanglingHeads \
-dontUseSoftClippedBases \
-stand_call_conf 20.0 \
-stand_emit_conf 20.0 \
-o $3.gatk.vcf 

# make copy of final bam file and remove various middle step bams 
cp $3.sorted.rg.ddup.split.realigned.recal.bam $3.final.bam &
rm $3.sorted.rg.ddup.split.realigned.recal.bam \
$3.sorted.rg.ddup.split.realigned.bam \
$3.sorted.rg.ddup.bam \
$3.sorted.rg.ddup.split.bam 

# convert aligned.out.sam to bam 
module load samtools-0.1.19
samtools view -bS Aligned.out.sam > $3.aligned.out.bam &
rm Aligned.out.sam 

mkdir ../RNAseqVariant
ln $3.final.bam $3.gatk.vcf $3.output.intervals ../RNAseqVariant
