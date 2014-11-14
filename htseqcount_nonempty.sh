#!/bin/bash


##################################################################################################################################
# by David Balli, Stanger laboratory  
#
#usage:
# bsub -q max_mem64 -n 4 -T 16 -e e.counts -o o.counts sh ~/scripts/htseqcount_nonempty.sh ARGVs
#
# will generate tab-delimited file of "GENEID\tCOUNT"
#
# output is used in DESeq2 differential analysis
#
##################################################################################################################################
#
# $1 - input sam file
# $2 - output base name
# this counts the reads from each sam file
# use genes.gtf annotation file in mm10
# using '-m intersection-nonempty' flag to reduce conservativness of count calls
#
##################################################################################################################################


htseq-count -m intersection-nonempty -i gene_name -a 10 --stranded=no $1 ~/bin/mm10_genes_gtf/genes.gtf > $2

mkdir ../counts_files
cp $2 ../counts_files