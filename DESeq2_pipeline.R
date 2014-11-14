####################################################################################  
# bioinformatic analysis of RNA-seq data using DESeq2
# by David Balli, Stanger laboratory 
# 
# combined from DESeq2 manual and vignette 
# http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf
# and from Dave Wheeler's blog at Massey University
# http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
#
# analysis per htseq count -intersections-nonempty
# controlling for batch (different RNAseq prep libraries)
#
# notes: 
#
####################################################################################

library("DESeq2")
library("gplots")
library("grDevices")
library("rgl")
library('ggplot2')
library("genefilter")
library("RColorBrewer")
library("vsn")

setwd("~/path/to/working/directory/")
directory <- "/path/to/counts/directory/"

# can merge individual sample files (i.e. control 1, control 2, etc.)
sampleFiles <- grep("PD"list.files(directory),value=T)

# view sampleFiles
sampleFiles

# can designate different batches of samples (i.e. different sequencers,
# PE vs SE, different library preps (BATCH1 vs BATCH2))
sampleBatch <- c("Batch1","Batch1","Batch1","Batch1","Batch1","Batch1",
                 "Batch2","Batch2","Batch2","Batch2")

# set sampleConditions and sampleTable for experimental conditions
sampleCondition <- c("Control, Experimental")
sampleTable <- data.frame(sampleName = sampleFiles, 
                          fileName = sampleFiles, 
                          condition = sampleCondition,
                          Batch = sampleBatch)

# view sampleTable
sampleTable 

ddsHTseq <- DESeqDataSetFromHTSeqCount( sampleTable = sampleTable, 
                                       directory = directory, 
                                       design = ~condition) # + batch)

## view ddsHTseq - should give summary of class, data, etc.
ddsHTseq

colData(ddsHTseq)$condition <- factor(colData(ddsHTseq)$condition, 
                                      levels = c('Control','Experimental'))

# gut of DESeq2 analysis
dds <- DESeq(ddsHTseq)
res <- results(dds)
# order results by padj value (most significant to least)
res <- res[order(res$padj),]

head(res)
# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj 

# save data 'res' to csv!
write.csv(as.data.frame(res), file = 'DATE_DESeq2_initial_analysis.csv')

# send raw combined counts to tab delimited file for GSEA, etc.
write.table(as.data.frame(counts(dds)), file = 'DATE_DESeq2_raw_counts.txt', sep = '\t')

# send normalized counts to tab delimited file for GSEA, gene X to gene Y comparisions
# normalization is done by taking the median of the gene-wise fold changes with respect to pseudo-reference 
# sample of geometric means across all samples.
# calculates sizeFactors across each sequence library
write.table(as.data.frame(counts(dds),normalized=T), file = "DATE_DESeq2_normalized.txt", sep = '\t')

# can open in Excel, and other programs

# produce DataFrame of results of statistical tests
# could way to record experimental design 
mcols(res, use.names = T)
write.csv(as.data.frame(mcols(res, use.name = T)),file = "DATE-DESeq2-test-conditions.csv")

# replacing outlier value with estimated value as predicted  by distrubution using 
# "trimmed mean" approach. recommended if you have several replicates per treatment
# DESeq2 will automatically do this if you have 7 or more replicates
ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj < 0.1,
             cleaned = results(ddsClean)$padj < 0.1)
addmargins(tab)
write.csv(as.data.frame(tab),file = 'DATE-DESeq2-replaceoutliers.csv')
resClean <- results(ddsClean)
resClean <- resClean[order(resClean$padj),]
head(resClean)
write.csv(as.data.frame(resClean),file = 'DATE-DESeq2-replaceoutliers-results.csv')

###############################################################
# Exploritory data analysis of RNAseq data with DESeq2
#
# these next R scripts are for a variety of visualization, QC and other plots to 
# get a sense of what the RNAseq data looks like based on DESEq2 analysis 
#
# 1) MA plot
# 2) rlog stabilization and variance stabiliazation
# 3) variance stabilization plot 
# 4) heatmap of clustering analysis 
# 5) PCA plot 
#
#
################################################################

# MA plot of RNAseq data for entire dataset 
# http://en.wikipedia.org/wiki/MA_plot
plotMA(dds, ylim=c(-8,8),main = "RNAseq experiment")
dev.copy(png, "DATE-DESeq2_MAplot_initial_analysis.png")
dev.off()

# transform raw counts into normalized values 
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization <- VST is good for heatmaps
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

# save normalized values
write.table(as.data.frame(assay(rld),file='DATE-DESeq2-rlog-transformed-counts.txt', sep='\t')
write.table(as.data.frame(assay(vsd),file='DATE-DESeq2-vst-transformed-counts.txt', sep='\t')


# plot to show effect of transformation
# axis is square root of variance over the mean for all samples
par(mai = ifelse(1:4 <= 2, par('mai'),0))
px <- counts(dds)[,1] / sizeFactors(dds)[1]
ord <- order(px)
ord <- ord[px[ord] < 150]
ord <- ord[seq(1,length(ord),length=50)]
last <- ord[length(ord)]
vstcol <- c('blue','black')
matplot(px[ord], cbind(assay(vsd)[,1], log2(px))[ord, ],type='l', lty = 1, col=vstcol, xlab = 'n', ylab = 'f(n)')
legend('bottomright',legend=c(expression('variance stabilizing transformation'), expression(log[2](n/s[1]))), fill=vstcol)
dev.copy(png,"DATE-DESeq2_variance_stabilizing.png")
dev.off()

# clustering analysis
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(condition, type, sep=" : "))
# OR
# if you want the conditions used
# rownames(mat) <- colnames(mat) <- with(colData(dds),condition)
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(13, 13))
dev.copy(png, "DATE-DESeq2-clustering.png")
dev.off()

# Principal components plot
# will show additional clustering of samples
# showing basic PCA function in R from DESeq2 package 
# this lacks sample IDs and only broad sense of sample clustering
print(plotPCA(rld, intgroup = c("condition")))
dev.copy(png, "DATE-DESeq2_PCA_initial_analysis.png")
dev.off()

# or ggplot and rgl 2D and 3D PCA plots
# 3D PCA plot and 2D PCA

rv <- rowVars(assay(rld))
select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
pca <- prcomp(t(assay(rld)[select,]))

# set condition
condition <- c("condition1", 'condition2')

# using grouping to ID sub classes
scores$grouping <- c("group1", "group2")

scores <- data.frame(sampleFiles, pca$x, condition)
scores
(pca <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(grouping))))  # or factor(condition)
+ geom_point(size = 5)
+ ggtitle("Principal Components")
+ scale_colour_brewer(name = " ", palette = "Set1")
+ theme(
  plot.title = element_text(face = 'bold'),
  legend.position = 'bottom', 
  legend.key = element_rect(fill = 'NA'),
  legend.text = element_text(size = 10, face = "bold"),
  axis.text.y = element_text(colour = "Black"), 
  axis.text.x = element_text(colour = "Black"),
  axis.title.x = element_text(face = "bold"), 
  axis.title.y = element_text(face = 'bold'),
  panel.grid.major.x = element_blank(), 
  panel.grid.major.y = element_blank(),
  panel.grid.minor.x = element_blank(), 
  panel.grid.minor.y = element_blank(), 
  panel.background = element_rect(color = 'black',fill = NA)
))

save("2014-9-30-2D-PCA_ggplot2.RData")
ggsave(pcaplot, file="2014-9-30-2D-PCA_ggplot2.pdf")

# optional K-means based clustering
km <- kmeans(assay(rld),3) # not sure how to integrate into plot
pca$cluster <- as.factor(km$cluster)

# 3D plot with rgl 
plot3d(pca$x, size = 10, lit = T)
text3d(pca$x, text = sampleFiles, adj = -0.2, font = 2, col = colorRampPalette(c("blue",'red'))(2))
plot3d(pca$x, col = pca$cluster)
rgl.snapshot("DATE-3D-PCA.png",'png')

# scatter plot of rlog transformations between Sample conditions
# nice way to compare control and experimental samples 
head(assay(rld))
plot(log2(1+counts(dds,normalized=T)[,1:2]),col='black',pch=20,cex=0.3, main='Log2 transformed')
plot(assay(rld)[,1:2],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,3:4],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,5:6],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,7:8],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,9:10],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,11:12],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,13:14],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")
plot(assay(rld)[,15:16],col='#00000020',pch=20,cex=0.3, main = "rlog transformed")

#heatmap of data
# ordering by most variable rows (e.g. what makes PCA plot)
# best kind of map by ordering by most variable rows (e.g. what makes PCA plot)
# can add side colors to delineate different histologies, samples, etc.  
# with heatmap.2 function 
select <- order(rowVars(assay(vsdMF)), decreasing=T)[seq_len(min(1000,length(rv)))] 
my_palette <- colorRampPalette(c("blue", "white", "red"))(2048)
par(cex.main=0.8)
heatmap.2(assay(vsdMF)[select,], col=my_palette,
          scale="row", 
          key=T,
          symkey=T,
          density.info="density",
          densadj=0.3, 
          ColSideColors= c("darkgreen","darkgreen", 
                           "mediumpurple2","mediumpurple2",
                           "mediumpurple2","mediumpurple2",
                           "darkgreen","darkgreen",
                           "mediumpurple2","mediumpurple2",
                           "darkgreen","darkgreen",
                           "mediumpurple2","mediumpurple2",
                           "mediumpurple2","mediumpurple2"),
          trace="none",
          cexCol=0.6, 
          labRow=F, 
          margins=c(6,6),
          main="Transcriptome of \nEMT populations")
dev.copy(pdf,"DATE-DESeq2-heatmap2.pdf", width = X, height = Y)
dev.off()


# heatmap with heatmap3 function
# same row variance across dataset
select <- order(rowVars(assay(vsdMF)), decreasing=T)[seq_len(min(1000,length(rv)))] 
my_palette <- colorRampPalette(c("blue", "white", "red"))(1024)
heatmap3(assay(vsdMF)[select,], col=my_palette,
         labRow = F,
         cexCol = 0.8,
         ColSideColors= c("darkgreen","darkgreen", 
                           "mediumpurple2","mediumpurple2",
                           "mediumpurple2","mediumpurple2",
                           "darkgreen","darkgreen",
                           "mediumpurple2","mediumpurple2",
                           "darkgreen","darkgreen",
                           "mediumpurple2","mediumpurple2",
                           "mediumpurple2","mediumpurple2"),
         ColSideLabs = F,
         balanceColor=T)
dev.copy(pdf,"DATE-DESeq2-heatmap3.pdf", width = X, height = Y)
dev.off()

# how to specify individual genes from heatmap
sample <- rowMeans(counts(dds))
select[1:50] # will give you vector of numbers
sample[select[1:50]] # will give you gene name and normalized value for top 50 genes


sessionInfo()


###############################################################
#
# Optional analyses 
# 
###############################################################

# multifactor designs
# can analysis with more than one factor influencing the counts  (e.g., different library preps, PE vs. SE)
# from manual section 1.5 for more information and rationale

# make copy of DESeq2 results
ddsMF <- dds

# change design formulate controlling for Batch
design(ddsMF) <- formula(~ Batch + condition)

# rerun DESeq analysis
ddsMF <- DESeq(ddsMF)
resMF <- results(ddsMF)

# order by padj values
resMF <- resMF[order(resMF$padj),]

head(resMF)
# should see DataFrame of baseMean, log2Foldchange, stat, pval, padj 
# save data 'resMF' to csv!
write.csv(as.data.frame(resMF),file='DATE-DESeq2_batchcontroll_initial_analysis.csv')

# visualize data wihtout independent filtering
resNoFilt <- results(dds, independentFiltering=F)
filter_table <- table(filtering=(res$padj <.1), noFiltering=(resNoFilt$padj < .1))
resNoFilt <- resNoFilt[order(resNoFilt$padj),]
head(resNoFilt)
write.table(as.data.frame(filter_table),file='filter_table')
write.csv(as.data.frame(resNoFilt),file='DATE-DESeq2_results_NOFILTERING.CSV')

# visualize pvales discarded by filtering
use <- res$baseMean > attr(res,"filterThreshold")
table(use)
h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=F)
h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=F)
colori <- c('do not pass'="khaki", 'pass'="powderblue")
barplot(height = rbind(h1$counts, h2$counts), beside = F,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))

plot(seq(along=which(showInPlot)), resFilt$pvalue[orderInPlot][showInPlot],
     pch=".", xlab = expression(rank(p[i])), ylab=expression(p[i]))
abline(a=0, b=alpha/length(resFilt$pvalue), col="red3", lwd=2)

# filtering threshold
attr(res, 'filterThreshold')
plot(attr(res,'filterNumRej'), type = 'b', ylab='number of rejection')
dev.copy(png,'DATE-DESeq2_filtering_threshold.png')
dev.off()

# plot of variance stabilization and rlog transformation versus untransformed data
par(mfrow=c(1,3))
notAllZero <- (rowSums(counts(dds))>0)
meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1), ylim = c(0,4))
meanSdPlot(assay(rld[notAllZero,]), ylim = c(0,4))
meanSdPlot(assay(vsd[notAllZero,]), ylim = c(0,4))