#source("http://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("BSgenome")
library("DESeq2")
library("ggplot2")
library("vsn")
library("RColorBrewer")
library("gplots")
library("pheatmap")
library("genefilter")
library("Category")
library(GOstats)
library(org.Mm.eg.db)
#Reading the data=============> Step 1
################################################################################################
getwd()
setwd("/Users/bachum/Dropbox/Myco-IRF8-Htseq-Counts/")
basedir <- "/Users/bachum/Dropbox/Myco-IRF8-Htseq-Counts/" #setwd("/Users/bachum/Desktop/OneDrive/NIH_Mac/Deseq2_Analysis/")
basedir
Countdirectory <- "/Users/bachum/Dropbox/Myco-IRF8-Htseq-Counts/Ht-seq-Counts-T96/"
#Build a design table of experimental parameters and save as a .txt file
metadata <- read.table("TargetsT96Desq.txt", header = TRUE)
#Read in the Sample Count files using grep command
#sampleFiles <- grep("tab",list.files(Countdirectory),value=TRUE)
# load metadata
#metadata$sampleFiles <- paste(metadata$files)
metadata
sampleTable <- data.frame(sampleName = metadata$samples, fileName= metadata$files, 
                          Genotype = metadata$Genotype,
                          Treatmemt = metadata$Treatmemt, group= metadata$group)
sampleTable
sampleTable$Days <- as.factor(sampleTable$Days)
sampleTable$Treatmemt <- as.factor(sampleTable$Treatmemt)
sampleTable$Genotype <- as.factor(sampleTable$Genotype)
sampleTable
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, 
                                  directory = Countdirectory,
                                  design = ~ Genotype + Treatmemt + Genotype:Treatmemt)
# review the created object
dds
# column information
colData(dds)
colnames(dds)
# relevel to get 0 h as a reference
dds$Treatment <- relevel(dds$Treatment, "NT")
dds$Genotype <- relevel(dds$Genotype, "hIRF8")
dds$Days <- relevel(dds$Days, "0")


ddsTC <- DESeqDataSet(dds, ~ genotype + time + genotype:time)
ddsTC <- DESeq(dds, test="LRT", reduced = ~ genotype + time)
##################################################################################
#Pre-filtering the dataset
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]

nrow(dds)

lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)
library("vsn")
meanSdPlot(cts, ranks = FALSE)



#The rlog transformation

rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)


par( mfrow = c( 1, 2 ) )
dds <- estimateSizeFactors(dds)
plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1),
     pch=16, cex=0.3)
plot(assay(rld)[,1:2],
     pch=16, cex=0.3)


#Sample distances
sampleDists <- dist(t(assay(rld)))
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
library("pheatmap")
library("RColorBrewer")
rownames(sampleDistMatrix) <- paste( rld$Genotype, rld$Treatment, sep="." )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#PCA plot
plotPCA(rld, intgroup = c("Genotype", "Treatmemt"))
colData(dds)

#PCA plot using the rlog-transformed values
(data <- plotPCA(rld, intgroup = c( "Genotype", "Treatmemt"), returnData=TRUE))

percentVar <- round(100 * attr(data, "percentVar"))
library("ggplot2")

s <- ggplot(data, aes(PC1, PC2, color=Genotype, shape=Treatmemt)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 
s


library("ggthemes")

s+ theme_bw() 


dds <- DESeq(dds)
dds
res <- results(dds)
colnames(res)
deseq2.fc <- res$log2FoldChange
names(deseq2.fc)=rownames(res)
exp.fc=deseq2.fc
exp.fc
out.suffix="deseq2"
require(gage)
data(kegg.gs)
fc.kegg.p <- gage(exp.fc, gsets = kegg.gs, ref = NULL, samp = NULL)
sel <- fc.kegg.p$greater[, "q.val"] < 0.1 & !is.na(fc.kegg.p$greater[, "q.val"])
path.ids <- rownames(fc.kegg.p$greater)[sel]
sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 & !is.na(fc.kegg.p$less[,"q.val"])
path.ids.l <- rownames(fc.kegg.p$less)[sel.l]
path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
require(pathview)
pv.out.list <- sapply(path.ids2[1:3], 
                      function(pid) pathview(gene.data =  exp.fc, 
                                             pathway.id = pid,species = "hsa", 
                                             out.suffix=out.suffix))



resultsNames(ddsTC)
mcols(resTC, use.names=TRUE)

resTC <- results(ddsTC, name="genotype_G388S_vs_Wt", test="Wald")
resTC[which.min(resTC$padj),]
summary(resTC)
resOrderedDF<- resTC[ order(resTC$log2FoldChange),]
resOrderedDF <- resOrderedDF[!is.na(resOrderedDF$pvalue),]
resSig <- subset(resOrderedDF, padj < 0.01 & abs(resOrderedDF$log2FoldChange)>=1)
dim(resSig)
summary(resSig)
norm.counts <- counts(ddsTC, normalized=TRUE)
log.norm.counts <- log2(norm.counts + 1)
resdifferential <- merge(as.data.frame(resSig), as.data.frame(log.norm.counts), by="row.names", sort=FALSE)
colnames(resdifferential)
names(resdifferential)[1] <- "Gene"
head(resdifferential)
write.csv(as.data.frame(resdifferential), file="genotype_G388S_vs_Wt.csv")
