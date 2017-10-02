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
setwd("/Users/bachum/Desktop/OneDrive/NIH_Mac/Oda_Whsc/")
basedir <- "/Users/bachum/Desktop/OneDrive/NIH_Mac/Oda_Whsc/" #setwd("/Users/bachum/Desktop/OneDrive/NIH_Mac/Deseq2_Analysis/")
basedir
Countdirectory <- "/Users/bachum/Desktop/OneDrive/NIH_Mac/Oda_Whsc/Htseq-Counts/Zmynd11/"
#Build a design table of experimental parameters and save as a .txt file
metadata <- read.table("Targets_Zmynd11.txt", header = TRUE)
#Read in the Sample Count files using grep command
sampleFiles <- grep("txt",list.files(Countdirectory),value=TRUE)
# load metadata
metadata$sampleFiles <- paste(metadata$files)
metadata
sampleTable <- data.frame(sampleName = metadata$samples, fileName= metadata$sampleFiles, 
                          genotype = metadata$genotype, 
                          Day = metadata$Treatment)
sampleTable
sampleTable$time <- as.factor(sampleTable$time)
sampleTable$Day <- as.factor(sampleTable$Day)
sampleTable$genotype <- as.factor(sampleTable$genotype)
sampleTable
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = Countdirectory,design = ~ Day)
dds
noint <- rownames(dds) %in% c("N_unmapped", "N_multimapping","N_noFeature","N_ambiguous")
head(noint)
tail(noint)
#Filtering
keep <- rowSums(cpm(d) > 1) >= 3 & !noint
d
d <- d[keep,]

# review the created object
dds
# column information
colData(dds)
colnames(dds)
# relevel to get 0 h as a reference
dds$Day <- relevel(dds$Day, "NT")
dds$genotype <- relevel(dds$genotype, "Wt")
dds$time <- relevel(dds$time, "0")
dds$Treatment <- relevel(dds$Treatment, "NT")


se <- DESeqDataSet(dds, design = ~ genotype + Day)
countdata <- assay(se)
head(countdata, 3)
coldata <- colData(se)
(ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                  colData = coldata,
                                  design = ~ genotype + Day))
nrow(dds)

## [1] 64102

dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)
sampleDists <- dist( t( assay(rld) ) )
sampleDists
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$genotype, rld$Day, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(rld, intgroup = c("genotype", "Day"))
(pcaData <- plotPCA(rld, intgroup = c( "genotype", "Day"), returnData=TRUE))
percentVar <- round(100 * attr(pcaData, "percentVar"))
library("ggplot2")

ggplot(pcaData, aes(PC1, PC2, color=Day, shape=genotype)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()


dds <- DESeq(dds)
(res <- results(dds))
results(dds, contrast=c("Day", "Treat", "NT"))

colnames(res)
resOrdered <- res[order(res$padj),]
summary(res)
resOrderedDF<- res[ order(res$log2FoldChange),]
resOrderedDF <- resOrderedDF[!is.na(resOrderedDF$pvalue),]
resSig <- subset(resOrderedDF, padj < 0.01 & abs(resOrderedDF$log2FoldChange)>=1)
dim(resSig)
summary(resSig)
resSig <- subset(resOrderedDF, padj < 0.1 & abs(resOrderedDF$log2FoldChange)>=1)
dim(resSig)
summary(resSig)
norm.counts <- counts(dds, normalized=TRUE)
log.norm.counts <- log2(norm.counts + 1)
resdifferential <- merge(as.data.frame(resOrderedDF), as.data.frame(log.norm.counts), by="row.names", sort=FALSE)
colnames(resdifferential)
names(resdifferential)[1] <- "Gene"
head(resdifferential)
write.csv(as.data.frame(resdifferential), file="KO.IFN2h.csv")
