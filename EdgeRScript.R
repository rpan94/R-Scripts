library(RColorBrewer)
library(genefilter)
library(edgeR)
library(ggplot2)
library(knitr)
getwd()
setwd("/Users/bachum/Desktop/OneDrive/Edger/Count_ucscmm9_STAR/")
#Reading the data
targets <- readTargets("/Users/bachum/Desktop/OneDrive/Edger/Targets.txt")
targets
targets$time <-factor(targets$time)
targets$genotype <-factor(targets$genotype)
targets$Treatment <- factor(targets$Treatment)
targets
head(targets)
d <- readDGE(targets, skip=5, comment.char = "!")
noint <- rownames(d) %in% c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
tail(noint)
d$samples
dim(d)
d
#Filtering
keep <- rowSums(cpm(d) > 1) >= 3 & !noint
d
d <- d[keep,]
dim(d)
d[10361,]
colnames(d)
###EdgeR
y <- DGEList(d)
y
genotype <- c("Wt", "Wt","Wt", "Wt","Wt", "Wt","Wt","Wt", "Wt","Wt", "Wt","KO","KO","KO","KO","KO","KO","KO","KO","KO","KO")
genotype
Treatment <- c("Ctrl", "IFN","IFN","IFN","IFN","Ctrl", "IFN","IFN","IFN","IFN","Ctrl", "IFN","IFN","IFN","IFN","Ctrl", "IFN","IFN","IFN","IFN")
time <- c("0", "6", "12", "24", "24","0", "6", "12", "48", "48","0", "6", "12", "24", "48","0", "6", "12", "24", "48")
#time <- c("0", "1", "3", "6", "24","0", "1", "3", "6", "24","0", "1", "3", "6", "24","0", "1", "3", "6", "24")

grouping <- factor(paste(genotype, time, sep="."))
grouping
design <- model.matrix(~0 + grouping)
colnames(design) <- levels(grouping)
colnames(design)

#Annotation
#library(org.Mm.eg.db)
#ls("package:org.Mm.eg.db")
#columns(org.Mm.eg.db)
#class(y)
#colnames(y)
#row.names(y)
#idfound <- y$%in% mappedRkeys(org.Mm.egENSEMBL)
#idfound
#y <- y[idfound,]
dim(y)
RPKM <- rpkm(y)
### plot an MDS plot to see how all the samples scale in multidemensional space
plotMDS.DGEList( d , main = "MDS Plot for Count Data", labels = colnames( d$counts ) )
###Create Groups represententing individual sample properties so they can be contrasted
#Group <- factor(paste(targets$genotype,targets$time,sep="."))
targets$time <- as.factor(targets$time)
targets$genotype <- as.factor(targets$genotype)
#targets$Memory <- as.factor(targets$Memory)

targets$time <- relevel(targets$time, ref="0")
#targets$Memory <- relevel(targets$Memory, ref="Nai")
targets$genotype <- relevel(targets$genotype, ref="Wt")
targets$Treatment <- relevel(targets$Treatment, ref="Ctrl")

#Nested interaction formulas
#Then form the design matrix:
design <- model.matrix(~0 + grouping)
colnames(design) <- levels(grouping)
design
colnames(design)
#If you want to compare the KO to Wt genotype, you can do this by comparing the KO-control group to the Wt-control group, 
#or by comparing the KO-treatment group to KO-treatment group. 
#This can be done by running:
KOvWtcontrol <- makeContrasts(KO.Ctrl- Wt.Ctrl, levels=design)
KOvWtcontrol
KOvWtTreatment <- makeContrasts(KO.IFN - Wt.IFN, levels=design)
KOvWtTreatment

Wt_IFN <- makeContrasts(Wt.IFN-Wt.Ctrl, levels = design)
KO_IFN <- makeContrasts(KO.IFN-KO.Ctrl, levels = design)

#Interaction at any time
design <- model.matrix(~genotype,data = targets)
colnames(design)


###Calculate Dispersions
y <- calcNormFactors(y,method="TMM")
y <- estimateCommonDisp(y)
#design <- model.matrix(~0+Group)
#colnames(design) <- levels(Group)
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
colnames(fit)
#Create contrasts
colnames(design)
my.contrasts <- makeContrasts(NT=KO.0-Wt.0,
                              Six=KO.6-KO.0,
                              Twelve=KO.12-KO.0,
                              TwentyFour=KO.24-KO.0,
                              FortyEight=KO.48-KO.0,
                              levels=design)
colnames(my.contrasts)

fit <- glmFit(y,design)
# Run the likelihood ratio test on all contrasts
KO_6h <- glmLRT(fit, contrast=my.contrasts[,"Six"])
KO_12h <- glmLRT(fit, contrast=my.contrasts[,"Twelve"])
KO_24h <- glmLRT(fit, contrast=my.contrasts[,"TwentyFour"])
KO_48h <- glmLRT(fit, contrast=my.contrasts[,"FortyEight"])

#Now most of the above contrasts are directly available as coecients:
#is the baseline KO vs Wt comparison
lrt_NT <- glmLRT(fit, coef=2)
lrt_6h <- glmLRT(fit, coef=3)
lrt_12h <- glmLRT(fit, coef=4)
lrt_24h <- glmLRT(fit, coef=5)
lrt_48h <- glmLRT(fit, coef=6)
#lrt_All <- glmLRT(fit, coef=c(3:6))
#lrt_All

l.control <- glmLRT(fit,contrast=KOvWtcontrol)
l.treatment <- glmLRT(fit,contrast=KOvWtTreatment)
Wt.IFN <- glmLRT(fit,contrast=Wt_IFN) 
KO.IFN <- glmLRT(fit,contrast=KO_IFN)
## Make a table of results
etable <- topTags(lrt_NT, n=nrow(keep))$table
colnames(etable)
etable <- etable[order(etable$logFC), ]

sigGenes <- etable[abs(etable$logFC)>=0.58 & etable$FDR < 0.01, ]
summary(sigGenes)
tail(sigGenes)
dim(sigGenes)
head(sigGenes)
write.csv(sigGenes, file = "KOvsWt2replicate_Refseq_0h_1.5_Fold_0.01_FDR.csv")

#row.names(sigGenes)
#Add Entrz_IDs to the significant genes
library("AnnotationDbi")
library("org.Mm.eg.db")
columns(org.Mm.eg.db)
sigGenes$RESEQ <- mapIds(org.Mm.eg.db,
                           keys=row.names(sigGenes),
                           column="REFSEQ",
                           keytype="ENTREZID",
                           multiVals="first")
#Add Ens_IDs to the significant genes
sigGenes$Ens_ID <- mapIds(org.Mm.eg.db,
                            keys=row.names(sigGenes),
                            column="ENSEMBL",
                            keytype="SYMBOL",
                            multiVals="first")

#Add MGI to the significant genes
sigGenes$Symbol <- mapIds(org.Mm.eg.db,
                          keys=row.names(sigGenes),
                          column="SYMBOL",
                          keytype="ENTREZID",
                          multiVals="first")

#Add MGI to the significant genes
sigGenes$GENENAME <- mapIds(org.Mm.eg.db,
                          keys=row.names(sigGenes),
                          column="GENENAME",
                          keytype="ENTREZID",
                          multiVals="first")
etable$Symbol <- mapIds(org.Mm.eg.db,
                          keys=row.names(etable),
                          column="SYMBOL",
                          keytype="ENTREZID",
                          multiVals="first")
colnames(sigGenes)
dim(etable)
dim(sigGenes)
count_data <- cpm(d, prior.count=2, log=TRUE)
colnames(count_data)
Results <- merge(as.data.frame(etable), as.data.frame(count_data), by="row.names", sort=FALSE)
colnames(Results)
#row.names(Results)
#names(Results)[1] <- "Gene"
head(Results)
write.csv(as.data.frame(Results), file="KOvsWt_Basal_All_FDR_1replicate.csv")

dev.off()
#For Heatmap

#Open the csv file remoce the first numerical column
#Prepare expression data
#Assume you have already installed heatmap3 package. And the 3 example file were download in current working directory. You can use the following codes in R to generate the figures.
#Prepare expression data
counts<-read.csv("counts.csv",header=T,row.names=1)
head(counts)
colnames(counts)
#Prepare column side annotation
clinic<-read.csv("clinic.csv",header=T,row.names=1)
clinic
colnames(clinic)
#Prepare row side color bar annotation
edgeR_result<-read.csv("edgeR_results.csv",header=T,row.names=1)
colnames(edgeR_result)
temp1<-(edgeR_result$logFC)
temp2<--log10(edgeR_result$FDR)
temp1<-colByValue(as.matrix(temp1),range=c(-4,4),col=colorRampPalette(c('chartreuse4','white','firebrick'))(1024))
temp2<-colByValue(as.matrix(temp2),range=c(0,5),col=heat.colors(1024))
colGene<-cbind(temp1,temp2)
row.names(colGene)<-row.names(edgeR_result)
colnames(colGene)<-c("log2FC","-Log10P")

#Generate Figure1
#counts, colGene and clinic were read throught the csv file
##Assume counts has counts information for each gene, colGene has the colors for each gene, clinic has the clinic information for each sample
temp<-apply(counts,1,sd)
temp
selectedGenes<-rev(order(temp))[1:242]
selectedGenes
heatmap3(counts[selectedGenes,],labRow="",margin=c(4,20),RowSideColors=colGene[selectedGenes,],ColSideCut=0.85,ColSideAnn=clinic,ColSideFun=function(x) showAnn(x),ColSideWidth=0.6,balanceColor=T)
dev.off()



#Gene ontology analysis
go <- goana(sigGenes)
topGO(go, ont="BP")
keg <- kegga(sigGenes)
topKEGG(keg,sort = NULL, number = 20L, truncate.path = NULL)
