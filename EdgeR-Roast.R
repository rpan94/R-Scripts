library(edgeR)
library(limma)
source("https://bioconductor.org/biocLite.R")
library('EDASeq')
library('Glimma')
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)

# Read the data into R
seqdata <- read.delim("/Users/bachum/Desktop/OneDrive/Edger/RsubreadCounts.txt", stringsAsFactors = FALSE)
# Read the sample information into R
targets <- read.delim("/Users/bachum/Desktop/OneDrive/Edger/Samples.txt")
Genotype <- factor(targets$Genotype)
design <- model.matrix(~Genotype)
targets

head(seqdata)
dim(seqdata)
# Remove first 2 columns from seqdata
countdata <- seqdata[,-(1:2)]
# Look at the output
head(countdata)
# Store EntrezGeneID as rownames
rownames(countdata) <- seqdata[,1]
head(countdata)
colnames(countdata)
# using substr, you extract the characters starting at position 1 and stopping at position 7 of the colnames
#colnames(countdata) <- substr(colnames(countdata),start=1,stop=10)
head(countdata)
table(colnames(countdata)==targets$FileName)
# Obtain CPMs
myCPM <- cpm(countdata)
# Have a look at the output
head(myCPM)
# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 1
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)
# Summary of how many TRUEs there are in each row
# There are 11433 genes that have TRUEs in all 12 samples.
table(rowSums(thresh))
# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- countdata[keep,]
summary(keep)
# Let's have a look and see whether our threshold of 0.5 does indeed correspond to a count of about 10-15
# We will look at the first sample
plot(myCPM[,1],countdata[,1])
# Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(myCPM[,1],countdata[,1],ylim=c(0,50),xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.5)
y <- DGEList(counts.keep)
# have a look at y
y
# See what slots are stored in y
names(y)
# Library size information is stored in the samples slot
y$samples
y$samples$lib.size
# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
barplot(y$samples$lib.size,names=colnames(y),las=3)
# Add a title to the plot
title("Barplot of library sizes")
# Get log2 counts per million
logcounts <- cpm(y,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")
# You need to give plotRLE a matrix of counts
plotRLE(y$counts,las=2)
title("RLE plot of unnormalised counts")
plotMDS(y)
# We specify the option to let us plot two plots side-by-sde
par(mfrow=c(1,2))
# Let's set up colour schemes for CellType
# How many cell types and in what order are they stored?
levels(targets$Genotype)
## Let's choose purple for basal and orange for luminal
col.cell <- c("purple","orange")[targets$Genotype]
col.cell
# Redo the MDS with cell type colouring
plotMDS(y,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend("bottomright",fill=c("purple","orange"),legend=levels(targets$Genotype))


# Let's combine status and cell type info on one plot using plotting characters and colours
# Make a new group variable that joins together cell type and status with the paste command
par(mfrow=c(1,1))
group <- factor(paste(targets$genotype,targets$condition,sep="."))
group
# the table function will tell us how many samples we have in each group
table(group)
points <- c(0,1,2,15,16,17)
colors <- rep(c("purple","orange","blue","red","black"), 2)
plotMDS(y, col=colors[group], pch=points[group])
legend("bottomright", legend=levels(group), pch=points, col=colors, ncol=2,cex=0.7)
plotMDS(y,dim=c(3,4),col=colors[group], pch=points[group])
legend("topright", legend=levels(group), pch=points, col=colors, ncol=2,cex=0.8)

pca <- prcomp(logcounts)
plot(pca$sdev)
labels <- paste(targets$name, targets$condition, targets$Treatment)
glMDSPlot(y, labels=labels, groups=group, folder="mds")
# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)
# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)
## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Plot the heatmap
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")
# Save the heatmap
png(file="High_var_genes.heatmap.png")
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")
dev.off()
# Apply normalisation to DGEList object
y <- calcNormFactors(y)
y$samples
par(mfrow=c(1,2))
plotMD(logcounts,column = )
abline(h=0,col="grey")
plotMD(logcounts,column = 4)
abline(h=0,col="grey")
par(mfrow=c(1,2))
plotMD(y,column = 4)
abline(h=0,col="grey")
plotMD(y,column = 4)
abline(h=0,col="grey")
# Look at group variable again
group
# Specify a design matrix without an intercept term
design <- model.matrix(~ 0 + group)
design
## Make the column names of the design matrix a bit nicer
colnames(design) <- levels(group)
design
par(mfrow=c(1,1))
v <- voom(y,design,plot = TRUE)
v

# What is contained in this object?
names(v)
par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(v$E),col="blue")
# Fit the linear model
fit <- lmFit(v)
colnames(design)
cont.matrix <- makeContrasts(Six = KO.0h - Wt.0h, levels=design)
cont.matrix
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
fit.cont
dim(fit.cont)
summa.fit <- decideTests(fit)
summary(summa.fit)
topTable(fit.cont,coef="Six",sort.by="p")
## This will give the same output
topTable(fit.cont,coef=1,sort.by="p")
columns(org.Mm.eg.db)


ann1 <- mapIds(org.Mm.eg.db,
       keys=row.names(fit.cont),
       column="ENTREZID",
       keytype="ENSEMBL",
       multiVals="first")

ann1 <- mapIds(org.Mm.eg.db,
              keys=row.names(fit.cont),
              column="GENENAME",
              keytype="ENSEMBL",
              multiVals="first")
ann1 <- mapIds(org.Mm.eg.db,
              keys=row.names(fit.cont),
              column="SYMBOL",
              keytype="ENSEMBL",
              multiVals="first")
str(ann1)
ann <- select(org.Mm.eg.db,keys=rownames(fit.cont),columns=c("ENTREZID","SYMBOL","GENENAME"))


limma.res <- topTable(fit.cont,coef=1,sort.by="p",n="Inf")
# We want to highlight the significant genes. We can get this from decideTests.
par(mfrow=c(1,2))
plotMD(fit.cont,coef=1,status=summa.fit[,"Six"])

# For the volcano plot we have to specify how many of the top genes to hightlight.
# We can also specify that we want to plot the gene symbol for the highlighted genes.
# let's highlight the top 100 most DE genes
volcanoplot(fit.cont,coef=1,highlight=100,names=fit.cont$genes$SYMBOL)
# Let's look at the first gene in the topTable, Wif1, which has a rowname 24117
?stripchart
stripchart(v$E["10",]~group,vertical=TRUE,las=2,cex.axis=0.8,pch=16,col=1:6,method="jitter")
# Let's decide that we are only interested in genes that have a absolute logFC of 1.
# This corresponds to a fold change of 2, or 0.5 (i.e. double or half).
# We can perform a treat analysis which ranks our genes according to p-value AND logFC.
# This is easy to do after our analysis, we just give the treat function the fit.cont object and specify our cut-off.
fit.treat <- treat(fit.cont,lfc=1)
res.treat <- decideTests(fit.treat)
summary(res.treat)
topTable(fit.treat,coef=1,sort.by="p")
# Notice that much fewer genes are highlighted in the MAplot
plotMD(fit.treat,coef=1,status=res.treat[,"Six"])
abline(h=0,col="grey")


limma.res.pval <- topTable(fit.cont,coef="Six",n="Inf",p.val=0.01)
head(limma.res.pval)
dim(limma.res.pval)
o <- order(limma.res.pval$logFC,decreasing = TRUE)
lfc.ordered <- limma.res.pval[o, ]
head(lfc.ordered)


up50 <- head(lfc.ordered, 50)
dim(up50)

mu50 <- match(rownames(up50), rownames(v))
up50.counts.voom <- v$E[mu50, ]
rownames(up50.counts.voom) <- up50$SYMBOL
heatmap.2(up50.counts.voom, Rowv = F, Colv = F, trace="none", ColSideColors=nice.col[group],col=bluered(100),labCol=group,margins=c(8,7))


go <- goana(fit.cont, coef=1,species = "Mm")
topGO(go, n=10)
topGO(go, ont="CC", sort="Up", n=10)


# Load in the mouse c2 gene sets
# The R object is called Mm.c2
load("mouse_c2_v5.rdata")
# Have a look at the first few gene sets
names(Mm.c2)[1:5]
# Number of gene sets in C2
length(Mm.c2)
c2.ind <- ids2indices(Mm.c2, rownames(Y))
gst.camera <- camera(v,index=c2.ind,design=design,contrast = cont.matrix[,1],inter.gene.cor=0.05)
gst.camera[1:10,]
table(gst.camera$FDR < 0.05)


topTable(fit.cont,coef=1,sort.by="p")
grep("IMM",names(c2.ind))
# Let's save these so that we can subset c2.ind to test all gene sets with WNT in the name
IMM <- grep("IMM",names(c2.ind))
IMM
# What are these pathways called?
names(c2.ind)[IFN]

#Let’s use ROAST to see if these WNT related gene sets tend to be differentially expressed. Note that the syntax for camera and roast is almost identical.
wnt.rst <- roast(v,index=c2.ind[IMM],design=design,contrast=cont.matrix[,1],nrot=9999)
wnt.rst

# Have a look at the logFCs and t-statistics in fit.cont
names(fit.cont)
head(fit.cont$coefficients)
head(fit.cont$t)
par(mfrow=c(2,1))
# barcode plot with logFCs
barcodeplot(fit.cont$coeff[,1], index=c2.ind[["REACTOME_INTERFERON_ALPHA_BETA_SIGNALING"]], main="REACTOME_INTERFERON_SIGNALING")
# barcode plot using t-statistics
barcodeplot(fit.cont$t[,1], index=c2.ind[["REACTOME_INTERFERON_GAMMA_SIGNALING"]], main="REACTOME_INTERFERON_GAMMA_SIGNALING")





library(GO.db)
help(GO.db)
ls("package:GO.db")
cyt.go <- c("GO:0032465", "GO:0000281", "GO:0000920")
cyt.go
Rkeys(org.Mm.egGO2ALLEGS) <- cyt.go
ind <- ids2indices(as.list(org.Mm.egGO2ALLEGS), fit.cont$genes$ENTREZID)
library(statmod)
fr <- fry(fit.cont, index=ind, design=design)
