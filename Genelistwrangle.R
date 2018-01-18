library(ggplot2)
library(gridExtra)
library(tidyverse)
mtc <- read_tsv(file = "/Users/bachum/Desktop/OneDrive/Pan-Richard/tables/KOP2vsWtP2.complete.txt", col_names = T)
colnames(mtc)
glimpse(mtc)
dim(mtc)
library(tidyr)
#remove all NA rows
mtc <- mtc %>% drop_na()
glimpse(mtc)
dim(mtc)
colnames(mtc)
#mtc <- read.csv(file="Scatter_Plot_NT.csv", header = T)
tail(mtc)
#Select only the required columns
mtc <- dplyr::select(mtc,Id, WtP2, KOP2, KOP3, FC, log2FoldChange,pvalue,padj)
head(mtc)
dim(mtc)
#Remove dot from ensemble annotation using command and add an additional variable Gene_ID
mtc <- mutate(mtc, Gene_Id = gsub("\\..*","",mtc$Id))
mtc <- dplyr::select(mtc,Gene_Id, WtP2, KOP2, KOP3, FC, log2FoldChange,pvalue,padj)
dim(mtc)
head(mtc)
min(mtc$WtP2)
max(mtc$WtP2)
library("AnnotationDbi")
library("org.Mm.eg.db")
columns(org.Mm.eg.db)
mtc$symbol <- mapIds(org.Mm.eg.db,
                     keys=mtc$Gene_Id,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
mtc$entrez <- mapIds(org.Mm.eg.db,
                     keys=mtc$Gene_Id,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

head(mtc)
mtc <- dplyr::select(mtc,symbol, WtP2, KOP2, KOP3, FC, log2FoldChange,pvalue,padj,entrez,Gene_Id)
write.csv(mtc, file = "/Users/bachum/Desktop/OneDrive/Pan-Richard/tables/KOP2vsWtP2.complete.csv", row.names = T)

damup <- read_csv("/Users/bachum/Desktop/OneDrive/Pan-Richard/tables/Dam-UP.csv")
glimpse(damup)
head(damup)
class(damup)

#################################################################################################
#How to subset a column in data frame based on another data frame/list
#A better option would be data.table
library(data.table)
setDT(mtc)[symbol %chin% damup$ID]
#or
subset(mtc, symbol %in% damup$ID)
#or
mtc %>%
  filter(symbol %in% damup$ID)
#################################################################################################
