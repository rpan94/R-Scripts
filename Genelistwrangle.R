library(ggplot2)
library(gridExtra)
library(tidyverse)
mtc <- read_tsv(file = "C:/Users/panrl/Documents/RNASeq/1-9-2018/KOP3vsWtP2.complete.txt", col_names = T)
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
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("org.Mm.eg.db")
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
dim(mtc)
mtc <- dplyr::select(mtc,symbol, WtP2, KOP2, KOP3, FC, log2FoldChange,pvalue,padj,entrez,Gene_Id)
write.csv(mtc, file = "C:/Users/panrl/Documents/RNASeq/1-9-2018/KOP3vsWtP2.complete.csv", row.names = T)

damup <- read_csv("C:/Users/panrl/Documents/RNASeq/DAM_Up.csv")
glimpse(damup)
dim(damup)
head(damup)
damup <- dplyr::filter(damup, change == 1)
dim(damup)

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
mtc1 <- mtc %>%
  filter(symbol %in% damup$ID)
dim(mtc1)
write.csv(mtc1, file = "C:/Users/panrl/Documents/RNASeq/1-9-2018/KOP3vsWtP2.DAM.csv", row.names = T)

##################################################################################################
#check row for row if a combination exists in another dataframe and add annotation saying yes if found and no if not
mtc$Dam.Gene <- ifelse(is.na(match(paste0(mtc$symbol), 
                                   paste0(damup$ID))),"No", "Yes")
filter(mtc, Dam.Gene == "Yes")
colnames(mtc)
p <- ggplot(mtc, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(col=Dam.Gene), size =1, alpha = 1, color = 'grey') + 
  geom_point(data = subset(mtc, Dam.Gene == 'Yes'), color = 'red')
p
#ggplot(df) +
#geom_point(aes(x = x, y = y, color = label,  size = size)) +
#  geom_point(data = subset(df, label == 'point'),
            # aes(x = x, y = y, color = label, size = size))
library(ggrepel)
#Select your gene of interest using the row number of excel sheet
r <- p+geom_text_repel(data=dplyr::filter(mtc, symbol %in% c("Ctsz",
                                                             "Ctss",
                                                             "Apoe",
                                                             "Apoc1",
                                                             "Apoc4",
                                                             "Npc2",
                                                             "Ch25h",
                                                             "Lpl",
                                                             "Trem2",
                                                             "Axl",
                                                             "Tyrobp",
                                                             "Igf1",
                                                             "Spp1",
                                                             "Gpnmb",
                                                             "Itgax",
                                                             "Csf1",
                                                             "Clec7a",
                                                             "Lyz2",
                                                             "Lgals3bp",
                                                             "Hifa",
                                                             "Plin2",
                                                             "Dpp7",
                                                             "Hexa",
                                                             "Ank",
                                                             "Acaca",
                                                             "Soat1")), aes(label= symbol))
r + geom_vline(xintercept = c(-1,1)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

###############################################################################################################################
