res <- read.csv(file = "/Users/bachum/Desktop/KO.4vsWt.4.complete.csv",
                  header = T, check.names = T,blank.lines.skip = T)
res <- as.data.frame(res)
res<-na.omit(res)
colnames(res)
library("AnnotationDbi")
library("org.Mm.eg.db")
columns(org.Mm.eg.db)
colnames(res)
rownames(res) <- res[, 1] ## set rownames

row.names(res)
library("AnnotationDbi")
library("org.Mm.eg.db")
columns(org.Mm.eg.db)
res$entrez = mapIds(org.Mm.eg.db,
                    keys=row.names(res), 
                    column="ENTREZID",
                    keytype="SYMBOL",
                    multiVals="first")
colnames(res)
head(res, 50)
res<-na.omit(res)
head(res, 50)
library(pathview)
library(gage)
library(gageData)
data(kegg.sets.mm)
data(sigmet.idx.mm)
kegg.sets.mm = kegg.sets.mm[sigmet.idx.mm]
head(kegg.sets.mm, 3)
#The gage() function requires a named vector of fold changes, 
#where the names of the values are the Entrez gene IDs. 
foldchanges <- res$log2FoldChange
head(foldchanges)
names(foldchanges) = res$entrez
head(foldchanges)
# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.mm, same.dir=TRUE)

# Look at both up (greater), down (less), and statatistics.
lapply(keggres, head)
keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater)
head(keggrespathways,50)
# Get the IDs.
keggresids <- keggrespathways[,1]
keggresids <- keggresids[1:2]
keggresids
getwd()
pv.out<- pathview(gene.data=foldchanges,gene.idtype="kegg",
         pathway.id = "mmu04350",species="mmu", 
         out.suffix = "Kegg-TGFbeta-Day4.kegg",
         kegg.native = T, same.layer=T)

data(go.sets.mm)
data(go.subs.mm)
gobpsets = go.sets.mm[go.subs.mm$BP]

gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)

lapply(gobpres, head)


