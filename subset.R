mm10 <- read.csv(file = "/Users/bachum/Desktop/mm10_annotations.csv", header = T, row.names = NULL)
head(mm10)
colnames(mm10)
#Subset a set of rows based on the row names
#Reas CSV first
selected_genes <- read.csv(file = "/Users/bachum/Desktop/Workbook2.csv", header = F, row.names = NULL)
class(selected_genes)
dim(selected_genes)
colnames(selected_genes)
#subset based on matching patternes using %in%
#Don't forget to use the coulmn name during subset
mm10_Diff <- mm10[mm10$txName %in% selected_genes$V1,]
dim(mm10_Diff)
mm10_Diff
write.csv(mm10_Diff,file = "/Users/bachum/Desktop/mm10_high_annotations.csv")
