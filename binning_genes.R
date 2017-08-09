library(binr)
ls("package:binr")
?bins
?bins.quaRPKMiles
dt <- read.csv(file = "/Users/bachum/Desktop/Chao-NT-MEF-RPKM.csv", header = T)
colnames(dt)
barplot(dt$RPKM)
min(dt$RPKM)
max(dt$RPKM)
bins=seq(min(dt$RPKM),max(dt$RPKM),1)
hist(dt$RPKM, breaks=bins)
library(ggplot2)
ggplot(dt,aes(x=RPKM)) + geom_histogram(binwidth=50)

library("plotrix")
gap.barplot(dt$RPKM, c(1000,12000), dt$Transcript_ID)

library(reshape)
gap.barplot(dt$RPKM_FPKM,
            gap=c(200,12000),
            xlab="Genes",
            ytics=c(0,1000,2000,3000,4000:12000),
            ylab="FPKM",
            xaxlab=dt$X0h_FPKM,
            xaxt="n")
# xaxt="n" is eseRPKMiall to remove everything from x axis (e.g. a clean x axis)
# then define a axis using the following 
axis(side = 1, at = seq_along(mdata$Animal),mdata$Animal,tick = FALSE)
abline(h=seq(200,205,.001), col="white")  # hiding vertical lines
axis.break(axis=2,breakpos=202.5,style="slash") # break the left Y axis


cuts <- bins(dt$RPKM, target.bins = 10, minpts = 500)
head(cuts)
cuts$breaks <- bins.getvals(cuts)
cuts$binct
str(cuts)

barplot(cuts$binct)

