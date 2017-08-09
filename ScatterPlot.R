library(ggplot2)
library(dplyr)
library(ggrepel)
mydat <- read.csv(file = "/Users/bachum/Desktop/IFN-BetavsCtrl.complete.csv",
        header = T)
colnames(mydat)
library(ggplot2)
colnames(mydat)
results <- mutate(mydat, sig=ifelse(mydat$log2FoldChange>1, "Up", "Down"))
ggplot(mydat, aes (x= NT, y = IFN)) + geom_point()
#change shape size
ggplot(mydat, aes (x= NT, y = IFN)) + geom_point(shape = 16, size = 1)
#Color points by some other column

sub.data<-subset(mydat, log2FoldChange > 3, select=c(Id,log2FoldChange))


p <- ggplot(mydat, aes (x= NT, y = IFN, color = Diff.Genes, show.legend = FALSE)) +
  geom_point(shape = 16, size = 1) +
  theme_bw() + scale_color_gradient(low = "#999999", high = "red")
p + geom_text(data=filter(mydat, abs(log2FoldChange)>2), aes(label=Id))
library(ggrepel)
png(filename="PanelFigure5.png",width=4,height=3,units="in",res=300)
par(mfrow=c(1,1), tcl=-0.5, family="Arial", mai=c(0.3,0.3,0.3,0.3))
p + geom_text_repel(data=filter(mydat, abs(log2FoldChange)>10), 
                    aes(label=Id), size =2) +
  theme(panel.border = element_rect(fill=NA, colour = "black", 
                                    size=1)) 
dev.off()
getwd()
####################################################################################
plot(mydat$NT, mydat$IFN, pch =20)
colors()[1:25]
par(family="arial")
#Lets change the color of points
valcol <- mydat$Diff.Genes > 1 
#I Diff.Genes defined using excel differential gene column by a cut off abs(log2foldcahnge) > 0.1
head(valcol)
plot(mydat$NT, mydat$IFN, pch =20, col =rgb(0, 0, valcol))
#or using a specific conditions like below I cose differential genes
myvalues <- abs(mydat$log2FoldChange) > 0.8 #& mydat$padj < 0.1
valcol <- myvalues
#Plot Two Plots Side By Side Using Basic R
# Set graphical parameter `mfrow`
#par(mfrow = c(2, 2)) 
#png("Desktop/figa.png", width = 400, height = 350, res =100)
plot(mydat$NT, mydat$IFN, pch =20, col =rgb(valcol,0,0),
     ann=FALSE, cex=0.5)
# Define the position of tick marks
v1 <- c(0:8)
# Define the labels of tick marks
v2 <- c(0:8)
# Add axis to the plot 
axis(side = 1, at = v1, labels = v2,
     tck=-.1,tcl = -0.5,cex.axis=1.05,
     col.axis="black",font.axis=5)
axis(side = 2, at = v1, labels = v2,
     tck=-.1,tcl = -0.5,cex.axis=1.05,
     col.axis="black",font.axis=5)
library(Hmisc)
# Set minor tick mark
minor.tick(nx = 3.0, 
           ny = 3.0, 
           tick.ratio=0.5)
#How to change the box type on an R plot (bty ="n") or 
#using title() option we can label titles and x -alis
title(xlab="NT(log10(RPKM)+1)", ylab="IFN(log10(RPKM)+1)",cex.lab=.7) 
# add solid horizontal lines at y=1,5,7
#abline(h=2)
# add dashed blue verical lines at x = 1,3,5,7,9
#abline(v=2) 
#lm(mydat$IFN ~ mydat$NT)
# Grid
#grid(5,5)
# 10x6 cm
####################################################################################
