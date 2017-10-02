library(ggplot2)
library(dplyr)
library(ggrepel)
#VOLCANOPLOT SCRIPT USING ggplot2
###########################################################################################
results <- read.csv("/Users/bachum/Desktop/Massspec.csv", header=TRUE)
head(results)


results <- mutate(results, sig=ifelse(results$padj<0.01 & abs(log2FoldChange)>1,"1.5-Fold & FDR<0.01", "Not Sig"))

results <- mutate(results, sig=ifelse(results$FDR<0.1,"Sig","Not Sig"))


colnames(results)

p <- ggplot(results, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(col=sig)) + scale_color_manual(values=c("red", "black")) + 
  theme_bw(base_size = 10) 


p <- ggplot(results, aes(IFN, -log10(FDR))) +
  geom_point(aes(col=sig)) + scale_color_manual(values=c("red", "black")) + 
  theme_bw(base_size = 10) 



p
#X-Axis limit
p + xlim(0.01,100 ) + ylim(0, 10)

#Size and Transpeency
q <- p + xlim(0.01,100) +  ylim(0,5) + geom_point(alpha = 0.01,size=0.01)
q
q+geom_text(data=filter(results, FDR<0.1), aes(label= Gene_ID))
library(ggrepel)
#Select your gene of interest using the row number of excel sheet
#Size of the graph coord_fixed(ratio = 0.2)

r <- q+geom_text_repel(data=results[c(52,136,338,73,30,84,45,111,344,130), ], aes(label= Gene_ID))
r
r <- q+geom_text_repel(data=filter(results, pvalue<0.01 & abs(log2FoldChange)>10), aes(label= Id))
r
s <- r + geom_text_repel(data=results[c(1212), ], aes(label= Id))
s
library("ggthemes")
s+ theme_bw() 
s + theme_economist()
s + theme_wsj()
r+geom_text_repel(data=filter(results, FDR<0.01 & abs(logFC)>10), aes(label= X)) + 
  theme_bw() + labs(title = "Hira-KO Vs Wt IFN-Gamma", size = 10) + coord_fixed(ratio = 0.2)



###########################################################################################

# Hide panel borders and grid lines
# But change axis line
p + theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(size = 0.5, linetype = "solid",
                                   colour = "black"))


# Change the colors of plot panel background to lightblue
# and the color of grid lines to white
p + theme(
  panel.background = element_rect(fill = "lightblue",
                                  colour = "lightblue",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                  colour = "white"), 
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "white")
)

