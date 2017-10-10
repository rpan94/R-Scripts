library(ggplot2)
library(reshape2)
getwd()
setwd("/Users/bachum/Desktop/")
load("/Users/bachum/Desktop/HA-2/avgprof.RData")
myData <- as.data.frame(regcovMat)
glimpse(myData)
#subset required columns for the plot
colnames(myData)
Const <- myData[,10:19] 
glimpse(Const)

#class(myData)/Users/bachum/Desktop/Zmynd11/
#View(myData)
#Add an additional column variable describing the x-axis and divide into 1 to 101 row names undername Bins
#use Bin variable for reshaping the dataframe
Const$Bins <- 1:101
Const[1:3,]
colnames(Const)
#plot(myData[,1])
glimpse(Const)
#use Bin variable for reshaping the dataframe
d <- melt(Const, id.vars="Bins")
glimpse(d)
# Every variable on the ggplotsame plot
ggplot(d, aes(Bins,value, col=variable)) + 
  geom_line()
# Use span to control the "wiggliness" of the default loess smoother
# The span is the fraction of points used to fit each local regression:
# small numbers make a wigglier curve, larger numbers make a smoother curve.
ggplot(d, aes(Bins,value, col=variable)) + 
  geom_line() + geom_smooth(span = 0.3, se = F) + theme_bw() + theme(panel.grid = element_blank())
# rename the x-ais By adding the scale_x_discrete(limits=c(0,20,40,60,80,100),labels = c("-5000","TSS","33%","66%", "TES", "5000"))
ggplot(d, aes(Bins,value, col=variable)) + 
  geom_line() + geom_smooth(span = 0.3, se = F) + theme_bw() + theme(panel.grid = element_blank()) +
  scale_x_discrete(limits=c(0,20,40,60,80,100),labels = c("-5000","TSS","33%","66%", "TES", "5000"))

#making the alpha component zero retains only the smoothened version of line plot
plot <- ggplot(d, aes(Bins,value, col=variable)) + 
  geom_line(alpha=0) + geom_smooth(span = 0.3, se = F) + theme_bw() + theme(panel.grid = element_blank()) +
  scale_x_discrete(limits=c(0,20,40,60,80,100),labels = c("-5000","TSS","33%","66%", "TES", "5000"))
plot
#Add new labels to x-and y-axis
Bins<- plot + labs(x = "GeneBody", y = 'Reads Per Million (RPM)')
Bins
ISG
All
Constitu
library(gridExtra)
comboplot <- grid.arrange(All, Constitu,Bins,ISG, fig2, ncol =2,widths = c(1,1))
ggsave("combo_plot.png", comboplot, width = 15, dpi = 300)

# Change font options:
# X-axis label: bold, red, and 20 points
# X-axis tick marks: rotate 90 degrees CCW, move to the left a bit (using vjust,
#   since the labels are rotated), and 16 points
#plot1 + theme(axis.title.x = element_text(face="bold", colour="#990000", size=10),
#           axis.text.x  = element_text(angle=270, vjust=0.5, size=10))

#plot1+ theme(axis.text.y   = element_text(size=10),
#        axis.text.x   = element_text(size=10),
#        axis.title.y  = element_text(size=10),
#        axis.title.x  = element_text(size=10),
#        panel.background = element_blank(),
#        panel.grid.major = element_blank(), 
#        panel.grid.minor = element_blank(),
#        axis.line = element_line(colour = "black"),
#        panel.border = element_rect(colour = "black", fill=NA, size=1))

# Separate plots
fig <- ggplot(d, aes(Bins,value)) + 
  geom_line(alpha=0) + geom_smooth(span = 0.3, se = F) +
  facet_wrap(~variable) + theme_bw() + 
  theme(panel.grid = element_blank()) + 
  scale_x_discrete(limits=c(0,20,40,60,80,100),labels = c("-5000","TSS","33%","66%", "TES", "5000"))

fig2 <- fig + theme(axis.text.x=element_text(angle=65, hjust=1, size = 6), 
                                     axis.text.y=element_text(size = 8))
ggsave("Fig2.png", fig2, width = 8, dpi = 300)


