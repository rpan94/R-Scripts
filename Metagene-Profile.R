getwd()
setwd("/Users/bachum/Desktop/")
load("/Users/bachum/Desktop/HA-2/avgprof.RData")
myData <- as.data.frame(regcovMat)
#class(myData)/Users/bachum/Desktop/Zmynd11/
#View(myData)
myData[1:3,]
colnames(myData)
#plot(myData[,1])
x <- c(1:101)

#########################################################################################
#How to fit a smooth curve to my data in R?
#http://stackoverflow.com/questions/3480388/how-to-fit-a-smooth-curve-to-my-data-in-r
#make a vector of values on x-axis for smmothing calculation
x <- c(1:101)

# 8 figures arranged in 4 rows and 2 columns
#par(mfrow=c(4,2))
#to remove x-axis values used xaxt='n'
#to make font normal used font=1
#to control font sixe of the axis and labels used cex
png(filename="ISG-MEFs.png",width=5,height=4,units="in",res=300)
par(mfrow=c(2,2), tcl=-0.5, family="Arial", mai=c(0.3,0.3,0.3,0.3))
#par(mgp=c(axis.title.position, axis.label.position, axis.line.position))
par(mgp=c(2,0.5,0))
par(family="Arial")
B1 <- smooth.spline(x, myData$`NT-B1`,spar=0.3)
B2 <- smooth.spline(x, myData$`NT-B2`,spar=0.3)
B3 <- smooth.spline(x, myData$`NT-B3`,spar=0.3)
B4 <- smooth.spline(x, myData$`NT-B4`,spar=0.3)
B5 <- smooth.spline(x, myData$`NT-B5`,spar=0.3)
B6 <- smooth.spline(x, myData$`NT-B6`,spar=0.3)
B7 <- smooth.spline(x, myData$`NT-B7`,spar=0.3)
B8 <- smooth.spline(x, myData$`NT-B8`,spar=0.3)
B9 <- smooth.spline(x, myData$`NT-B9`,spar=0.3)
B10 <- smooth.spline(x, myData$`NT-B10`,spar=0.3)



plot(B1,type = "l", lwd = 1.5, col = "black",
     ylab = "", xlab = "", xaxt='n',font=1, cex.axis = 0.8,
     cex.lab =0.3, ylim=c(0.02, 0.18))
#lines(Six,type = "l", lwd = 3, col = "red",ylab = "RPM", xlab = "GeneBody", xaxt='n')
#lines(Tfour,type = "l", lwd = 3, col = "blue",ylab = "RPM", xlab = "GeneBody", xaxt='n')
#axis(1, at=c(0,20,40,60,80,100), labels=c("-5000","TSS","33%","66%", "TES", "5000"),font=0.5, cex.axis = 0.5)
box(lwd=1.5)
#lo <- loess(myData$`Wt-NT-ZMYND11` ~ as.numeric(x), myData)
#lines(predict(lo))
#After plotting add simply the seconf data point by using line command
lines(B2,type = "l", lwd = 1.5, col = "red",ylab = "", xlab = "", xaxt='n')
lines(B3,type = "l", lwd = 1.5, col = "green",ylab = "", xlab = "", xaxt='n')
lines(B4,type = "l", lwd = 1.5, col = "pink",ylab = "", xlab = "", xaxt='n')
lines(B5,type = "l", lwd = 1.5, col = "orange",ylab = "", xlab = "", xaxt='n')
lines(B6,type = "l", lwd = 1.5, col = "violet",ylab = "", xlab = "", xaxt='n')
lines(B7,type = "l", lwd = 1.5, col = "brown",ylab = "", xlab = "", xaxt='n')
lines(B8,type = "l", lwd = 1.5, col = "cyan",ylab = "", xlab = "", xaxt='n')
lines(B9,type = "l", lwd = 1.5, col = "darkmagenta",ylab = "RPM", xlab = "", xaxt='n')
lines(B10,type = "l", lwd = 1.5, col = "blue",ylab = "", xlab = "", xaxt='n')
library(Hmisc)
minor.tick(nx=3,ny=3,tick.ratio=0.5)
#Add specific label information to X-axis
axis(1, at=c(0,20,40,60,80,100), labels=c("-5000","TSS","33%","66%", "TES", "5000"),font=0.8, cex.axis = 0.8)
#legend("topleft",c("B1","B2","B3","B4","B5","B6","B7","B8","B9","B10"),cex=.8, 
#       col=c("black","red","green","pink","orange","violet","brown","cyan","darkmagenta","blue"), 
#       horiz = T, inset = c(0.01, 0.01),bty = "n")
dev.off()
NT <- smooth.spline(x, myData$`NT-ISGConst`,spar=0.3)
Six <- smooth.spline(x, myData$`IFN—6h-ISGConst`,spar=0.3)
#par(mfrow=c(3,2))
plot(NT,type = "l", lwd = 2, col = "black",
     ylab = "", xlab = "", xaxt='n',font=1, cex.axis = 0.8,cex.lab =0.8, ylim = c(0.02,0.15))
lines(Six,type = "l", lwd = 2, col = "red",ylab = "RPM", xlab = "GeneBody", xaxt='n')
axis(1, at=c(0,20,40,60,80,100), labels=c("-5000","TSS","33%","66%", "TES", "5000"),font=0.8, cex.axis = 0.8)
#legend("topright",c("NT","IFN"),lty=c(1,1),lwd=c(2.5,2.5),col=c("black","red"))
box(lwd=1.5)
library(Hmisc)
minor.tick(nx=3,ny=3,tick.ratio=0.5)
#Add specific label information to X-axis
axis(1, at=c(0,20,40,60,80,100), labels=c("-5000","TSS","33%","66%", "TES", "5000"),font=0.8, cex.axis = 0.8)
dev.off()


#########################################################################################
