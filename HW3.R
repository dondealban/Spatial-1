#https://github.com/ScouserInTrousers/Spatial.git
# Install packages and mask functions --------------------
rm(list=ls())
list.of.packages <- c("dplyr", "stringr","ggmap","spatstat","fields","gstat","lattice")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
select <- dplyr::select
rm(list.of.packages,new.packages)
# Functions to be used --------------------
svgm <- function(df,bin=NULL,cut=1){
  par(mfrow=c(1,2))
  x <- df[,1]
  y <- df[,2]
  d <- as.matrix(dist(cbind(x,y)))
  z <- df[,3]
  p <- gstat::variogram(z~1, loc = ~x+y, data=df,cloud=T,cutoff=max(d)/cut,width=bin,alpha=c(0,90))
  plot(p,ylab=expression((Z(s[i])-Z(s[j]))^2))
}
# Loading data --------------------
WIPP_transmissivity <- read.table("E:/MATH532/HW3/WIPP_transmissivity.txt", header=TRUE, quote="\"")
# Part 1 --------------------
#Plot the data, comment on any trends
ggplot(data=WIPP_transmissivity,aes(x=x,y=y)) + 
  geom_point() +
  geom_text(data = WIPP_transmissivity, aes(x,y, label = seq(1:length(x))), hjust = 2)
#Plot the semivariogram cloud for ALL pairs of points 
svgm(WIPP_transmissivity)

#How does the semivariogram cloud change when the large value of -10.12 is removed?
wipp <- WIPP_transmissivity[-which(WIPP_transmissivity$logT < -10),]
svgm(wipp)

#Investigate anisotropy in the data 
d <- as.matrix(dist(cbind(wipp$x,wipp$y)))
vg1=variogram(logT~1, ~x+y,data=wipp,cutoff=max(d)/2,width=4)
plot(vg1,pch=19,col=1,ylab=expression(paste("Estimated ",gamma(h),sep="")))


