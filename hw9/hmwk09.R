## Spatial Statistics: Spring 2016
## Code to Accompany Hmwk#9

#Clear memory
rm(list=ls())
library(fields)
library(spatstat)


#####################################################
##Points collected from a plain empty square
#####################################################
data=read.csv(file="random_points_plain_square.csv",header=T)
head(data)
attach(data)

plot(x,y,col=rgb(0,0,0,0.35),pch=19)

x1=x/15
y1=y/15
plot(x1,y1,col=rgb(0,0,0,0.35),pch=19,xlab="",ylab="")
title("Open Empty Square",cex.main=1.5)


#Creating an object of type ppp in spatstat
win=owin(c(0,15),c(0,15))
data=ppp(x,y,window=win)  #Note, you will get some error message about some of the points being duplicates.  Just ignore this.

plot(data)

detach(data)

#####################################################
##Points collected from a square with a single point in the center
#####################################################
data=read.csv(file="random_points_center_dot.csv",header=T)
head(data)
attach(data)

summary(x)
summary(y)

plot(x,y,col=rgb(0,0,0,0.35),pch=19)


x2=x/13.8
y2=y/13.8
plot(x2,y2,col=rgb(0,0,0,0.35),pch=19,xlab="",ylab="")
title("Dot in Center of Square",cex.main=1.5)
#points(.5,.5,col=2,pch=19)



#Creating an object of type ppp in spatstat
win=owin(c(0,13.8),c(0,13.8))
data=ppp(x,y,window=win)  #Note, you will get some error message about some of the points being duplicates.  Just ignore this.

plot(data)

detach(data)

#####################################################
##Points collected from a square with a single point off center
#####################################################
data=read.csv(file="random_points_offcenter_dot.csv",header=T)
head(data)
attach(data)

plot(x,y,col=rgb(0,0,0,0.35),pch=19)


#Creating an object of type ppp in spatstat
win=owin(c(0,13.8),c(0,13.8))
data=ppp(x,y,window=win)  #Note, you will get some error message about some of the points being duplicates.  Just ignore this.

plot(data)

detach(data)