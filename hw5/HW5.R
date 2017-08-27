# Install packages and mask functions --------------------
rm(list=ls())
load(paste(getwd(),"/.RData",sep=""))
list.of.packages <- c("sp", "dplyr", "stringr","spatstat","fields","gstat","lattice","ggplot2","knitr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)
select <- dplyr::select
rm(list.of.packages,new.packages)
# Get and clean WIPP_Transmissivity if not in existence -------------------- 
if(!exists("wipp",.GlobalEnv)){
  W <- read.table(paste(getwd(),"/WIPP_transmissivity.txt",sep=""), header=TRUE, quote="\"")
  wipp <- W %>% filter(logT != min(logT))
  rm(W)
}
# PROBLEM 1 --------------------
## Part (a) Fit the spherical svgm w/o nugget using the gstat package 
if(class(wipp) != "SpatialPointsDataFrame"){
  sp::coordinates(wipp) = ~x+y
}
vg <- gstat::variogram(logT~1,data=wipp)
initial <- gstat::vgm(psill=2.4,model="Sph",range=11)
fit.vg <- gstat::fit.variogram(vg, initial, fit.method=2)

# Functions to facilitate plotting
spherical=function(t2,t3,h){
  less.than=which(h<=t3)
  more.than=which(h>t3)
  gam <- rep(NA,length(h))
  gam[less.than] <- t2*(1.5*h[less.than]/t3-0.5*(h[less.than]/t3)^3)
  gam[more.than] <- (t2)*rep(1,length(more.than))
  return(gam)
}
sph_points <- function(x){
  spherical(t2=fit.vg$psill,t3=fit.vg$range,h=x)
}
# Creating dataframe to be passed in to ggplot
df <- vg %>% select(2,3) %>% as.data.frame %>% mutate(g = sph_points(dist)) %>% rename(distance = dist, gam.hat = gamma, gam.sph = g)

# Ggplot object
gp <- ggplot(df, aes(x=distance,y=gam.hat)) + 
  theme_classic()+
  geom_point(size=3) + stat_function(colour="red",fun = sph_points) +
  ggtitle(label = "Empirical Svgm with Overlaid Spherical (No Nugget) Fit") +
  ylab(expression(hat(gamma)))

# Plot the object
gp
## Part (b) 
k <- krige(logT~1,locations=wipp,newdata=wipp,model=fit.vg)
z <- wipp$logT
z.hat <- k$var1.pred
all.equal(z, z.hat)

## Part (c) 
x.grid=seq(min(wipp$x),max(wipp$x),len=100)
y.grid=seq(min(wipp$y),max(wipp$y),len=100)
GRID=expand.grid(x.grid,y.grid,KEEP.OUT.ATTRS = F) %>% `colnames<-`(c("gx","gy"))
coordinates(GRID) = ~gx+gy

kg <- krige(logT~1,locations=wipp,newdata=GRID,model=fit.vg)
par(mfrow=c(1,2))
quilt.plot(x=as.matrix(as.data.frame(GRID)), y = kg$var1.pred,xlab="x",ylab="y")
title("Gridded Predictions")

quilt.plot(x=as.matrix(as.data.frame(GRID)), y = kg$var1.var,xlab="x",ylab="y")
title("Variability of Predictions")

## Part (d) 
loocv_seq <- seq(1:40)
cv <- function(x, sdata){
  if(class(sdata) != "SpatialPointsDataFrame"){
    sp::coordinates(sdata) = ~x+y
  }
  t <- sdata[-x,]
  vg <- gstat::variogram(logT~1,data=t)
  fit.vg <- gstat::fit.variogram(vg, initial, fit.method=2)
  kcv <- krige(formula = logT~1,locations = t,newdata=sdata[x,],model=fit.vg,debug.level=0)
  return(kcv$var1.pred)
}
zhat <- sapply(loocv_seq,FUN = cv,sdata=wipp)
zobs <- wipp$logT
(rmse.cv <- sqrt( mean( (zobs - zhat)^2 ) ) )


kcv <- krige.cv(logT~1,locations=wipp,model=fit.vg,nfold=nrow(wipp))
logt_hat <- kcv$var1.pred
# Sanity check
if(all.equal(sqrt(mean((logt_hat - wipp$logT)^2)), sqrt( mean( ( kcv$residual )^2 ) ))) {
  rmse.cv <- sqrt( mean( ( kcv$residual )^2 ) )
}

## Part (e) 
par(mfrow=c(1,1))
sep <- (logt_hat - wipp$logT)^2
quilt.plot(x=as.data.frame(wipp) %>% select(-3), y = sep)
title("Absolute Error of Prediction")

# Plot predicted values NEXT to observed values. Identify worst prediction
# ddf <- as.data.frame(cbind(wipp$logT, logt_hat)) %>% 
#   mutate(sq.err = (logt_hat - wipp$logT)^2) %>% 
#   rename(Z.obs=V1, Z.pred=logt_hat,Sq.Err=sq.err) 
# ex <- subset(ddf, Sq.Err ==max(Sq.Err))
# ggplot(data=ddf, aes(x=Z.obs,y=Z.pred),size=5) + geom_point() +
#   theme_classic() +
#   geom_point(data=ex,color="red") +
#   geom_text(data= ex,mapping=aes(x=Z.obs,y=Z.pred,label="Most Wrong",color="red"),hjust=0.8,vjust=1.5) + guides(color=F) +
#   ggtitle(label = "Predicted vs. Observed logT with \n Worst Prediction Highlighted")

# Problem 2 ----------
# This is a sanity check
all.equal(dist(wipp@coords)[1], sqrt((wipp@coords[1,1] - wipp@coords[2,1])^2+(wipp@coords[1,2] - wipp@coords[2,2])^2))

# Part b
loocv_idw_seq <- seq(1:40)
cv.idw <- function(x, sdata,p){
  if(class(sdata) != "SpatialPointsDataFrame"){
    sp::coordinates(sdata) = ~x+y
  }
  t <- sdata[-x,]
  vg <- gstat::variogram(logT~1,data=t)
  fit.vg <- gstat::fit.variogram(vg, initial, fit.method=2)
  zhat_idw <- idw(formula = logT~1,locations = t,newdata=sdata[x,],idp=p,debug.level=0)
  return(zhat_idw$var1.pred)
}
loocv_idw_p1_resids <- wipp$logT - sapply(loocv_idw_seq,FUN = cv.idw,sdata=wipp,p=1)
loocv_idw_p2_resids <- wipp$logT - sapply(loocv_idw_seq,FUN = cv.idw,sdata=wipp,p=2)

(rmse.cv.idw1 <- sqrt( mean( loocv_idw_p1_resids^2 ) ) )
(rmse.cv.idw2 <- sqrt( mean( loocv_idw_p2_resids^2 ) ) )