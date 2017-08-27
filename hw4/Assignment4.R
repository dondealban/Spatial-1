# Install packages and mask functions --------------------
rm(list=ls())
list.of.packages <- c("dplyr", "stringr","ggmap","spatstat","fields","gstat","lattice","geoR")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
select <- dplyr::select
rm(list.of.packages,new.packages)
# Loading and processing data --------------------
WIPP_transmissivity <- read.table(paste(getwd(),"/WIPP_transmissivity.txt",sep=""), header=TRUE, quote="\"")
wipp <- WIPP_transmissivity[-which(WIPP_transmissivity$logT < -10),]
# PROBLEM 1 --------------------
##Part (a)
(vg <- variogram(logT~1, ~x+y, data=wipp, cutoff=20,width=1))
(vg.cressie <- variogram(logT~1, ~x+y, data=wipp, cutoff=20,width=1,cressie=TRUE))
par(mfrow=c(1,2))
plot(vg,pch=19,col=1,ylab=expression(hat(gamma)(h)[M]),main="Empirical Semivariogram")
plot(vg.cressie,pch=19,col=2,ylab=expression(hat(gamma)(h)[C]),main="Cressie's Empirical Semivariogram")

#Defining the exponential (no nugget) function
exp.no_nug=function(t2,t3,h){
  t2*(1-exp(-h/t3))
}
#Defining the weighted regression sum of squares corresponding to the above semivariogram model
WRSS.exp_no_nug=function(param.vector,svar){
  gamma.mod=exp.no_nug(param.vector[1],param.vector[2],svar[,2])
  sum((svar[,1]/(gamma.mod^2))*((svar[,3]-gamma.mod)^2))
}
#Initial parameter guesses
exp.initial=c(2,2)
#Store result of optim
(param.exp_no_nug=optim(par = exp.initial, fn = WRSS.exp_no_nug, svar = vg,method="L-BFGS-B", lower=c(0,0),upper=c(Inf,Inf)))
#Store parameter estimates for theta that minimizes WRSS
theta.exp_no_nug <- param.exp_no_nug$par
#Store WRSS at the theta values for computation of BIC
wrss.exp_no_nug <- param.exp_no_nug$value
#Create data frame to allow for plotting of empirical svgm with overlaid svgm model 
df <- as.data.frame(cbind(distance=vg[,2],gamma.hat = vg[,3],gamma.mod = theta.exp_no_nug[1]*(1-exp(-vg[,2]/theta.exp_no_nug[2]))))
#Ggplot object
gp <- ggplot(df, aes(x=distance,y=gamma.hat)) + 
  theme_classic()+
  geom_point(size=3) + stat_function(colour="red",fun = exp.no_nug, args = c(t2=theta.exp_no_nug[1],t3=theta.exp_no_nug[2],h=df$x)) +
  ggtitle(label = "Empirical Svgm with Overlaid Exponential (No Nugget) Fit") +
  ylab(expression(hat(gamma)))
#Plot the empircial svgm
gp 

##Part (b)
#Defining the gaussian (no nugget) function
gauss.no_nug=function(t2,t3,h){
  t2*( 1-exp( -(h/t3)^2 ) )
}
#Defining the weighted regression sum of squares corresponding to the above semivariogram model
WRSS.gauss_no_nug=function(param.vector,svar){
  gamma.mod=gauss.no_nug(t2=param.vector[1],t3 = param.vector[2],h=svar[,2])
  sum((svar[,1]/(gamma.mod^2))*((svar[,3]-gamma.mod)^2))
}
#Initial parameter guesses
p.initial=c(2,2)	
#Store result of optim
param.gauss_no_nug=optim(par = p.initial, fn = WRSS.gauss_no_nug, svar = vg,method="L-BFGS-B", lower=c(0,0),upper=c(Inf,Inf))
#Store parameter estimates for theta that minimizes WRSS
theta.gauss_no_nug <- param.gauss_no_nug$par
#Store WRSS at the theta values for computation of BIC
wrss.gauss_no_nug <- param.gauss_no_nug$value
#Create data frame to allow for plotting of empirical svgm with overlaid svgm model 
df.gauss_no_nug <- df %>% select(-3) %>% mutate(gamma.mod =gauss.no_nug(t2=theta.gauss_no_nug[1],t3 = theta.gauss_no_nug[2],h=distance) )
#Ggplot object 
gp.gauss_no_nug <- ggplot(df.gauss_no_nug, aes(x=distance,y=gamma.hat)) + 
  theme_classic()+
  geom_point(size=3) + stat_function(colour="red",fun = gauss.no_nug, args = c(t2=theta.gauss_no_nug[1],t3=theta.gauss_no_nug[2],h=df.gauss_no_nug$x)) +
  ggtitle(label = "Empirical Svgm with Overlaid Gaussian (No Nugget) Fit") +
  ylab(expression(hat(gamma)))
#Plot the empirical svgm 
gp.gauss_no_nug 

##Part (c) 
#Defining the exponential (yes nugget) function
exp.nug <- function(t1,t2,t3,h){
  t1 + t2*(1-exp(-h/t3))
}
#Defining the weighted regression sum of squares corresponding to the above semivariogram model
WRSS.exp_nug=function(param.vector,svar){
  gamma.mod=exp.nug(param.vector[1],param.vector[2],param.vector[3],svar[,2])
  sum((svar[,1]/(gamma.mod^2))*((svar[,3]-gamma.mod)^2))
}
#Initial parameter guesses
p.initial=c(0,2,2)
#Store result of optim
param.exp_nug <- optim(par = p.initial, fn = WRSS.exp_nug, svar = vg,method="L-BFGS-B", lower=c(0,0,0),upper=c(Inf,Inf,Inf))
#Store parameter estimates for theta that minimizes WRSS
theta.exp_nug <- param.exp_nug$par
#Store WRSS at minimum theta for calculation of BIC 
wrss.exp_nug <- param.exp_nug$value
#Create data frame to allow for plotting of empirical svgm with overlaid svgm model 
df.exp_nug <- df %>% select(-3) %>% mutate(gamma.mod = exp.nug(theta.exp_nug[1], theta.exp_nug[2], theta.exp_nug[3],h=distance))
#Ggplot object 
gp.exp_nug <- ggplot(df.exp_nug, aes(x=distance,y=gamma.hat)) + 
  theme_classic()+
  geom_point(size=3) + stat_function(colour="red",fun = exp.nug, args = c(t1=theta.exp_nug[1],t2=theta.exp_nug[2],t3=theta.exp_nug[3],h=df.exp_nug$x)) +
  ggtitle(label = "Empirical Svgm with Overlaid Exponential (With Nugget) Fit") +
  ylab(expression(hat(gamma)))
gp.exp_nug

#Defining the gaussian (yes nugget) function
gauss.nug=function(t1,t2,t3,h){
  t1 + t2*( 1-exp( -(h/t3)^2 ) )
}
#Defining the weighted regression sum of squares corresponding to the above semivariogram model
WRSS.gauss_nug=function(param.vector,svar){
  gamma.mod=gauss.nug(param.vector[1],param.vector[2],param.vector[3],svar[,2])
  sum( (svar[,1]/(gamma.mod^2))*((svar[,3]-gamma.mod)^2) )
}
#Initial parameter guesses
p.initial=c(0,2,2)
#Store results of optim
param.gauss_nug=optim(par = p.initial, fn = WRSS.gauss_nug, svar = vg,method="L-BFGS-B", lower=c(0,0,0),upper=c(Inf,Inf,Inf))
#Store parameter estimates for theta that minimizes WRSS
theta.gauss_nug <- param.gauss_nug$par
#Store WRSS at minimum theta for calculation of BIC
wrss.gauss_nug <- param.gauss_nug$value
#Create data frame to allow for plotting of empirical svgm with overlaid svgm model 
df.gauss_nug <- df %>% select(-3) %>% mutate(gamma.mod = gauss.nug(theta.gauss_nug[1], theta.gauss_nug[2], theta.gauss_nug[3],h=distance))
#Ggplot object 
gp.gauss_nug <- ggplot(df.gauss_nug, aes(x=distance,y=gamma.hat)) + 
  theme_classic()+
  geom_point(size=3) + stat_function(colour="red",fun = exp.nug, args = c(t1=theta.gauss_nug[1],t2=theta.gauss_nug[2],t3=theta.gauss_nug[3],h=df.gauss_nug$x)) +
  ggtitle(label = "Empirical Svgm with Overlaid Gaussian (With Nugget) Fit") +
  ylab(expression(hat(gamma)))
gp.gauss_nug

##Part (d) 
K=dim(vg)[1]
calc.bic <- function(x){
  p <- length(x$par)
  K*log(x$value/K)+p*log(K)
}
bic.list <- sapply(X = list(param.exp_no_nug,param.gauss_no_nug,param.exp_nug,param.gauss_nug),FUN = calc.bic)
# PROBLEM 2 --------------------
##Part (a): Use gstat to fit exponential svgm
initial <- vgm(psill=2,model="Exp", range=2)
fvgm <- fit.variogram(vg,initial,fit.method = 2)

##Part (b) 
exp.no_nug_reml <- fit.variogram(vg, initial, 5)

#So, let's use fit.variogram. However, according to the help page, only the sill parameter is fitted. 
#How to go about the range? 

##Part (c) 
model_list <- levels(initial$model)
#Create function to take the models and calculate the BIC for each. 
fitter <- function(name, emp.svgm,m){
  v <- list()
  f <- list()
  w <- c()
  r=15
  k=0
  for(i in 1:length(name)){
    if(name[i] %in% c("Err","Int","Nug","Hol")){r=0}
    if(name[i] =="Mat"){k=1}
    v[[i]] <- vgm(psill=2.5,model = name[i], range=r,kappa=k)
    f[[i]] <- fit.variogram(emp.svgm,v[[i]],fit.method=2)
    w[i] <- attributes(f[[i]])$SSErr
  }
  K=dim(emp.svgm)[1]
  p <- 2
  if(m==T){f=NULL}
  return(list(BIC=K*log(w/K)+p*log(K), MODEL=f))
}
#Return BIC for a list of models 
fitter(c("Exp","Gau","Sph","Mat","Cir","Wav"),emp.svgm = vg,m=F)
#Create a function to return spherical semivariogram model
spherical=function(t2,t3,h){
  less.than=which(h<=t3)
  more.than=which(h>t3)
  gam <- rep(NA,length(h))
  gam[less.than] <- t2*(1.5*h[less.than]/t3-0.5*(h[less.than]/t3)^3)
  gam[more.than] <- (t2)*rep(1,length(more.than))
  return(gam)
}
#Storing the fitted parameters for the spherical model for these data 
params <- fitter(c("Sph"),emp.svgm = vg,m=F)$MODEL[[1]]
#Create a function to plot the above svgm model in the forthcoming ggplot object

plot_best_bic <- function(x){
  spherical(t2 = params$psill, t3 = params$range,h=x)
}
df_best_bic <- as.data.frame( cbind(distance=vg[,2],gamma.hat = vg[,3],gamma.mod = spherical(params$psill,params$range,vg[,2]) ) )
gp_best_bic <- ggplot(df_best_bic, aes(x=distance,y=gamma.hat)) + 
  theme_classic()+
  geom_point(size=3) + stat_function(colour="red",fun = plot_best_bic) +
  ggtitle(label = "Empirical Svgm with Overlaid Spherical (No Nugget) Fit") +
  ylab(expression(hat(gamma)))
#Plot the object
gp_best_bic
# PROBLEM 3 --------------------
vg1 <- variogram(logT~1, ~x+y,data=wipp,cutoff=18,width=1.5, alpha=c(0,90))
plot(vg1,pch=19,col=1,ylab=expression(hat(gamma)(h),sep=""), main="Empirical Semivariogram \n North/South and East/West")
##Part (b) 
#Creating separate svgm models for the different directions 
vg0 <- vg1[vg1$dir.hor==0,]
vg90 <- vg1[vg1$dir.hor==90,]

initial0 <- vgm(psill=2,model="Exp", range=5)
initial90 <- vgm(psill=4,model="Exp",range=15)

v0 <- fit.variogram(vg0,initial0,fit.method=2)
v90 <- fit.variogram(vg90,initial90, fit.method=2)
plotfun0 <- function(x){
  exp.no_nug(t2=v0$psill,t3=v0$range,h=x)
}
plotfun90 <- function(x){
  exp.no_nug(t2=v90$psill,t3=v90$range,h=x)
}
df0 <- as.data.frame(cbind(distance=vg0[,2],gamma.hat = vg0[,3],gamma.mod = v0$psill*(1-exp(-vg0[,2]/v0$range))))
gp0 <- ggplot(df0, aes(x=distance,y=gamma.hat)) + 
  theme_classic()+
  geom_point(size=3) + stat_function(colour="red",fun = plotfun0) +
  ggtitle(label = "Empirical Svgm with Overlaid Exponential (No Nugget) Fit\n North/South Direction") +
  ylab(expression(hat(gamma)))
#Plot the object
gp0
#Repeat with the data of the other direction 
df90 <- as.data.frame(cbind(distance=vg90[,2],gamma.hat = vg90[,3],gamma.mod = v90$psill*(1-exp(-vg90[,2]/v90$range))))
gp90 <- ggplot(df90, aes(x=distance,y=gamma.hat)) + 
  theme_classic()+
  geom_point(size=3) + stat_function(colour="red",fun = plotfun90) +
  ggtitle(label = "Empirical Svgm with Overlaid Exponential (No Nugget) Fit\n East/West Direction") +
  ylab(expression(hat(gamma)))
#Plot the object
gp90
##Part (c)
(eff.range0 <- -log(0.05)*v0$range)
(eff.range90 <- -log(0.05)*v90$range)


##Part(e)
#Geometric anisotropy is on pg48 of notes. Probably going to use geoR::coords.aniso
new_coords <- geoR::coords.aniso(wipp[,1:2],aniso.pars = c(90,v90$range/v0$range))
isot.data <- as.data.frame(cbind(new_coords,wipp$logT)) %>% rename(x=V1,y=V2,logT=V3)
ggplot(data=wipp,aes(x,y)) + 
  theme_classic()+
  geom_point(size=3) + 
  ggtitle(label="Original Coordinates")
ggplot(data=isot.data,aes(x,y)) + 
  theme_classic() + 
  geom_point(size=3,col="red") + 
  ggtitle(label="Transformed Coordinates")
#Fit semivariograms to the transformed data, see if they level off
vg2 <- variogram(logT~1, ~x+y,data=isot.data,cutoff=18,width=1.5)
d <- as.data.frame(cbind(distance=vg2[,2],gamma.hat = vg2[,3]))
ggplot(d, aes(x=distance,y=gamma.hat)) + 
  theme_classic()+
  geom_point(size=3) + 
  ggtitle(label = "Empirical Semivariogram of Transformed Coordinates") +
  ylab(expression(hat(gamma)))