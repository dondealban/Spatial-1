# Install packages and mask functions --------------------
rm(list=ls())
list.of.packages <- c("dplyr", "stringr","ggmap","spatstat","fields","gstat","lattice")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
select <- dplyr::select
rm(list.of.packages,new.packages)
# Functions to be used --------------------
svgm <- function(df){
  x <- df %>% select(one_of("x","lon","long"))
  y <- df %>% select(one_of("y","lat","latt"))
  z <- dff[,3]
  gstat::variogram(z~1, loc = ~x+y, data=df,cloud=T)
}
# Loading data --------------------
WIPP_transmissivity <- read.table("E:/MATH532/HW3/WIPP_transmissivity.txt", header=TRUE, quote="\"")
# Part 1 --------------------
#Plot the data, comment on any trends
ggplot(data=WIPP_transmissivity,aes(x=x,y=y)) + 
  geom_point() +
  geom_text(data = WIPP_transmissivity, aes(x,y, label = seq(1:length(x))), hjust = 2)
#Plot the semivariogram cloud for ALL pairs of points 
vgm.cloud <- variogram(logT~1, loc = ~x+y, data=WIPP_transmissivity,cloud=T)

#How does the semivariogram cloud change when the large value of -10.12 is removed?
wipp <- WIPP_transmissivity[-which(WIPP_transmissivity$logT < -10),]
