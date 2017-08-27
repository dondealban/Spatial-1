# Install packages and mask functions --------------------
rm(list=ls())
load(paste(getwd(),"/.RData",sep=""))
list.of.packages <- c("sp", "dplyr", "stringr","spatstat","fields","gstat","lattice","ggplot2",
                      "knitr","spdep","maptools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)
select <- dplyr::select
rm(list.of.packages,new.packages)
data(columbus)
attach(columbus)
colum <- readShapePoly(system.file("etc/shapes/columbus.shp", package="spdep"[1]))
# Problem 1 ----------

## Part (a)

plot(colum)
spplot(colum,"CRIME",main="Burglary/Thefts per Thousand Households")
spplot(colum,"HOVAL",main="Housing Value (in $1,000's)")
spplot(colum,"INC",main="Household Income (in $1,000's)")
spplot(colum,"OPEN",main="Open Space in Neighborhood")
spplot(colum,"PLUMB",main="Percentage of Housing Units without Plumbing")
spplot(colum,"DISCBD",main="Distance to Central Business District")
