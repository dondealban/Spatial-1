# Installing packages and masking functions --------------------
rm(list=ls())
list.of.packages <- c("dplyr", "stringr","ggmap","spatstat","fields")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)
select <- dplyr::select
rm(list.of.packages,new.packages)
# Loading data --------------------
WIPP_transmissivity <- read.table("E:/MATH532/HW3/WIPP_transmissivity.txt", header=TRUE, quote="\"")