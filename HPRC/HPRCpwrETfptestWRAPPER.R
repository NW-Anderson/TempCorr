install.packages(c("phytools","diversitree","geiger",'doSNOW'))
install.packages("Rmpi")
install.packages("doMPI", dependencies=TRUE)
library(R.utils) 
library(doSNOW)
library(foreach)
cl<-makeCluster(28)
on.exit(stopCluster(cl))
library(phytools)
library(diversitree)
library(geiger)
library(profvis)

library(doMPI)
source("fnx.R")
source('HPRCpwrETfptestFUNCTION.R')

results <- foreach{}