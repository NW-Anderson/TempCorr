# ##### PARAMS #####
# install.packages(c("phytools","diversitree","geiger",'doSNOW'))
# # install.packages("Rmpi")
# # install.packages("doMPI", dependencies=TRUE)
# #  library(doMPI)
# library(R.utils) 
# library(doSNOW)
# library(foreach)
# cl<-makeCluster(3)
# on.exit(stopCluster(cl))
# library(phytools)
# library(diversitree)
# library(geiger)
# # library(doMC)
# library(profvis)
# # registerDoMC(10)
# source("fnx.R")
# 
# statistic = "mean"
# # statistic = "min"
# nul.model = "ind"
# obs.model = "dep" # not supported
# h1 = "greater"
# # h1 = "lesser"
# # h1 = 'twotail'
# nul.iter = 100 # number of null data points desired
# stoc.iter = 20 # number of stochastic maps per data point desired
# message = T
# plot = T
# TO <- NULL
# ##### PARAMS #####
# # setting the rate that is varied. .1 == no correlation
# higher values == stronger correlation
# oddrate <- seq(from=.1, to=1, length.out=3)
# different tree sizes
# ntaxa <- c(50, 100, 200)

# smp.size <- 100
# final.result <- matrix(,length(ntaxa), length(oddrate))

# for(j in 1:length(ntaxa)){
#   # different rates 
#   for(i in 1:length(oddrate)){
#     # repeated tests for each 
HPRC.pwr.fp.test <- function(ntaxa, oddrate){
  ##### PARAMS #####
  #install.packages("phytools")
  # nstall.packages("diversitree")
  # install.packages("geiger")
  library(R.utils) 
  library(doSNOW)
  library(foreach)
  cl<-makeCluster(28, type="SOCK")
  registerDoSNOW(cl)
  on.exit(stopCluster(cl))
  library(phytools)
  library(diversitree)
  library(geiger)
  library(profvis)
  source("fnx.R")
  
  statistic = "mean"
  # statistic = "min"
  nul.model = "ind"
  obs.model = "dep" 
  h1 = "greater"
  # h1 = "lesser"
  # h1 = 'twotail'
  nul.iter = 100 # number of null data points desired
  stoc.iter = 20 # number of stochastic maps per data point desired
  message = T
  plot = T
  TO <- NULL
  smp.size <- 100
  ##### PARAMS #####
  opts <- list(preschedule=FALSE)
  sig.tests <- foreach(h = 1:smp.size, .options.multicore=opts, .combine = 'c', .packages=c("phytools","diversitree","geiger")) %dopar% {
    source("fnx.R", local = TRUE)
    # ~9 hr
    # start <- Sys.time()
    # creating a random tree
    tree <- trees(pars=c(1,.1), type="bd", max.taxa=ntaxa)[[1]]
    tree$edge.length <- tree$edge.length/max(branching.times(tree))
    # creating simulated tip state and organizing into a data frame
    par <- rbind(c(0, .1, .5, 0), c(.1, 0, 0 ,.5), c(.5, 0, 0, oddrate), c(0, .5, .1, 0))
    # tests if each state occurs in the stochastic mapping
    good.data2 <- F
    while(good.data2 == F){ 
      # tests if the simulated tip states have an acceptable number of each state
      good.data1 <- F
      while(good.data1 == F){
        data <- as.integer(rTraitDisc(phy = tree, 
                                      model = par, 
                                      k = 4,
                                      states = 1:4,
                                      root.value = sample(1:4, 1)))
        # changes the 1 to 4 from rTraitDisc to 11,12,21,22 for entering testcorr function
        data1 <- data2 <- c()
        for(g in 1:length(data)){
          if(data[g] == 1){
            data1 <- c(data1, 1)
            data2 <- c(data2, 1)
          }else if(data[g] == 2){
            data1 <- c(data1, 1)
            data2 <- c(data2, 2)
          }else if(data[g] == 3){
            data1 <- c(data1, 2)
            data2 <- c(data2, 1)
          }else if(data[g] == 4){
            data1 <- c(data1, 2)
            data2 <- c(data2, 2)
          }
        }
        if(sum(data1 == 1)>.3*ntaxa && sum(data1 == 1)<.7*ntaxa &&
           sum(data2 == 1)>.3*ntaxa && sum(data2 == 1)<.7*ntaxa) good.data1 <- T
      }
      tip.dat <- data.frame(tree$tip.label,data1,data2)
      tester <- VecCom(tip.dat)
      if('12' %in% tester && '11' %in% tester && '21' %in% tester && '22' %in% tester) good.data2 <- T
    }
    # performing the tempcorr test on the simulated tree and data
    
    rap <- TestCorr(tree, 
                    tip.dat, 
                    statistic, 
                    h1, 
                    nul.model, 
                    obs.model,
                    nul.iter,
                    stoc.iter,
                    message,
                    plot)
    # interested in the number of tests that return a significant result when none 
    # are supposed to
    if(oddrate == .1){
      rap$`p-values` < .05/8
    }else{
      # only interested in the times it returns a significant result
      # in the type of transition with the oddrate
      rap$`p-values`[[3]] < .05
    }
    # end <- Sys.time()
  }
  sig.tests <- unlist(sig.tests)
  final.result <- sum(sig.tests) / smp.size
  fle.nme <- paste('finalresult_odd=',oddrate,'_ntaxa=',ntaxa,'.RData', sep = '')
  save(final.result, file = fle.nme)
  
  # final.result[i, j] <- sum(sig.tests)/smp.size
  # cat(' (',i,',',j,') ')
}
#   }
#   cat('\n\n')
# }
# colnames(final.result) <- c('50','100', '200')
# rownames(final.result) <- c('NULL(.1)', 'SLOW(.55)', 'FAST(1)')
# final.result

