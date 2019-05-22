# p-values that for transitions that dont occur should be NA
# or better yet we could ask what proportion of datasets
# dont show a dataset.  For example, if we see know F1F2 changes
# but 10% of simulated datasets dont have F1F2 changes then not very
# suprising but if only 1% of simulated datasets lack F1F2 changes then 
# that is suprising and should be a significant p-value

# Primary wrapper function is TestCorr 
# internal functions start around line 200

# the 'load' of TestCorr is iter^2 + iter. ie if iter = 100
# then there will be 10,100 stochastic mappiings total
# 100 per the 100 simulated tip states + 100 for the observed

TestCorr <- function(tree, 
                     tip.dat, 
                     statistic, 
                     h1, 
                     nul.model,
                     obs.model,
                     nul.iter,
                     stoc.iter,
                     message,
                     plot,
                     TO = NULL){
  
  # validation testing
  if(!h1 %in% c("greater", "twotail", "lesser")){stop("h1 should be: greater, twotail, or lesser indicating your hypothesis.")}
  if(!nul.iter == round(nul.iter)){stop('nul.iter should be an integer.')}
  if(!stoc.iter == round(stoc.iter)){stop('stoc.iter should be an integer.')}
  if(!message %in% c(T, F)){stop("Message should be true or false.")}
  # if(model != 'ind'){stop('Only the ind model is supported.')}
  if(!plot %in% c(T,F)){stop('Plot must be true or false.')}
  if(!statistic %in% c('mean', 'min')){stop('Statistic should be mean or min for the choice for the time between transitions.')}
  
  
  if(is.null(TO)){
    if(length(tree$tip.label) <= 50){
      TO <- 3.883202 * 2
    }else if(length(tree$tip.label) > 50 && length(tree$tip.label) <= 100){
      TO <- 8.987760 * 2
    }else if(length(tree$tip.label) > 100 && length(tree$tip.label) <= 200){
      TO <- 15.400414 * 2
    }else if(length(tree$tip.label) > 200 && length(tree$tip.label) <= 500){
      TO <- 19.31102 * 2
    }else{
      cat('Size of tree exceeds testing, defaulting TO to 100 seconds.')
      cat('\n\n')
      TO <- 100
    }
  }else{
    if(!TO == round(TO)){stop('TO should be an integer value for number of seconds until a simulated data set is abandoned
                                or left blank (NULL)')}
  }
    # creates a list of two matrices filled with 
    # the deltaM of each type for the desired number 
    # of simulations
    res <- TempCorr(tree = tree,
                    tip.dat = tip.dat,
                    statistic = statistic,
                    nul.model = nul.model,
                    obs.model = obs.model,
                    h1 = h1,
                    nul.iter = nul.iter,
                    stoc.iter = stoc.iter,
                    message = message,
                    TO = TO)
    
    
    # creates the histogram of the distribution of 
    # the null with the obs stat line overlayed
    if(plot == TRUE){
      for(i in 1:ncol(res$sim)){
        if(sum(!is.na(res$sim[,i])) != 0){
          hist(res$sim[,i],main = colnames(res$sim)[i], 
               breaks = 20, xlab = 'time between mutations')
          switch(statistic,
                 "mean" = abline(v = mean(res$obs[,i],
                                          na.rm = TRUE),
                                 col = 'red'),
                 "min"  = abline(v = min(res$obs[,i],
                                         na.rm = TRUE),
                                 col = 'red'))
        }
      }
    }
    # creates a vector containing the p val of finding the obs statistic on the null distribution
    # for each of the 8 mutation sequences of interest
    p <-  c()
    # keeps the warning messages from displaying multiple times if there are multiple problems
    two <- F
    three <- F
    # This will contain a vector of values indicating what flags to print for wierd cases
    caseflags <- list()
    for(i in 1:ncol(res$sim)){
      # splitting the possible p val computations into cases 
      # case 1: The observed trait value is NA while the simulated data set contains non NA values
      if(is.nan(mean(res$obs[,i],na.rm = TRUE)) == TRUE && sum(!is.na(res$sim[, i])) != 0){
        caseflags[[i]] <- c(1, sum(is.na(res$sim[, i]))/nul.iter, i)
        ptemp <- 'inc*'
        # case 2: The observed trait is Na and the simulated data set is completely NA
      }else if(is.nan(mean(res$obs[,i],na.rm = TRUE)) == TRUE && sum(!is.na(res$sim[, i])) == 0){
        caseflags[[i]] <- c(2)
        ptemp <- 'inc**'
        # Case 3: The transition occured in the observed but none of the null
      }else if(is.nan(mean(res$obs[,i],na.rm = TRUE)) == FALSE && sum(!is.na(res$sim[, i])) == 0){
        caseflags[[i]] <- c(3)
        ptemp <- 'inc***'
        # Case 4: The transition occured in the oberved and in some of the null (the valid cases)
      }else if(is.nan(mean(res$obs[,i],na.rm = TRUE)) == FALSE && sum(!is.na(res$sim[, i])) != 0){
        switch(paste(h1, statistic, collapse=""),
               'greater mean' = {
                 above <- sum(res$sim[, i] >= mean(res$obs[, i], na.rm = TRUE), na.rm = TRUE)
                 above <- above + sum(is.na(res$sim[, i]))
               },
               'greater min' =  {
                 above <- sum(res$sim[, i] >= min(res$obs[, i], na.rm = TRUE), na.rm = TRUE)
                 above <- above + sum(is.na(res$sim[, i]))
               },
               'lesser mean' =  above <- sum(res$sim[, i] <= mean(res$obs[, i], na.rm = TRUE), na.rm = TRUE),
               'lesser min' =   above <- sum(res$sim[, i] <= min(res$obs[, i], na.rm = TRUE), na.rm = TRUE),
               'twotail mean' = {
                 above <- sum(res$sim[, i] >= mean(res$obs[, i], na.rm = TRUE), na.rm = TRUE)
                 above <- above + sum(is.na(res$sim[, i]))
               },
               'twotail min' =  {
                 above <- sum(res$sim[, i] >= min(res$obs[, i], na.rm = TRUE), na.rm = TRUE)
                 above <- above + sum(is.na(res$sim[, i]))
               }
        )
        ptemp <- above / nul.iter
        if(h1 == 'twotail'){
          if(ptemp > .5){
            ptemp <- 1 - ptemp
          }
          ptemp <- 2 * ptemp
        }
      }
      p[i] <- ptemp
    }
    names(p) <- colnames(res$sim)
    if(message == TRUE){
      if(length(caseflags) > 0){
        cat('\n\n')
        cat('Inconclusive Warnings:')
        for(i in 1:length(caseflags)){
          if(!is.null(caseflags[[i]])){
            if(caseflags[[i]][1] == 1){
              cat('\n\n')
              cat('*: The ')
              cat(names(p)[caseflags[[i]][3]])
              cat(' transition did not occur in the observed mappings and did not in ')
              cat(caseflags[[i]][2] * 100)
              cat(' percent of the null simulations.')
            }else if(caseflags[[i]][1] == 2 && two == F){
              cat('\n\n')
              cat('**: This transition did not occur in the observed mappings nor any of the null.')
              two <- T
            }else if(caseflags[[i]][1] == 3 && three == F){
              cat('\n\n')
              cat('***: This transition occured in the observed but in none of the null.')
              three <- T
            }
          }
        }
      }
      cat('\n\n')
      cat('P val of observed statistics versus the simulated null')
      cat('\n\n')
      print(p)
      cat('\n\n')
    }
    final.result <- list()
    final.result[[1]] <- p
    final.result[[2]] <- res$obs
    final.result[[3]] <- res$sim
    names(final.result) <- c("p-values", "observed","null")
    return(final.result)
  }
  
  
  
  
  ### INTERNAL FUNCTIONS
  # Returns either the min or mean time of each type of trait change given 
  # ProcMap is the main function that takes a stochastic map
  # and the users choice of statistic and returns the estimate
  # of that statistic on that stochastic map
  ProcMap <- function(history, measure){
    ## Get paths from a tree
    paths <- getPaths(tree=history, type="branch")
    
    ## Get times from paths
    # copying the blank list from the GetTimes funciton
    all.times <- vector("list", length=8)
    names(all.times) <- c("F1F2", "F1R2", "R1F2", "R1R2",
                          "F2F1", "F2R1", "R2F1", "R2R1")
    # here a loop that runs one for each pathway
    for (i in (1:length(paths))){
      if(length(history$maps[paths[[i]]]) == 1){
        tempvec <- unlist(history$maps[paths[[i]]])
      }else{
        tempvec <- CollapseVector(unlist(history$maps[paths[[i]]]))
      }
      # GetTimes current path
      if(length(tempvec) > 1){
        new.times <- GetTimes(tempvec)
        # add values for current path to collection
        all.times <- mapply(c, all.times, new.times, SIMPLIFY=FALSE)
      }
    }
    ## get statistic for each type of change pair
    switch(measure,
           "mean" = times <- lapply(FUN = Mean, X = all.times),
           "min"  = times <- lapply(FUN = Min, X = all.times))
    names(times) <- c("F1F2", "F1R2", "R1F2", "R1R2",
                      "F2F1", "F2R1", "R2F1", "R2R1")
    return(times)
  }
  
  ## Combines two two state data vectors into a single 4 state vector
  VecCom <- function(data){
    ## Combining the two two state trait vectors into a single 4 state vector
    new.data <- apply(X=data[, 2:3], MARGIN=1, FUN=paste, collapse="")
    names(new.data) <- as.character(data[, 1])
    return(new.data)
  }
  
  # altered Mean and Min function to deal with empty vectors
  Mean <- function(x){
    if(length(x)>0){
      mean(x)
    }else{
      NA
    }
  }
  Min <-function(x){
    if(length(x)>0){
      min(x)
    }else{
      NA
    }
  }
  
  # provides a list of all the paths from root to tip
  getPaths <- function(tree, type="branch") {
    # extract the table
    tab <- tree$edge
    # id the tips - occur in right column but not left
    tips <- tab[, 2][!tab[, 2] %in% tab[, 1]]
    
    # could look at vector of node labels
    # should be just ntips
    
    # id the root - occurs in left column but not right
    root <- (tab[, 1][!tab[, 1] %in% tab[, 2]])[1]
    # create a list to store results in
    paths <- list()
    # a loop to go through all tips
    for (i in 1:length(tips)) {
      # pull the current tip
      x <- tips[i]
      # vector to store branches in
      if(type == "branch") bl <- vector()
      # checks to see if we found root yet
      while (x[length(x)] != root) {
        # which node leads to the most recently sampled node
        y <- tab[which(x[length(x)] == tab[, 2]), 1]
        # the row of tab is equivelant to the branch index so
        # we get our branch id this way
        if(type == "branch"){
          bl <- c(bl, which(x[length(x)] == tab[, 2]))
        }
        # store our new node
        x[(length(x) + 1)] <- y
      }
      # lets reverse it because I think root to tip
      if (type == "node") paths[[i]] <- rev(x)
      if (type == "branch") paths[[i]] <- rev(bl)
    }
    return(paths)
  }
  
  # gets time interval between changes
  GetTimes <- function(x){
    # create a list to hold results
    results <- vector("list", length = 8)
    names(results) <- c("F1F2", "F1R2", "R1F2", "R1R2",
                        "F2F1", "F2R1", "R2F1", "R2R1")
    # travel down vector
    for (i in 2:(length(x) - 1)){
      switch(paste(names(x)[(i-1):(i+1)], collapse=" "),
             "11 21 22" = results[[1]] <- c(results[[1]], as.numeric(x[i])),
             "12 22 21" = results[[2]] <- c(results[[2]], as.numeric(x[i])),
             "21 11 12" = results[[3]] <- c(results[[3]], as.numeric(x[i])),
             "22 12 11" = results[[4]] <- c(results[[4]], as.numeric(x[i])),
             "11 12 22" = results[[5]] <- c(results[[5]], as.numeric(x[i])),
             "21 22 12" = results[[6]] <- c(results[[6]], as.numeric(x[i])),
             "12 11 21" = results[[7]] <- c(results[[7]], as.numeric(x[i])),
             "22 21 11" = results[[8]] <- c(results[[8]], as.numeric(x[i]))
      )
    }
    return(results)
  }
  
  # collapses a stochastic hirstory
  CollapseVector <- function(x){
    z <- x[1]
    for(i in 2: length(x)){
      if(names(z)[length(z)] != names(x)[i]){
        z <- c(z, x[i])
      }else{
        z[length(z)] <- z[length(z)] + x[i]
      }
    }
    return(z)
  } 
  
  # creates the mappings for the observed and simulated data
  # returns two matrices containing the mutation times for obs and sim
  TempCorr <- function(tree, 
                       tip.dat, 
                       statistic, 
                       h1, 
                       nul.model,
                       obs.model,
                       nul.iter,
                       stoc.iter,
                       message,
                       TO){
    # this sets up our matrix that allows for state dependent rates
    # we will need to populate this change we want to allow for 
    # rates to be correlated or not while still accounting for closeness
    # on the tree.  
    switch(obs.model,
           "dep" =   rmat <- matrix(c(0,1,2,0,
                                      3,0,0,4,
                                      5,0,0,6,
                                      0,7,8,0), 4, 4, byrow=T),
           "ind"  =   rmat <- matrix(c(0,1,2,0,
                                       3,0,0,2,
                                       4,0,0,1,
                                       0,4,3,0), 4, 4, byrow=T))
    # convert the data to 4 state
    new.dat <- VecCom(data = tip.dat)
    # now lets make the stochastic map
    results <- vector(length=9)
    for(j in 1:stoc.iter){
      if(message == TRUE){
        cat("\014")
        cat(paste(j, "stochastic maps analyzed"))
      }
      recons <- make.simmap(tree, 
                            x = new.dat, 
                            model = rmat,
                            message = F)
      x <-ProcMap(history=recons, measure=statistic)
      if(j == 1) results <- unlist(x)
      if(j > 1) results <- rbind(results, unlist(x))
    }
    obs <- results
    # now we create our null distribution
    res.data <- c()
    trait1 <- tip.dat[,2]
    trait2 <- tip.dat[,3]
    names(trait1) <- names(trait2) <- tip.dat[,1]
    # creating the rate matrices from the data file
    Q1 <- ace(x = trait1, 
              phy = tree, 
              type = "discrete", 
              model = "ARD")
    Q2 <- ace(x = trait2, 
              phy = tree, 
              type = "discrete", 
              model = "ARD")
    # Begin simulating tip data for the tree based on the Q matrices found earlier
    # making the q matrix the right format for the r trait disc funtion
    qmat1 <- Q1$index.matrix
    qmat1[2,1] <- Q1$rates[1]  
    qmat1[1,2] <- Q1$rates[2] 
    qmat2 <- Q2$index.matrix
    qmat2[2,1] <- Q2$rates[1]  
    qmat2[1,2] <- Q2$rates[2]  
    # begin the simulations
    # this loop is for the i=1...nul.iter rows of the null distribution
    for (i in 1:nul.iter){
      if(message == TRUE){
        cat("\014")
        cat(paste(i, "null dataset analyzed"))
      }
      # val is the validity of the current simulation,  
      # val becomes TRUE when stoc.iter number of stochastic maps for the current
      # simulation's data have been analyzed
      val <- FALSE
      rm(recons)
      # this allows only valid mean simulations to be added to the result data structure
      while(val == FALSE){
        # valid is the validity of the simulated traits. We do no want a 
        # mapping where too few changes occur
        valid <- FALSE
        # create two sets of simulated tip characteristics 
        # based on the qmatrix of our observed data
        while (valid == FALSE){
          sim.trait1 <- as.integer(rTraitDisc(phy = tree, 
                                              model = qmat1, 
                                              k = 2,
                                              states = 1:2,
                                              root.value = sample(1:2, 1)))
          x <- sum(sim.trait1 == 1)/length(sim.trait1)
          if(x > .1 & x < .9) valid <- TRUE
        }
        valid <- FALSE
        while (valid == FALSE){
          sim.trait2 <- as.integer(rTraitDisc(phy = tree, 
                                              model = qmat2, 
                                              k = 2,
                                              states = 1:2,
                                              root.value = sample(1:2, 1)))
          x <- sum(sim.trait2 == 1)/length(sim.trait2)
          if(x > .1 & x < .9) valid <- TRUE
        }
        # Creating the combined 4 state data vector for the current simmulation
        sim.tip.data <- VecCom(cbind(names(trait1), sim.trait1, sim.trait2))
        # Creating the q matrix format going into the make.simmap function
        switch(nul.model,
               "dep" =   rmat <- matrix(c(0,1,2,0, # not coded above yet
                                          3,0,0,4,
                                          5,0,0,6,
                                          0,7,8,0), 4, 4, byrow=T),
               "ind"  =   rmat <- matrix(c(0,1,2,0,
                                           3,0,0,2,
                                           4,0,0,1,
                                           0,4,3,0), 4, 4, byrow=T))
        colnames(rmat) <- rownames(rmat) <- c("11", "12","21","22")
        keep <- colnames(rmat) %in% sim.tip.data
        rmat <- rmat[keep, keep]
        # sim data holds all of the data from the stochastic mappings per null data point
        sim.data <- c()
        # count is the current number of timeouts in a row 
        count <- 0
        # j stochastic mappings per each i row of the null data set
        for(j in 1:stoc.iter){
          if(message == TRUE && count < 5){cat('.')}
          # is the validity of the current stochastic mapping it changes 
          # to TRUE
          v <- FALSE
          # only adds valid stochastic mappings 
          while(v == FALSE){
            if(count < 5){
              # here we set recons = NA for later when it will be sorted out if the trycatch catches an error
              recons <- NA
              # create the stochastic mapping of the simulated tip data
              tryCatch({
                # we include a timeout to deal with the fringe case in which the simulated data would cause the following 
                # function to enter a loop
                suppressWarnings(recons <- withTimeout({make.simmap(tree, x=sim.tip.data, model=rmat, message = FALSE)}, 
                                                       timeout = TO, onTimeout = "error"))}, 
                error = function(e){if(message == TRUE){cat(' Timeout Reached ')} })
              count <- count + 1
              # Here we sort the data from this simulation before adding it to the results
              # allows the first iteration to work
              suppressWarnings(if(is.null(nrow(sim.data > 1)) && is.na(recons) == FALSE){
                # processes the simulation mapping and adds it to the results
                temp.data <- ProcMap(history = recons, measure = statistic)
                sim.data <- rbind(sim.data, unlist(temp.data))
                v <- TRUE
                count <- 0
                # this accounts for the case where there is a timeout on the first iteration and a couple other cases where 
                # try catch is thrown
              } else if(is.na(recons)){
                v <- FALSE
                # this is for valid simulations after the first
              } else{
                temp.data <- ProcMap(history = recons, measure = statistic)
                # checks that the current data is not a duplication of the previous
                if(sum(as.numeric(unlist(temp.data)), na.rm = TRUE) != 
                   sum(as.numeric(unlist(sim.data[nrow(sim.data),])), na.rm = TRUE)){
                  sim.data <- rbind(sim.data, unlist(temp.data))
                  v <- TRUE
                  count <- 0
                } else{
                  v <- FALSE
                  cat(' Duplicated Data Point ')
                }
              })
            }else{v <- TRUE
            } # if count
          } # while v
        } # for j

        if(count >= 5 && message == TRUE){
          cat('\n\n')
          cat(' Frequent Problems Encountered. Abandoning Simulated Data Set ')
          cat('\n\n')
        }
        # processing all of the stochastic mappings
        if(count < 5){
          tempy <- c()
          for(i in 1:ncol(sim.data)){
            switch(statistic,
                   "mean" =   wop <- mean(sim.data[,i],na.rm = TRUE),
                   "min"  =   wop <- min(sim.data[,i],na.rm = TRUE))
            tempy <- cbind(tempy, wop)
          }
          # adding the new null data point to the result data structure
          res.data <- rbind(res.data,tempy)
          val <- TRUE
        }
      } # while val
    }# for loop
    
    colnames(res.data) <- c(names(temp.data))
    sim <- res.data
    
    results <- list(obs, sim)
    names(results) <- c("obs", "sim") 
    return(results)
  }
  