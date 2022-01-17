

#' @name MC_CRUDE
#' @title Crude MonteCarlo Simulation
#' @description Method to calculate failure probability for structural engineering
#'
#' @param lsf objective function with limit state function in form of function(x) {x[1]+x[2]...}
#' @param lDistr list ob distribiutions regarding the distribution object of TesiproV
#' @param cov_user The Coefficent of variation the simulation should reach
#' @param n_max maximum of iteration the MC should do - its like a stop criterion
#' @param n_batch Size per batch for parallel computing
#' @param use_threads TRUE for parallel computing, false for single core
#' @param dataRecord If True all single steps are recorded and available in the results file after on
#' @param debug.level If 0 no additional info if 2 high output during calculation
#'
#' @return The results will be provided within a list with the following objects. Acess them with "$"-accessor
#' @return pf probablity of failure
#' @return pf_FORM probablity of failure of the FORM Algorithm
#' @return var variation
#' @return cov_mc coefficent of the monteCarlo
#' @return n_mc number of iterations done
#'
#'
#' @references Spaethe, G.: Die Sicherheit tragender Baukonstruktionen, 2. Aufl. Wien: Springer, 1991. – ISBN 3-211-82348-4
#'
#' @import parallel
#' @author (C) 2021 - M. Ricker, K. Nille-Hauf, T. Feiri - Hochschule Biberach, Institut fuer Konstruktiven Ingenieurbau
#'
#' @export
#'

MC_CRUDE<-function(lsf,lDistr,cov_user=0.025,n_batch=100,n_max=1e7,use_threads = 64,dataRecord=TRUE,debug.level=0){

  debug.TAG <- "MC_Crude"
  debug.print(debug.level,debug.TAG,c(TRUE), msg="Crude Monte-Carlo Simulation started...")
  tic<-Sys.time()
  pf<-1
  I_n <- 0
  pf_i <- 0
  var_i <- 0
  n_sim<-0
  n_vars<-length(lDistr)
  v<-matrix(nrow=n_batch,ncol=n_vars)
  cov_mc <- +Inf
  k <- 0

  data.pf <- vector("numeric",n_max)
  data.var <- vector("numeric",n_max)
  data.cov <- vector("numeric",n_max)
  data.n_sim <- vector("numeric",n_max)
  data.time <- vector("numeric",n_max)

  #mc_local is used for parallelisation purpose
  mc_local <- function(x){
    for(i in 1:n_vars){
      v[,i]<-lDistr[[i]][[1]]$r(n_batch) #do.call(paste("r",lDistr[[i]][[1]], sep = ""),list(n_batch,lDistr[[i]][[2]][1],lDistr[[i]][[2]][2]))
    }
    I<-sum(as.numeric(apply(v,1,lsf)<0))
    return(I)
  }

  while(1){

    # #create Realisations
    if(use_threads>1){
      I_n <- I_n + sum(unlist(parallel::mclapply(seq(1,use_threads),mc_local)))
      n_sim <- n_sim + use_threads*n_batch
    }else{
      # This case is for comparsion between multicore and single core (mclapply vs apply...)
      for(i in 1:n_vars){
        v[,i]<-lDistr[[i]][[1]]$r(n_batch) #do.call(paste("r",lDistr[[i]][[1]], sep = ""),list(n_batch,lDistr[[i]][[2]][1],lDistr[[i]][[2]][2]))
      }
      #Evaluate LSF
      I<-as.numeric(apply(v,1,lsf)<0)
      n_sim<-n_sim+n_batch
      I_n <- I_n+sum(I)
    }

    k <- k+1
    if(dataRecord){
      data.n_sim[k] <- n_sim
      data.time[k] <- Sys.time()-tic
    }


    #Calculate stochastics
    debug.print(debug.level,debug.TAG, I_n, "I_n: ")
    if(I_n>0){
      pf <- I_n/n_sim
      var <- 1/(n_sim-1)*((1/n_sim)*I_n-pf^2)
      cov_mc <- sqrt(var)/pf;

      if(dataRecord){
        data.pf[k] <- pf
        data.var[k] <- var
        data.cov[k] <- cov_mc
      }

      info.print(debug.TAG,debug.level,c("I_n","pf", "cov", "nsim"),c(I_n, pf,cov_mc,n_sim))
      debug.print(debug.level,debug.TAG, pf, "Pf: ")
      debug.print(debug.level,debug.TAG, var, "var: ")
      debug.print(debug.level,debug.TAG, cov_mc, "CoV_mc: ")
      debug.print(debug.level,debug.TAG, cov_user, "cov_user: ")
      debug.print(debug.level,debug.TAG, n_sim, "n_sim: ")
      debug.print(debug.level,debug.TAG, n_max, "n_max: ")
    }else{
      info.print(debug.TAG,debug.level,c("I_n","pf", "cov", "nsim"),c(I_n, "NA","NA",n_sim))
    }

    if(cov_mc<cov_user){break;}
    if(n_sim>n_max){break;}
  }

  cat("\n")
  duration<-Sys.time()-tic

  if(dataRecord){
    range <- 1:k
    df <- data.frame(
      "n_sim"=data.n_sim[range],
      "pf"=data.pf[range],
      "var"=data.var[range],
      "cov"=data.cov[range],
      "time"=data.time[range]
    )
  }else{
    df <- data.frame()
  }

  output<-list(
    "method"="MCC",
    "beta"=-stats::qnorm(pf),
    "pf"=pf,
    "var"=var,
    "cov_mc"=cov_mc,
    "cov_user"=cov_user,
    "n_mc"=n_sim,
    "n_max"=n_max,
    "n_batch"=n_batch,
    "n_threads"=use_threads,
    "data"=df,
    "runtime"=duration
  )

  debug.print(debug.level,debug.TAG,c(duration), msg="Crude Monte-Carlo Simulation finished in [s]:  ")
  return(output)
}
