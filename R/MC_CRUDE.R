

#' @name MC_CRUDE
#' @title Crude MonteCarlo Simulation
#' @description Method to calculate failure probability for structural engineering
#'
#' @param lsf objective function with limit state function in form of function(x) {x[1]+x[2]...}
#' @param lDistr list ob distribiutions regarding the distribution object of TesiproV
#' @param cov_user The Coefficent of variation the simulation should reach
#' @param n_max maximum of iteration the MC should do - its like a stop criterion
#' @param n_batch Size per batch for parallel computing
#' @param use_threads Number of threads for parallel computing, use_threds=1 for single core. Doesnt work on windows!
#' @param dataRecord If True all single steps are recorded and available in the results file after on
#' @param debug.level If 0 no additional info, if 2 high output during calculation
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

MC_CRUDE<-function(lsf,lDistr,cov_user=0.05,n_batch=400,n_max=1e7,use_threads = 6,dataRecord=TRUE,debug.level=0){

  debug.TAG <- "MC_Crude"
  debug.print(debug.level,debug.TAG,c(TRUE), msg="Crude Monte-Carlo Simulation started...")
  tic<-proc.time()
  pf<-1
  I_n <- 0
  pf_i <- 0
  var_i <- 0
  n_sim<-0
  n_vars<-length(lDistr)
  v<-matrix(nrow=n_batch,ncol=n_vars)
  cov_mc <- +Inf
  k <- 0

  if(dataRecord){
    l_size=n_max/(use_threads*n_batch)
    data.pf <- vector("numeric",l_size)
    data.var <- vector("numeric",l_size)
    data.cov <- vector("numeric",l_size)
    data.n_sim <- vector("numeric",l_size)
    data.time.user <- vector("numeric",l_size)
    data.time.sys <- vector("numeric",l_size)
    data.time.elapsed <- vector("numeric",l_size)
    data.time.user.child <- vector("numeric",l_size)
    data.time.sys.child <- vector("numeric",l_size)
  }

  #mc_local is used for parallelisation purpose
  if(use_threads>1){
    RNGkind(kind="L'Ecuyer-CMRG")
  }


  mc_local <- function(x){
    for(i in 1:n_vars){
      v[,i]<-lDistr[[i]][[1]]$r(n_batch)
    }
    I<-sum(as.numeric(apply(v,1,lsf)<0))
    return(I)
  }

  while(1){

    # create Realisations
    if (Sys.info()[[1]]=="Windows"){ #Windows is not able to fork - parallelisation under windows os not efficient
      I_n <- I_n + sum(unlist(parallel::mclapply(seq(1,use_threads),mc_local, mc.cores=1)))
    }else{ # parallelisation for unix based platforms (macOS, Linux etc.)
      I_n <- I_n + sum(unlist(parallel::mclapply(seq(1,use_threads),mc_local, mc.set.seed = TRUE, mc.cores=use_threads)))
      stats::runif(1)
    }
    n_sim <- n_sim + use_threads*n_batch

    k <- k+1
    if(dataRecord){
      data.n_sim[k] <- n_sim
      p_stamp <- proc.time()-tic
      data.time.user[k] <- p_stamp[1]
      data.time.sys[k] <- p_stamp[2]
      data.time.elapsed[k] <- p_stamp[3]
      data.time.user.child[k] <- p_stamp[4]
      data.time.sys.child[k] <- p_stamp[5]
    }


    #Calculate stochastics
    debug.print(debug.level,debug.TAG, I_n, "I_n: ")
    if(I_n>0){
      pf <- I_n/n_sim
      var <- 1/(n_sim-1)*((1/n_sim)*I_n-pf^2)
      cov_mc <- sqrt(var)/pf

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

    if(!is.numeric(cov_mc)){
      cov_mc <- 1
      warning("The simulation broke. Something went wrong with the calculation of pf or var.")
      break;
    }

    if(cov_mc<cov_user){break;}
    if(n_sim>n_max){break;}
  }

  cat("\n")
  duration<-proc.time()-tic

  if(dataRecord){
    range <- 1:k
    df <- data.frame(
      "n_sim"=data.n_sim[range],
      "pf"=data.pf[range],
      "var"=data.var[range],
      "cov"=data.cov[range],
      "time.user"=data.time.user[range],
      "time.sys"=data.time.sys[range],
      "time.elapsed"=data.time.elapsed[range],
      "time.user.child"=data.time.user.child[range],
      "time.user.sys"=data.time.sys.child[range]
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
    "runtime"=duration[1:5]
  )

  debug.print(debug.level,debug.TAG,c(duration), msg="Crude Monte-Carlo Simulation finished in [s]:  ")
  return(output)
}
