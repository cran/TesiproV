#' @name MC_IS
#' @title MonteCarlo Simulation with importance sampling
#' @description Method to calculate failure probability for structural engineering using a simulation method with
#' importance sampling (a method to reduce the amount of needed samples)
#'
#' @param lsf objective function with limit state function in form of function(x) {x[1]+x[2]...}
#' @param lDistr Distributions in input space
#' @param cov_user The Coefficent of variation the simulation should reach
#' @param n_batch Size per batch for parallel computing
#' @param n_max maximum of iteration the MC should do - its like a stop criterion
#' @param use_threads determine how many threads to split the work (1=singlecore, 2^n = multicore)
#' @param sys_type Determine if parallel or serial system (in case MCIS calculates a system)
#' @param dataRecord If True all single steps are recorded and available in the results file afteron
#' @param beta_l In Systemcalculation: LSF´s with beta higher than beta_l wont be considered
#' @param densityType determines what distributiontype should be taken for the h() density
#' @param dps Vector of design points that sould be taken instead of the result of a FORM analysis
#' @param debug.level If 0 no additional info if 2 high output during calculation
#'
#' @return The results will be provided within a list with the following objects. Acess them with "$"-accessor
#' @return pf probablity of failure
#' @return pf_FORM probablity of failure of the FORM Algorithm
#' @return var variation
#' @return cov_mc coefficent of the monteCarlo
#' @return n_mc number of iterations done
#' @references DITLEVSEN O, MADSEN H. Structural reliability methods, vol. 178. New York: Wiley; 1996.
#' @references Spaethe, G.: Die Sicherheit tragender Baukonstruktionen, 2. Aufl. Wien: Springer, 1991. – ISBN 3-211-82348-4
#'
#' @import parallel
#' @import edfun
#' @author (C) 2021 - K. Nille-Hauf, T. Feiri, M. Ricker - Hochschule Biberach, Institut fuer Konstruktiven Ingenieurbau
#' @export


MC_IS<-function(lsf,lDistr,cov_user=0.05,n_batch=16,n_max=1e6,use_threads=8,
                sys_type="parallel",dataRecord=TRUE,beta_l=100,densityType="norm",dps=NULL,debug.level=0){

  # SystemProbablity
  if(is.list(lsf)){
    debug.TAG <- "MC_IS_System"
    debug.print(debug.level,debug.TAG,c(TRUE), msg="Monte-Carlo Simulation with Importance Sampling started...")
    tic<-Sys.time()

    # Amount of LSF´s
    n_lsfs <- length(lsf)

    # System
    res_form <- vector("list",n_lsfs)
    res_form.beta <- vector("numeric",n_lsfs)
    res_form.dp <- vector("list",n_lsfs)
    res_form.beta_l <- vector("numeric")


    # Var mapping
    var.names <- vector("character")
    var.map <- vector("numeric")
    var.per.lsf <- vector("numeric")
    vars <- list()
    for(i in 1:n_lsfs){
      n_vars <- length(lDistr[[i]])
      var.per.lsf[i] <- n_vars
      for(j in 1:n_vars){
        var <- lDistr[[i]][[j]]
        if(!var[[1]]$name %in% var.names){
          k <- length(var.names)+1
          var.names[k] <- var[[1]]$name
          vars[[k]] <- var
          var.map <- c(var.map,k)
        }else{
          var.map <- c(var.map,match(var[[1]]$name,var.names))
        }
      }
    }
    n_vars_unique <- length(vars)

    # var.in.lsf ist eine Liste mit n_lsfs Einträgen. In jedem Eintrag ist stehen die Indizes der Vars die verwendet werden
    k<-1
    var.in.lsf <- list()
    for(i in 1:n_lsfs){
      var.in.lsf[[i]] <- var.map[k:(k+var.per.lsf[i]-1)]
      k <- k+var.per.lsf[i]
    }

    # lsf.in.var ist eine List mit n_vars_unique Einträgen. In jedem Eintrag stehen die LSF´s in dem die Variable verwendet wird
    lsf.in.var <- list()
    for(i in 1:n_vars_unique){
      loc_vars <- vector()
      for(j in 1:n_lsfs){
        if(i %in% var.in.lsf[[j]]){
          loc_vars <- c(loc_vars,j)
        }
      }
      lsf.in.var[[i]] <- loc_vars
    }


    # FORM Analysis
    for(i in 1:n_lsfs){

      # if(is.numeric(dps[i])){ # Vorgegebene Designpunkte, Beta wird dennoch mit FORM ermittelt
      #   res_form[[i]]<- TesiproV::FORM(lsf[[i]],lDistr[[i]],n_optim=20, loctol = 0.0001)
      #   res_form.dp[[i]] <- dps
      # }else{
      #   if(dps[i]=="SORM"){ # Betas und DPs sollen mit SORM ermittelt werden
      #     res_form[[i]]<- TesiproV::SORM(lsf[[i]],lDistr[[i]])
      #     res_form.dp[[i]] <- res_form[[i]]$design.point_X
      #
      #   }else{ # Betas und DPs werden mit FORM ermittelt
      #     res_form[[i]]<- TesiproV::FORM(lsf[[i]],lDistr[[i]],n_optim=20, loctol = 0.0001)
      #     res_form.dp[[i]] <- res_form[[i]]$x_points
      #   }
      # }

      res_form[[i]]<- TesiproV::FORM(lsf[[i]],lDistr[[i]],n_optim=20, loctol = 0.0001)
      res_form.dp[[i]] <- res_form[[i]]$x_points
      res_form.beta[[i]] <- res_form[[i]]$beta
      if(res_form.beta[[i]] < beta_l){
        res_form.beta_l <- append(res_form.beta_l, res_form.beta[[i]])
      }
      if (!is.null(dps)){
        if(is.list(dps) && length(dps)==n_lsfs){
          if(!is.null(dps[[i]])){
            if(length(dps[[1]]) == length(res_form.dp[[i]])){
              res_form.dp[[i]] <- dps[[i]]
            }else{
              warning(sprintf("The amount (%i) of defined designpoints for lsf #%i does not fit to the required amount (%i)",
                      length(dps[[i]]),i,length(res_form.dp[[i]])))
            }
          }
        }else{
          warning("Input for Designpoint is not valid! (Not a list or not of the same length than provided lsf)")
        }
      }

    }

    # Zur Vereinfachung wird bisher noch nicht mit beta_l gerechnet!


    # Wichtungsfaktoren
    ai <- matrix(nrow=n_vars_unique, ncol=n_lsfs)
    for(i in 1:n_vars_unique){
      sum_beta <- sum(1/unlist(res_form.beta[lsf.in.var[[i]]]))
      ai[i,lsf.in.var[[i]]] <- 1/res_form.beta[lsf.in.var[[i]]]/sum_beta
    }


    dp <- matrix(nrow=n_vars_unique, ncol=n_lsfs)
    for(i in 1:n_lsfs){
      varids <- var.in.lsf[[i]]
      k <- 1
        for(j in 1:length(varids)){
          varid <- varids[j]
          dp[varid,i] <- res_form.dp[[i]][k]
          k <- k+1
        }
    }
    dp.means <- rowMeans(dp,na.rm = TRUE)


    hv_func <- list()
    hv_func2 <- list()
    dp.mean <- vector("numeric",n_vars_unique)
    samplingArea <- list()
    for(p in 1:length(vars)){

      hv_func[[p]] <-  eval(substitute(function(...){
           hv_loc <- 0
          for(u in 1:length(lsfvars[[lp]])){


            if(lvars[[lp]][[1]]$DistributionType=="lnorm" || lvars[[lp]][[1]]$DistributionType=="lt" || lvars[[lp]][[1]]$DistributionType=="weibull"){
              dfun <- "stats::dlnorm"

              Mean_ <- ldp[lp,lsfvars[[lp]][u]]
              Sd_ <-  lvars[[lp]][[1]]$Sd #lvars[[lp]][[1]]$Cov*Mean_

              p2 <- sqrt(log(1+(Sd_/(Mean_))^2))
              p1 <- log((Mean_)/sqrt(1+(Sd_/(Mean_))^2))

            }else{
              dfun <- "stats::dnorm"
              p1 <- ldp[lp,lsfvars[[lp]][u]]
              p2 <- lvars[[lp]][[1]]$Sd
            }

            # if(densityType=="origin"){
            #   dfun <- paste(lvars[[lp]][[1]]$DistributionPackage,"::d",lvars[[lp]][[1]]$DistributionType,sep="")
            # }else{
            #   dfun <- paste("stats::d",densityType,sep="")
            # }
            # hv_loc <- hv_loc + lai[lp,u]*dnorm(...,ldp[lp,u],lvars[[lp]][[1]]$Sd)

            # hv_loc <- hv_loc + lai[lp,lsfvars[[lp]][u]]*do.call((eval(parse(text=dfun))),list(...,ldp[lp,lsfvars[[lp]][u]],lvars[[lp]][[1]]$Sd))
            hv_loc <- hv_loc + lai[lp,lsfvars[[lp]][u]]*do.call((eval(parse(text=dfun))),list(...,p1,p2))
          }
          hv_loc
        }), list(lp=p,ldp=dp,lvars=vars,lai=ai,lsfvars=lsf.in.var))

      info.print(debug.TAG,debug.level,c("Creating sampling density",". May took a while... Number:","of"),c(var.names[p],p,n_vars_unique))
      samplingWidth <- 4
      xmin <- min(dp[p,],na.rm=TRUE)-samplingWidth*vars[[p]][[1]]$Sd
      xmax <- max(dp[p,],na.rm=TRUE)+samplingWidth*vars[[p]][[1]]$Sd
      # samplingArea[[p]] <- seq(vars[[p]][[1]]$q(0.000001),vars[[p]][[1]]$q(1-0.000001),length.out=2e4)+(dp.means[p]-vars[[p]][[1]]$mean)
      # qfun <- paste(vars[[p]][[1]]$DistributionPackage,"::q",vars[[p]][[1]]$DistributionType,sep="")

      # xmin <- do.call((eval(parse(text=qfun))),list(x,min(dp[p,]),vars[[p]][[1]]$Sd))
      # xmax <- do.call((eval(parse(text=qfun))),list(1-x,min(dp[p,]),vars[[p]][[1]]$Sd))
#
#       x <- 0.000001
#       xmin <- vars[[p]][[1]]$q(x)-abs((vars[[p]][[1]]$mean-min(dp[p,],na.rm = TRUE)))
#       if (!(vars[[p]][[1]]$p(-1)>0)){
#         if(xmin<0){
#           xmin <- 0.001
#         }
#       }
#       xmax <- vars[[p]][[1]]$q(1-x)+abs(vars[[p]][[1]]$mean-max(dp[p,],na.rm = TRUE))

      samplingArea[[p]] <- seq(xmin,xmax,length.out=1e4)
      hv_func[[p]] <- edfun(dfun = hv_func[[p]],x=samplingArea[[p]])
    }


    # for(p in 1:length(vars)){
    #   hist(hv_func[[p]]$rfun(1e6),main=vars[[p]][[1]]$name)
    #   q <- seq(0,1,0.01)
    #   plot(x=q,y=hv_func[[p]]$qfun(q),type="l",main=vars[[p]][[1]]$name)
    #   x <- samplingArea[[p]]
    #   plot(x=x,y=hv_func[[p]]$dfun(x),type="l",main=vars[[p]][[1]]$name)
    # }



    I_sum <- 0
    pf_i <- 0
    var_i <- 0
    n_sim <- 0
    cov<-Inf
    k <- 0
    beta <- 0

    data.pf <- vector("numeric",n_max)
    data.I_sum <- vector("numeric",n_max)
    data.var <- vector("numeric",n_max)
    data.cov <- vector("numeric",n_max)
    data.n_sim <- vector("numeric",n_max)
    data.time <- vector("numeric",n_max)

    while(1){
      v <- vector("numeric",n_lsfs)
      fx_i  <- vector("numeric",n_lsfs)
      hv_i  <- vector("numeric",n_lsfs)


      for (i in 1:n_vars_unique) {
        v[i] <- hv_func[[i]]$rfun(1)
        fx_i[i] <- vars[[i]][[1]]$d(v[i])
        hv_i[i] <- hv_func[[i]]$dfun(v[i])
      }

      fx <- prod(fx_i)
      hv <- prod(hv_i)

      I_i <- vector("numeric",n_lsfs)
      for(j in 1:n_lsfs){
        lsf_i <- lsf[[j]]
        I_i[j]<-as.numeric(lsf_i(v[var.in.lsf[[j]]])<0)
      }

      I <- as.numeric(sum(I_i)>0)
      I_sum <- I_sum + I
      qt <- fx/hv
      if(!(is.na(I)|is.nan(I)) && !(is.na(qt)|is.nan(qt))){
        pf_i <- pf_i + I * qt
        var_i <- var_i + I * qt^2
        n_sim  <- n_sim + 1

      k <- k+1
      if(dataRecord){
        data.I_sum[k] <- I_sum
        data.n_sim[k] <- n_sim
        data.time[k] <- Sys.time()-tic
      }

        if((n_sim>1 && !(is.nan(pf_i)|is.na(pf_i)) && pf_i>0)){
            pf <- (1/n_sim)*pf_i
            var <- (1/(n_sim-1))*(((1/n_sim)*var_i)-pf^2)
            cov <- sqrt(var)/pf
            beta <- -stats::qnorm(pf)
            if(dataRecord){
              data.pf[k] <- pf
              data.var[k] <- var
              data.cov[k] <- cov
            }
        }
      }
      info.print(debug.TAG,debug.level,c("I_sum","pf_i","pf","beta","cov", "nsim"),c(I_sum,pf_i,pf,beta,cov,n_sim))

      if(cov<cov_user){break;}
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
        "time"=data.time[range],
        "I_sum"=data.I_sum[range]
      )
    }else{
      df <- data.frame()
    }

    output<-list(
      "method"="MCIS_System",
      "beta"=-stats::qnorm(pf),
      "pf"=pf,
      "FORM_beta"=res_form.beta,
      "design_points"=res_form.dp,
      "ai"=ai,
      "var"=var,
      "cov_mc"=cov,
      "cov_user"=cov_user,
      "n_mc"=n_sim,
      "n_max"=n_max,
      "data"=df,
      "runtime"=duration
    )
    debug.print(debug.level,debug.TAG,c(duration), msg="Monte-Carlo Simulation with Importance Sampling finished in [s]: ")
    return(output)




  ##########################################################################SingleProbablity####################################################################
  # Actually MC_IS_SYS is able to calculate a single LSF correctly but is a bit more unefficient since it samples the quantile of the density function
  # instead of using the implemented versions
  # thats why this second implementation is faster and still used
    }else if(is.function(lsf)){

    debug.TAG <- "MC_IS_Single"
    debug.print(debug.level,debug.TAG,c(TRUE), msg="Monte-Carlo Simulation with Importance Sampling started...")
    tic<-Sys.time()

    # Initializing vars
    n_vars <- length(lDistr)
    I<-0

    if(n_batch>n_max){n_batch<-n_max}

    # Get Designpoint by FORM Analysis

    if(is.numeric(dps)){ # Vorgegebene Designpunkte, Beta wird dennoch mit FORM ermittelt
      res_form<- list()
      dp <- dps
    }else{
        res_form<- TesiproV::FORM(lsf,lDistr,n_optim=20, loctol = 0.0001)
        debug.print(debug.level,debug.TAG,as.character(res_form),"Results of FORM: ")
        dp<-res_form$x_points
    }


    # Function is driven via multi threading
    # n_batch determines how many samples should be taken each thread, if its to small the computaional affort to organize the multi threading
    # is to big (overhead-problem) and if n_batch is to big the steps of COV are too big and will overshoot the destinated cov by far (more threads used than necessary)
    mc_local<-function(x){
      v<-matrix(nrow=n_batch, ncol=n_vars)
      hv_i <- matrix(nrow=n_batch, ncol=n_vars)
      fx_i <- matrix(nrow=n_batch, ncol=n_vars)

      for (i in 1:n_vars){
        m <- dp[i]
        sd <- lDistr[[i]][[1]]$Sd
        if(densityType=="origin"){
          dfun <- paste(lDistr[[1]][[1]]$DistributionPackage,"::d",lDistr[[1]][[1]]$DistributionTyp,sep="")
          rfun <- paste(lDistr[[1]][[1]]$DistributionPackage,"::r",lDistr[[1]][[1]]$DistributionTyp,sep="")
        }else{
          dfun <- paste("d",densityType,sep="")
          rfun <- paste("r",densityType,sep="")
        }

        v[,i] <- do.call((eval(parse(text=rfun))),list(n_batch,m,sd))
        fx_i[,i] <- sapply(v[,i],lDistr[[i]][[1]]$d)

        hvd <- function(x){do.call((eval(parse(text=dfun))),list(x,m,sd))}

        hv_i[,i] <- sapply(v[,i],hvd)
      }

      fx <- apply(fx_i,1,prod)
      hv <- apply(hv_i,1,prod)

      I<-as.numeric(apply(v,1,lsf)<0)

      qt <- fx/hv

      pf_i <- I * qt
      var_i <- I * qt^2
      return(list(sum(pf_i), sum(var_i),sum(I)))
    }


    n_sim <- 0
    cov<-1
    pf_mc<-0
    var_mc<-0
    I_mc<-0

    data.pf_mc <- vector("numeric",n_max)
    data.var_mc <- vector("numeric",n_max)
    data.I_mc <- vector("numeric",n_max)
    data.pf <- vector("numeric",n_max)
    data.var <- vector("numeric",n_max)
    data.cov <- vector("numeric",n_max)
    data.n_sim <- vector("numeric",n_max)
    data.time <- vector("numeric",n_max)
    k <- 0
    while (1)
    {
      k <- k+1
      r<-parallel::mclapply(seq(1:use_threads),mc_local)

      for(i in 1:use_threads){
        pf_mc<-pf_mc+r[[i]][[1]]
        var_mc<-var_mc+r[[i]][[2]]
        I_mc<-I_mc + r[[i]][[3]]
      }

      n_sim<-n_sim + (length(r)*n_batch)

      if(dataRecord){
        data.pf_mc[k] <- r[[i]][[1]]
        data.var_mc[k] <- r[[i]][[2]]
        data.I_mc[k] <- r[[i]][[2]]
        data.n_sim[k] <- n_sim
        data.time[k] <- Sys.time()-tic
      }

      if(!(is.nan(pf_mc)|is.na(pf_mc)) & pf_mc>0 & n_sim >1){

        pf <- (1/n_sim)*pf_mc
        var <- (1/(n_sim-1))*(((1/n_sim)*var_mc)-pf^2)
        cov <- sqrt(var)/pf

        if(dataRecord){
          data.pf[k] <- pf
          data.var[k] <- var
          data.cov[k] <- cov
        }

        info.print(debug.TAG,debug.level,c("I_mc","pf_mc","var_mc", "pf", "cov", "nsim"),c(I_mc,pf_mc,var_mc, pf,cov,n_sim))
        debug.print(debug.level,debug.TAG,cov,"cov")
        debug.print(debug.level,debug.TAG,pf,"pf")
        debug.print(debug.level,debug.TAG,n_sim,"n_sims")
      }else{
        info.print(debug.TAG,debug.level,c("I_mc","pf", "cov", "nsim"),c(I_mc, "NA","NA",n_sim))
      }


      if(cov<cov_user){break;}
      if(n_sim>n_max){break;}
    }


    if(dataRecord){
      range <- 1:k
      df <- data.frame(
        "pf_mc"=data.pf_mc[range],
        "var_mc"=data.var_mc[range],
        "I_mc"=data.I_mc[range],
        "n_sim"=data.n_sim[range],
        "pf"=data.pf[range],
        "var"=data.var[range],
        "cov"=data.cov[range],
        "time"=data.time[range]
      )
    }else{
      df <- data.frame()
    }

    cat("\n")
    duration<-Sys.time()-tic

    output<-list(
      "method"="MCIS_Single",
      "beta"=-stats::qnorm(pf),
      "pf"=pf,
      "FORM_beta"=res_form$beta,
      "FORM_pf"=res_form$pf,
      "design_points"=res_form$x_points,
      "var"=var,
      "cov_mc"=cov,
      "cov_user"=cov_user,
      "n_mc"=n_sim,
      "n_max"=n_max,
      "data"=df,
      "runtime"=duration
    )
    debug.print(debug.level,debug.TAG,c(duration), msg="Monte-Carlo Simulation with Importance Sampling finished in: ")
    return(output)
  }
}
