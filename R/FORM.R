#' @name FORM
#' @title First Order Reliablity Method
#' @description Method to calculate failure probability for structural engineering using approximation of limit state function
#' with linear part.
#' @param lsf objective function with limit state function in form of \code{function(R,E) {R-E}}. Supplied by a SYS_ object, do not supply yourself.
#' @param lDistr list ob distribiutions regarding the distribution object of TesiproV. Supplied by a SYS_ object, do not supply yourself.
#' @param n_optim number of opimaziationcycles (not recommended/need for lagrangian algorithms).
#' @param loctol Tolerance of the local solver algorithm
#' @param optim_type Optimaziationtypes. Available: Augmented Lagrangian Algorithm (use: "auglag"),
#' Rackwitz-Fissler Algorithm (use: "rackfies").
#' @param debug.level If 0 no additional info if 2 high output during calculation
#'
#' @return The results will be provided within a list with the following objects.
#' @return beta HasoferLind Beta Index
#' @return pf probablity of failure
#' @return u_points solution points
#' @return dy gradients
#'
#'
#' @references HASOFER AM, LIND NC. An exact and invarient first order reliability format. J Eng Mech Div Proc ASCE 1974;100(1):111–21.
#' @references Rackwitz-Fiessler: RACKWITZ R., FIESSLER B. Structural reliability under combined random load sequences. Comput Struct 1978;9(5), S. 489–94.
#' @references Optimised algorithm: YPMA, J., JOHNSON, S.G., BORCHERS, H.W., EDDELBUETTEL, D., RIPLEY, B., HORNIK K., CHIQUET, J., ADLER, A., nloptr: R Interface to NLopt. R package. 2020. Version 1.2.2.
#' @references Spaethe, G.: Die Sicherheit tragender Baukonstruktionen, 2. Aufl. Wien: Springer, 1991. – ISBN 3-211-82348-4
#' @import nloptr
#' @import pracma
#' @author (C) 2021 - K. Nille-Hauf, T. Feiri, M. Ricker - Hochschule Biberach, Institut fuer Konstruktiven Ingenieurbau
#' @export
#'
FORM<-function(lsf,lDistr,n_optim=10,loctol=1e-2,optim_type="rackfies",debug.level=0){

  debug.TAG <- "FORM_OP"
  debug.print(debug.level,debug.TAG,c(optim_type), msg="FORM Algorithm started with Solvingtype:")
  tic<-proc.time()

  n_vars <- length(lDistr)

  # transformation between normal domain to standard normal distributed space (gaussian domain)
  gausTrans <- function(a, toGaussian=TRUE){
    x<-vector(length=n_vars)
    y<-vector(length=n_vars)
    b<-vector(length=n_vars)

    if(toGaussian){
      #Entspricht Y zu X
      y<-a
      for(i in 1:n_vars){
        p <- stats::pnorm(y[i])-(1E-12)
        b[i]<- lDistr[[i]][[1]]$q(p)
      }
    }else{
      #Entspricht X zu Y
      x<-a
      for(i in 1:n_vars){
        x_ <- x[i] - lDistr[[i]][[1]]$X0
        b[i]<-stats::qnorm(lDistr[[i]][[1]]$p(x_))
      }
    }
    return(b)
  }


  # lsf in gaussian domain
  h<-function(y){
    return(lsf(gausTrans(y,TRUE)))
  }


  beta<-vector("numeric",n_optim)
  par_u<-matrix(ncol = n_vars,nrow = n_optim)
  par_x<-matrix(ncol = n_vars,nrow = n_optim)
  dy<-matrix(ncol = n_vars,nrow = n_optim)




    #Optimization
    if (optim_type=="auglag"){
      for(i in 1:n_optim){

        #Zufallsstartvektor erzeugen
        q_rand <- stats::qnorm(stats::runif(n_vars))

        #Standart-Optimierungsfunktion
        unorm <- function(u){sqrt(sum(u^2))}
        info.print(debug.TAG,debug.level,c("optim_run_i","/"),c(i,n_optim))
        res<-nloptr::auglag(x0 = q_rand,
                    fn = unorm,
                    heq=h,
                    localtol=loctol)

        beta[i]<-res$value
        par_u[i,]<-res$par
        par_x[i,]<-gausTrans(res$par,FALSE)
        dy[i,]<-grad(h,res$par)
      }
      if (res$convergence<0){
        warning("Error, Optimization failed!")
        return()
      }

      duration<-proc.time()-tic
      output<-list(
        "beta"=min(beta),
        "pf"=stats::pnorm(-abs(min(beta))),
        "u_points"=par_u[which.min(beta),],
        "x_points"=par_x[which.min(beta),],
        "dy"=dy[which.min(beta),],
        "optim_type"=optim_type,
        "runtime"=duration
      )

    }else if(optim_type=="rackfies"){

      n_sim <- 0
      x <- vector("numeric", n_vars)
      y <- vector("numeric", n_vars)
      dy <- vector("numeric", n_vars)
      m_s <- vector("numeric", n_vars)
      sd_s <- vector("numeric", n_vars)

      finite_diff <- function(lsf,x,h){
        dx <- vector("numeric",length(x))
        for (i in 1:length(x)) {
          x_lowerbound <- x
          x_upperbound <- x
          x_lowerbound[i] <- x[i]-h
          x_upperbound[i] <- x[i]+h

          dx_low <- lsf(x_lowerbound)
          dx_high <- lsf(x_upperbound)
          dx[i] <- (dx_high - dx_low)/(2*h)
        }
        return(dx)
      }

      for (i in 1:n_vars) {
        x[i] <- lDistr[[i]][[1]]$mean
      }
      debug.print(debug.level,debug.TAG,c("x_Start"),x)

      repeat{
        for (i in 1:n_vars) {

          x_ <- x[i]-lDistr[[i]][[1]]$X0
          sd_s[i] <- (1/lDistr[[i]][[1]]$d(x_)*stats::dnorm(stats::qnorm(lDistr[[i]][[1]]$p(x_))))
          m_s[i] <- x[i] - sd_s[i]*stats::qnorm(lDistr[[i]][[1]]$p(x_))

          y[i] <- (x[i] - m_s[i])/sd_s[i]

        }
        debug.print(debug.level,debug.TAG,c("m_stern"),m_s)
        debug.print(debug.level,debug.TAG,c("sd_stern"),sd_s)
        debug.print(debug.level,debug.TAG,c("y"),y)

        h_y <- lsf(x)
        debug.print(debug.level,debug.TAG,c("h(y) = lsf(x) ="),h_y)

        dy <- vector("numeric", n_vars)
        dy <- finite_diff(lsf,x,1E-5)*sd_s
        debug.print(debug.level,debug.TAG,c("dy"),dy)

        # alphas
        alpha <- vector("numeric",n_vars)
        for (i in 1:n_vars) {
          alpha[i] <- -dy[i]/((sum(dy^2))^(1/2))
        }
        debug.print(debug.level,debug.TAG,c("alphas"),alpha)

        beta <- (h_y-sum(y*dy))/((sum(dy^2))^(1/2))
        debug.print(debug.level,debug.TAG,c("beta"),beta)

        x_neu <- vector("numeric",n_vars)
        for (i in 1:n_vars) {
          x_neu[i] <- m_s[i] + alpha[i]*sd_s[i]*beta
        }
        debug.print(debug.level,debug.TAG,c("x_k+1"),x_neu)
        diff <- sum((x_neu - x)^2)^(1/2)

        x <- x_neu

        n_sim <- n_sim+1
        debug.print(debug.level,debug.TAG,c("n_sim","beta","diff","loctol"),c(n_sim,beta,diff,loctol))

        info.print(debug.TAG,debug.level,c("n_sim","beta","diff","loctol"),c(n_sim,beta,diff,loctol))
        if(diff <= loctol | n_sim > n_optim){ break }
      }


      duration<-proc.time()-tic
      output<-list(
        "method"="FORM",
        "beta"=beta,
        "pf"=stats::pnorm(-beta),
        "x_points"=x,
        "dy"=dy,
        "alpha"=alpha,
        "diff"=diff,
        "n_sim"=n_sim,
        "optim_type"=optim_type,
        "runtime"=duration[1:5]
      )
    }
  #Postprocessing


  debug.print(debug.level,debug.TAG,c(duration), msg="\nFORM Algorithm finished in  [s]: ")

  return(output)
}


