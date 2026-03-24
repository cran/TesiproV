#' @name MVFOSM
#' @title Mean-Value First-Order Second-Moment (MVFOSM)
#' @description
#' Mean-Value First-Order Second-Moment (MVFOSM) method for the
#' approximation of failure probabilities in structural reliability analysis.
#'
#' The MVFOSM method linearises the limit-state function at the mean values
#' of the basic random variables and estimates the reliability index
#' and probability of failure based on first-order moment information.
#'
#' This classical approach provides a computationally efficient
#' approximation but may be inaccurate for strongly nonlinear
#' limit-state functions.
#'
#' @param lsf LSF Definition, can be Expression or Function. Defined by the FLAG isExpression (see below)
#' @param lDistr List of Distributions
#' @param h If isExpression is False, than Finite Difference Method is used for partial deviation. h is the Windowsize
#' @param isExpression Boolean, If TRUE lsf has to be typeof expression, otherwise lsf has to be type of function()
#' @param debug.level If 0 no additional info if 2 high output during calculation
#'
#' @return MVFOSM returns an object containing the following elements:
#' \itemize{
#'   \item \code{beta}: Estimated reliability index.
#'   \item \code{pf}: Estimated probability of failure.
#'   \item \code{design.point}: Design point in the original \eqn{x}-space.
#'   \item \code{alphas}: Direction cosines (importance factors) in the standard normal space.
#'   \item \code{runtime}: Total runtime of the algorithm.
#' }
#' @references
#' Freudenthal, A. M. (1956).
#' Safety and the probability of structural failure.
#' \emph{Transactions of the American Society of Civil Engineers},
#' 121, 1337-1397.
#' @author (C) 2021-2026 K. Nille-Hauf, T. Feiri, M. Ricker, T. Lux -- Hochschule Biberach (until 2022),
#' TU Dortmund University - Chair of Structural Concrete (since 2023)
#' @export

MVFOSM <- function(lsf,
                   lDistr,
                   h=0.0001,
                   isExpression=FALSE,
                   debug.level){

  debug.TAG <- "MVFOSM"
  debug.print(debug.level,debug.TAG,c(TRUE), msg="Mean Value First Order Second Moment Method started...")
  tic<-proc.time()

  n_vars <- length(lDistr)
  m <- vector("numeric",n_vars)
  sd <- vector("numeric",n_vars)
  dx <- vector("numeric", n_vars)

  # Calculate Sd
  for (i in 1:n_vars){
    m[i] <- lDistr[[i]][[2]][1]
    sd[i] <- lDistr[[i]][[2]][2]
  }

  if(isExpression){
    s <- 0
    for(i in 1:n_vars){
      assign(lDistr[[i]][[1]]$name,m[i])
      d <- stats::D(lsf,lDistr[[i]][[1]]$name)
      dx[i] <- (eval(d))
      s <- s + (dx[i]^2)*(sd[i]^2)
    }
    sum_sd <- sqrt(s)
    lsf_mean <- eval(lsf)


  }else{
    finite_diff <- function(lsf,x,h){
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
    dx <- vector("numeric", n_vars)
    dx <- finite_diff(lsf,m,h)
    s<-0
    for (i in 1:n_vars) {
      s <- s + ((dx[i]^2) * (sd[i]^2))
    }
    sum_sd <- sqrt(s)

    # Evaluate LSF with means
    lsf_mean <- lsf(m)
  }

  # Estimate beta and pf
  beta <- lsf_mean/sum_sd
  pf <- stats::pnorm(-beta,0,1)

  alpha <- vector("numeric",n_vars)
  # Get alphas
  for(i in 1:n_vars)
  {
    alpha[i] <- -(dx[i]*sd[i])/sum_sd
  }

  designPoint <- vector("numeric",n_vars)
  # Get designpoints
  for(i in 1:n_vars){
    designPoint[i] <- m[i]-(alpha[i]*beta*sd[i])
  }

  cat("\n")
  duration<-proc.time()-tic

  output <- list(
    "method"="MVFOSM",
    "beta"=beta,
    "pf"=pf,
    "design.point.x"=designPoint,
    "alpha"=alpha,
    "runtime"=duration[1:5]
  )

  debug.print(debug.level,debug.TAG,c(duration), msg="Mean Value First Order Second Moment Method finished in: ")
  return(output)

  return(output)
}


