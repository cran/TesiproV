#' MVFOSM
#' @name MVFOSM
#' @param lsf LSF Definition, can be Expression or Function. Defined by the FLAG isExpression (see below)
#' @param lDistr List of Distributions
#' @param h If isExpression is False, than Finite Difference Method is used for partial deviation. h is the Windowsize
#' @param isExpression Boolean, If TRUE lsf has to be typeof expression, otherwise lsf has to be type of function()
#' @param debug.level If 0 no additional info if 2 high output during calculation
#'
#' @return beta, pf, design.point in x space, alphas, runtime
#' @references  FREUDENTHAL, A.M. Safety and the probability of structural failure. Am Soc Civil Eng Trans 1956; 121(2843):1337–97.
#' @author (C) 2021 - K. Nille-Hauf, T. Feiri, M. Ricker - Hochschule Biberach, Institut fuer Konstruktiven Ingenieurbau#'
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
      #lsf ist funktion mit lsf(x){x[1]...}
      #x vektor für
      #h ist bandbreite

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


