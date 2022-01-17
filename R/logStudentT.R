#' Density Function for logarithmic student T distritbution
#'
#' @param x quantiles
#' @param m mean (1. parameter)
#' @param s standard deviation (2. parameter)
#' @param n 3. paramter
#' @param nue degrees of freedom
#' @return density
#'
#' @examples
#' dlt(0.5,3,6,2,5)
#'
#' @author (C) 2021 - K. Nille-Hauf, T. Feiri, M. Ricker - Hochschule Biberach, Institut fuer Konstruktiven Ingenieurbau
#' @export
#'
dlt <- function(x, m, s, n, nue){


  #PDF
  #v = nue
  #
  #d <- (gamma((nue + 1)/2) / (x*gamma(nue/2)*sqrt(pi*nue)*s))*(1+1/nue*((log(x)-m)/s*sqrt(n/(n+1)))^2)^(-(nue+1)/2)
  if(x>0){
    s_trans <- s/sqrt(n/(n+1))
    x_trans <- (log(x)-m)/s_trans

    d <- 1/(x*s_trans)*stats::dt(df=nue, x=x_trans)
  }else{
    d <- 0
  }

  return(d)

}



#' Probablity Function for logarithmic student T distritbution
#'
#' @param q quantiles
#' @param m mean (1. parameter)
#' @param s standard deviation (2. parameter)
#' @param n 3. paramter
#' @param nue degrees of freedom
#'
#' @return density
#' @author (C) 2021 - M. Ricker, K. Nille-Hauf, T. Feiri - Hochschule Biberach, Institut fuer Konstruktiven Ingenieurbau
#' @export
#'
plt <- function(q, m, s, n, nue){
  if(q>0){
    s_trans <- s/sqrt(n/(n+1))
    q_trans <- (log(q)-m)/s_trans
    p <- stats::pt(q_trans, df=nue)
  }else{
    p <- 0
  }
  return(p)
}

#' Quantil Function for logarithmic student T distritbution
#'
#' @param p probablity
#' @param m mean (1. parameter)
#' @param s standard deviation (2. parameter)
#' @param n 3. paramter
#' @param nue degrees of freedom
#'
#' @author (C) 2021 - M. Ricker, K. Nille-Hauf, T. Feiri - Hochschule Biberach, Institut fuer Konstruktiven Ingenieurbau
#' @return quantile
#' @export
#'
qlt <- function(p, m, s, n, nue){
  s_trans <- s/sqrt(n/(n+1))
  q <- exp(stats::qt(df=nue, p=p)*s_trans+m)
  return(q)
}

#' Random Realisation-Function for logarithmic student T distritbution
#'
#' @param n_vals number of realisations
#' @param m mean (1. parameter)
#' @param s standard deviation (2. parameter)
#' @param n 3. paramter
#' @param nue degrees of freedom
#' @author (C) 2021 - M. Ricker, K. Nille-Hauf, T. Feiri - Hochschule Biberach, Institut fuer Konstruktiven Ingenieurbau
#' @return random number
#' @export
#'
rlt  <- function(n_vals, m, s, n, nue){
  q <- stats::runif(n_vals)
  r <- qlt(q, m, s,n,nue)

  # Implement socalled natural boundaries
  # 5 times sd left and right...
  if(r<0){
    r <- 0
  }
  return(r)
}



# m <- 3.85
# s <- 0.09
# n <- 3
# nue <- 6
# x <- seq(1,100,0.1)
# d_st <- dlstudentt(x,m,s,n,nue)
# plot(x=x, y=d_st, type="l")
#
#
# q <- seq(1,100,0.1)
# p_st <- plstudentt(q, m, s,n,nue)
# plot(x=q, y=p_st, type="l")
#
#
#
# q <- seq(0,1,0.01)
# q_st <- qlstudentt(q, m, s,n,nue)
# plot(x=q, y=q_st, type="l")
#
# n_vals <- 1000000
# r_st <- rlstudentt(n_vals, m, s,n,nue)
# mean(r_st)
# plot(r_st)
#
#
# dlstudentt(qlstudentt(0.5, m , s, n, nue),m,s,n,nue)
# plstudentt(qlstudentt(0.5, m , s, n, nue), m, s, n, nue)
# qlstudentt(0.5, m , s, n, nue)
#
# n_vals <- 1000000
# r_st <- rlstudentt(n_vals, m, s,n,nue)
# mean(r_st)
# plot(r_st)
