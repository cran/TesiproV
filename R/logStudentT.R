#' Log-Student t Distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the log-Student t distribution.
#'
#' The cumulative distribution function is given by
#' \deqn{
#' F_X(x) =
#' F_{t_{\nu}}\left(
#' \frac{\ln(X/m)}{s}
#' \sqrt{\frac{n}{n+1}}
#' \right),
#' }
#' where \eqn{F_{t_{\nu}}(\cdot)} denotes the cumulative distribution
#' function of a log-Student t distribution with \eqn{\nu} degrees of freedom.
#'
#' The log-Student t distribution arises as the marginal
#' distribution of a log-normal distributed variable with conjugate
#' log-normal-gamma prior uncertainty in \eqn{m} and \eqn{s}.
#'
#' The hyperparameter vector \code{hyper.param} is defined as
#' \deqn{(m, s, n, \nu)}
#' where
#' \itemize{
#'   \item \eqn{m} is the mean value of an equivalent sample of size \eqn{n},
#'   \item \eqn{s} is the empirical standard deviation of an equivalent sample
#'         of size \eqn{\nu + 1},
#'   \item \eqn{\nu} denotes the degrees of freedom of the Student t distribution.
#' }
#'
#' For details see:
#' \url{https://www.jcss-lc.org/publications/jcsspmc/concrete.pdf}
#'
#' @name LogStudentT
#' @aliases logStudentT lt
#'
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n_vals Number of random values to generate.
#' @param hyper.param Numeric vector \eqn{(m, s, n, \nu)}.
#'
#' @return
#' \itemize{
#'   \item \code{dlt()} returns the density.
#'   \item \code{plt()} returns the distribution function.
#'   \item \code{qlt()} returns the quantile function.
#'   \item \code{rlt()} generates random deviates.
#' }
#'
#' @examples
#' hp <- c(3.4, 0.14, 3, 10)
#' #' dlt(0.5, hp)
#'
#' @author
#' (C) 2021-2026 K. Nille-Hauf, T. Feiri, M. Ricker, T. Lux --
#' Hochschule Biberach (until 2022), TU Dortmund University - Chair of Structural Concrete (since 2023)

#' @rdname LogStudentT
#' @export
dlt <- function(x, hyper.param){
  m <- hyper.param[1]
  s <- hyper.param[2]
  n <- hyper.param[3]
  v <- hyper.param[4]

  s_trans <- s/sqrt(n/(n+1))
  x_trans <- (log(x)-m)/s_trans

  d <- 1/(x*s_trans)*stats::dt(df=v, x=x_trans)
  return(d)
}

#' @rdname LogStudentT
#' @export
plt <- function(q, hyper.param){
  m <- hyper.param[1]
  s <- hyper.param[2]
  n <- hyper.param[3]
  v <- hyper.param[4]

  s_trans <- s/sqrt(n/(n+1))
  q_trans <- (log(q)-m)/s_trans
  p <- stats::pt(q_trans, df=v)
  return(p)
}

#' @rdname LogStudentT
#' @export
qlt <- function(p, hyper.param){
  m <- hyper.param[1]
  s <- hyper.param[2]
  n <- hyper.param[3]
  v <- hyper.param[4]

  s_trans <- s/sqrt(n/(n+1))
  q <- exp(stats::qt(df=v, p=p)*s_trans+m)
  return(q)
}

#' @rdname LogStudentT
#' @export
rlt  <- function(n_vals, hyper.param){
  m <- hyper.param[1]
  s <- hyper.param[2]
  n <- hyper.param[3]
  v <- hyper.param[4]

  q <- stats::runif(n_vals)
  r <- qlt(q, hyper.param)
  return(r)
}

