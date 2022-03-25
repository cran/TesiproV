#' @name SORM
#' @title Reliability Analysis at Biberach University of applied sciences
#' @description
#' # S. Marelli, and B. Sudret, UQLab: A framework for uncertainty quantification in Matlab, Proc. 2nd Int. Conf. on Vulnerability, Risk Analysis and Management (ICVRAM2014), Liverpool (United Kingdom), 2014, 2554-2563.
#' S. Lacaze and S. Missoum, CODES: A Toolbox For Computational Design, Version 1.0, 2015, URL: www.codes.arizona.edu/toolbox.
#' X. Z. Wu, Implementing statistical fitting and reliability analysis for geotechnical engineering problems in R. Georisk: Assessment and Management of Risk for Engineered Systems and Geohazards, 2017, 11.2: 173-188.
#' @param lsf objective function with limit state function in form of function(x) {x[1]+x[2]...}
#' @param lDistr list ob distribiutions regarding the distribution object of TesiproV
#' @param debug.level If 0 no additional info if 2 high output during calculation
#'
#' @return The results will be provided within a list with the following objects. Acess them with "$"-accessor
#' @return beta ... HasoferLind Beta Index
#' @return pf ... probablity of failure
#' @return u_points ... solution points
#' @return dy ... gradients
#' @import pracma
#' @author (C) 2021 -  T. Feiri, K. Nille-Hauf, M. Ricker - Hochschule Biberach, Institut fuer Konstruktiven Ingenieurbau
#' @references Breitung, K. (1989). Asymptotic approximations for probability integrals. Probabilistic Engineering Mechanics 4(4), 187–190. 9, 10
#' @references Cai, G. Q. and I. Elishakoff (1994). Refined second-order reliability analysis. Structural Safety 14(4), 267–276. 9, 10
#' @references Hohenbichler, M., S. Gollwitzer, W. Kruse, and R. Rackwitz (1987). New light on first- and second order reliability methods. Structural Safety 4, 267–284. 10
#' @references Tvedt, L. (1990). Distribution of quadratic forms in normal space – Applications to structural reliability. Journal of Engineering Mechanics 116(6), 1183–1197. 10
#'
#' @export

SORM <- function(lsf, # the limit-state function.
                 lDistr,
                 debug.level=0){ #  description of distribution of random variables

  debug.TAG <- "SORM_RABU"
  debug.print(debug.level,debug.TAG,c(TRUE), msg="SORM started...")
  tic <- proc.time()

  # Isoprobabilistic transformation to the standard Gaussian space U
  # Tlsf, convenient for users and required by mistral::FORM. AZC-style Tinv could be an alternative, e.g. Tinv <- function(X) { lsf(qnorm(pnorm(X),0.25,1)) }
  Tlsf <- function(X) {
    for(i in 1:length(X)){
      XU <- stats::pnorm(X[i])
      if(XU == 1){ XU = 1-1E-9}
      X[i] <- lDistr[[i]][[1]]$q(XU)
    }
    return(lsf(X))
  }


  # Orthonormal rotation matrix (implementation inspired by uq_gram_schmidt)
  GRAM_RABU <- function(alphas) {
    dimension <- length(alphas)
    Base<-diag(dimension)
    Base[,dimension] <- alphas
    NewBase <- Base
    Rot <- matrix( rep( 0, dimension*dimension), nrow = dimension)
    Rot[dimension,dimension] <- 1
    for (k in dimension:1){
      Z <- matrix( rep( 0, dimension), nrow = dimension)
      j <- k+1
      while(j<=dimension) {
        Z <- Z + pracma::dot(Base[,k],NewBase[,j]) * NewBase[,j]
        j <- j+1
      }
      NewBase[,k] <- Base[,k] - Z;
      Rot[k,k] <- pracma::Norm(NewBase[,k],2);
      NewBase[,k] <- NewBase[,k] / Rot[k,k];
      if (k<dimension) {
        Rot[k,(k+1):dimension] <- t(Base[,k]) %*% NewBase[,(k+1):dimension];
      }
    }
    return(t(NewBase))
  }




  # FORM
  res <- TesiproV::FORM(lsf, lDistr,loctol = 0.01,n_optim=10)

  # Alpha's (direction cosines) as components of the unit gradient vector:
  # The alphas can be determined through the gradient vector in the U space at the design-point u*
  gradient <- res$dy
  alphas <- gradient/pracma::Norm(gradient)

  # Matrix R
  R <- GRAM_RABU(alphas)

  # Hessian Matrix H
  H <- pracma::hessian(Tlsf, res$x_points) # Second-derivative Hessian matrix at the design-point u*


  # Matrix A as a product of R,H,R'divided by the lenght of the gradient vector at the design-point u*
  A <- R %*% H %*% t(R)
  A <- A/pracma::Norm(gradient)

  # Curvatures (found through the eigenvalues of the reduced matrix A
  if(!is.nan(A)){
    curvatures <- eigen(A[1:nrow(A)-1,1:nrow(A)-1])$values
    warning("No corvature found. Maybe bad designpoints of FORM.")
  }else{
    curvatures <- rep(0,length(gradient))
  }
  if (curvatures >= 1) {warning("FORM did not converge to the design point")}
  debug.print(0,debug.TAG,c("curvatures"),curvatures)

  # SORM: Reliability Index (Beta) and Probability of failure (Pf)
  # Breitung's formula
  SORMpfBreitung <- res$pf/sqrt(prod(1+res$beta*curvatures)) # Breitung's Formula
  SORMbetaBreitung <- -stats::qnorm(SORMpfBreitung) # Beta_SORM   norminv() in AZC, icdf() in UQL
  # NOTE: Check if the starting point was negative or not based on the origin point from FORM (as in the UQLab)

  # HohenBichler's formula
  CDFminusBeta <- stats::pnorm(-res$beta, 0, 1);
  RatioCDFPDF <- stats::dnorm(-res$beta, 0, 1)/CDFminusBeta;
  HohenbichlerVec <- 1./sqrt(1 + RatioCDFPDF*curvatures);
  SORMpfHB <- CDFminusBeta*prod(HohenbichlerVec);
  SORMbetaHB <- (-1)*stats::qnorm(SORMpfHB)

  # Tvedt's formula
  A11 <- 1/sqrt(prod(1+res$beta*curvatures))
  A12 <- res$beta*res$pf-stats::dnorm(-res$beta)
  A1 <- res$pf*A11
  A2 <- A12*(A11-1/sqrt(prod(1+(res$beta+1)*curvatures)))
  A3 <- (res$beta+1)*A12*(A11-Re(1/sqrt(prod(1+(res$beta+1i)*curvatures))));
  SORMpfTvedt <- A1+A2+A3;
  SORMbetaTvedt <- -stats::qnorm(SORMpfTvedt)


  cat("\n")
  duration<-proc.time()-tic

  # Return the results
  output <- list(
            "method"="SORM",
            "beta"=SORMbetaBreitung,
            "pf"=SORMpfBreitung,
            "FORM_pf"= res$pf,
            "FORMbeta"=res$beta,
            "design.point_X"=res$x_points,
            "SORMpfTvedt"=SORMpfTvedt,
            "SORMbetaTvedt"=SORMbetaTvedt,
            "SORMpfHB"=SORMpfHB,
            "SORMbetaHB"=SORMbetaHB,
            "runtime"=duration[1:5]
            )

  debug.print(debug.level,debug.TAG,c(duration), msg="SORM finished in [s]: ")
  return(output)
  }



