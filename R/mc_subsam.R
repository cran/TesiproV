#' @name MC_SubSam
#' @title Monte Carlo Simulation with Subset-Sampling
#' @description
#' MonteCarlo with Subset-Sampling (MC_SubSam) method for the
#' estimation of small failure probabilities in structural reliability analysis.
#'
#' Subset Simulation is an advanced Monte Carlo technique that
#' expresses a rare failure event as a sequence of intermediate
#' conditional events with higher probabilities.
#'
#' Failure is defined as \eqn{g(x) < 0}. The intermediate events are
#' \eqn{g(x) \le q_j} with the \code{p0}-quantile \eqn{q_j} of the limit-state
#' values of level \eqn{j}. Only the last level uses the failure event
#' \eqn{g(x) < 0}.
#'
#' The method combines conditional sampling with Markov Chain Monte Carlo (MCMC)
#' techniques to efficiently estimate very small probabilities of failure,
#' which would be computationally expensive to obtain using crude Monte Carlo simulation.
#'
#'
#' @param lsf limit-state function
#' @param lDistr list of basevariables in input space
#' @param Nsubset number of samples in each simulation level
#' @param p0 level probability or conditional probability
#' @param MaxSubsets maximum number of simulation levels that are used to terminate the simulation procedure to avoid infinite loop when the target domain cannot be reached
#' @param Alpha Significance level of the confidence interval. The interval
#'   covers the probability 1 - Alpha (default 0.05, i.e. 95 percent).
#' @param variance gaussian, uniform
#' @param debug.level If 0 no additional info if 2 high output during calculation
#' @param seed Optional integer value. If specified, sets a fixed seed for reproducible random numbers (useful for testing).
#'
#' @return MC_SubSam returns an object containing the following elements:
#' \itemize{
#'   \item \code{beta}: Estimated reliability index.
#'   \item \code{pf}: Estimated probability of failure.
#'   \item \code{betaCI}: Confidence interval of the reliability index.
#'   \item \code{pfCI}: Confidence interval of the probability of failure.
#'   \item \code{CoV}: Coefficient of variation of the estimator.
#'   \item \code{NumOfSubsets}: Number of simulation levels.
#'   \item \code{NumOfEvalLSF_nom}: Nominal number of limit-state function evaluations
#'         (Markov chains / iterations).
#'   \item \code{NumOfEvalLSF_eff}: Effective number of limit-state function evaluations.
#'   \item \code{runtime}: Total runtime of the function.
#' }
#'
#' @references
#' Au, S. K., & Beck, J. L. (2001).
#' Estimation of small failure probabilities in high dimensions
#' by subset simulation.
#' \emph{Probabilistic Engineering Mechanics}, 16(4), 263-277.
#'
#' Marelli, S., & Sudret, B. (2014).
#' UQLab: A framework for uncertainty quantification in Matlab.
#' In \emph{Proceedings of the 2nd International Conference on
#' Vulnerability, Risk Analysis and Management (ICVRAM2014)},
#' Liverpool, United Kingdom, 2554-2563.
#'
#' @author (C) 2021-2026 K. Nille-Hauf, T. Feiri, M. Ricker, T. Lux -- Hochschule Biberach (until 2022),
#' TU Dortmund University - Chair of Structural Concrete (since 2023)
#' @import pracma
#' @export

MC_SubSam <- function(lsf, # limit-state function
                      lDistr, # list of marginal distributions (basic variables)
                      Nsubset=1e5, # number of samples in each simulation level
                      p0=0.1, # level probability or conditional probability
                      MaxSubsets=10, # maximum number of simulation levels that are used to terminate the simulation procedure to avoid infinite loop when the target domain cannot be reached
                      Alpha=0.05, # confidence level
                      variance="uniform",
                      debug.level=0, # verbosity level
                      seed = NULL) # seed for random number generator
{

  # set seed for random number generator Mersenne Twister
  # (RNGkind() must come first, because it re-initialises the seed)
  RNGkind("Mersenne-Twister", "Inversion")
  if (!is.null(seed)) {
    set.seed(seed)
    if (debug.level >= 1) message(sprintf("Fixed RNG seed set to %d", as.integer(seed)))
  }

  debug.TAG <- "MC_SUS"
  debug.print(debug.level,debug.TAG,c(TRUE), msg="Monte-Carlo Simulation with Subset Sampling started...")
  tic<-proc.time()

  # ceck if lDistr is properly defined
  .check_lDistr(lDistr, "MC_SubSam")

  # === Transformation to the standard Gaussian space U and initial algorithm parameters definition ==

  # Isoprobabilistic transformation to the standard Gaussian space U
  Tlsf <- function(X) {
    XU <- stats::pnorm(X)
    for(i in seq_along(X)){
      X[i] <- lDistr[[i]][[1]]$q(XU[i])
    }
    return(lsf(X))
  }

  # Algorithm parameters definition (they are used for the "minimum" length of the Markov Chains)
  chainlength <- floor(1/p0)
  chainnumber <- floor(Nsubset*p0)
  n <- length(lDistr)

  # ====  SUBSET SIMULATION WITH MONTE CARLO SIMULATION (MCS) ====

  # === Initial MCS Stage: generate samples and read the function values ==
  ii <- 1
  U <- matrix(ncol=n,stats::rnorm(Nsubset*n,0,1))  # Generate samples from a standard normal distribution
  LSF <- matrix(apply(U,1,Tlsf))  # Compute the model responses for the samples (in standard normal space)
  LFScounter = nrow(U)  # Indicator to count the number of times the LSF is called
  PDFeval <- matrix(apply(stats::dnorm(U, 0, 1),1,prod)) # Compute the PDF values in the standard normal space
  q <- c()
  Pf <- c()
  Pfcond <- c()
  LSFhistory <- c()

  # === MCS stage: iterations among the subsets ==
  while (1) {
    Pf[ii] <- sum((LSF[]<0))/Nsubset   # Estimate the failure probability
    LSFhistory[[ii]] <- matrix(LSF)     # Store the history

    # Terminate the iteration due to convergence or maximum number of subsets
    if ((Pf[ii] > p0) || (ii >= MaxSubsets)) {
      Pfcond[ii] <- Pf[ii]
      q[ii] <- 0
      break
    }

    # Estimate the quantile
    q[ii] <- stats::quantile(LSF, p0, type=5)
    # Conditional probability of the next level. For continuous responses this is p0,
    # for discrete responses (ties at the quantile) it is the share of samples with LSF <= q
    Pfcond[ii] <- sum(LSF <= q[ii]) / Nsubset

    # Estimate the seeds for the next MCMC
    # Seeds: all samples with LSF <= q. If there are more than chainnumber of them (ties at the
    # quantile), a random subset is used, otherwise the samples at the quantile would be
    # under-represented in the seeds.
    idxsort <- sort(LSF,index.return=TRUE)$ix
    idx_le <- which(LSF <= q[ii])
    if (length(idx_le) > chainnumber) {
      idx_seed <- idx_le[sort(sample.int(length(idx_le), chainnumber))]
      idx_seed <- idx_seed[order(LSF[idx_seed])]
    } else {
      idx_seed <- idxsort[1:chainnumber]
    }
    Q <- U[idx_seed, , drop = FALSE]
    QPDF <- matrix(PDFeval[idx_seed])

    # Update the MCMC options
    AcceptCritGi <- q[ii]
    AcceptCritLSF <- matrix(LSF[idx_seed]) # Acceptance criterion for the next subset
    RuntimeLSF = AcceptCritLSF;

    Xseed <- Q
    Yseed <- matrix(apply(stats::dnorm(Xseed),1,prod))
    Xsample <- matrix(ncol=ncol(Xseed),nrow=0)
    Ysample <- matrix(ncol=ncol(Yseed),nrow=0)

    while (1) {

      if (variance == "gaussian") {
        Xprop = Xseed + (matrix(stats::rnorm(prod(dim(Xseed))),nrow(Xseed),ncol(Xseed))*1)
      } else { # if (variance == "uniform") {
        Xprop = Xseed + ((matrix(stats::runif(prod(dim(Xseed))),nrow(Xseed),ncol(Xseed))*2-1)*1)
      }

      a = dim(Xseed);
      b = dim(Xprop);
      c = dim(AcceptCritLSF);

      # Metropolis-Hastings MCMC component-wise

      UU <- Xseed  # repmat(Xseed,a(2),1);
      for (i in seq_len(a[2])[-1]) { UU <- rbind(UU,Xseed) }
      UPROP <- Xprop # repmat(Xprop, a(2), 1);
      for (i in seq_len(a[2])[-1]) { UPROP <- rbind(UPROP,Xprop) }
      jj = pracma::kron(pracma::eye(a[2]), pracma::ones(a[1],1));
      UU[jj==1] = UPROP[jj==1];
      YPROP <- apply(stats::dnorm(UU),1,prod) # prod(normpdf(UU),2);
      YY <- Yseed # repmat(Yseed, a(2), 1);
      for (i in seq_len(a[2])[-1]) { YY <- rbind(YY,Yseed) }
      idx <- (matrix(stats::runif(b[1]*b[2])) < YPROP/YY); # rand(b(1)*b(2),1)
      IDX <- matrix(nrow=a[1], data=idx)

      # Write an updated proposed candidate sample
      Uprop2 <- Xseed;
      Uprop2[IDX] = Xprop[IDX];
      Yprop <- apply(stats::dnorm(Uprop2),1,prod)
      Xprop <- Uprop2;

      idx <- rowSums(IDX)
      idx[rowSums(IDX) != 0] = 1;

      # When X is in a range of PDF=0
      idx[Yseed == 0] = 1;

      # Assigning the new samples of the Markov Chain
      LSFacc = NaN*zeros(a[1], c[2]);
      for (i in 1:a[1]) {
        if (idx[i]) {
          LSFacc[i] <- Tlsf(Xprop[i,])
          LFScounter = LFScounter + 1
        }
      }

      for (i in 1:nrow(Xseed)) {
        if (idx[i] && LSFacc[i]<=AcceptCritGi) {
          Xseed[i,] <- Uprop2[i,]
        }
      }
      for (i in 1:nrow(Yseed)) {
        if (idx[i] && LSFacc[i]<=AcceptCritGi) {
          Yseed[i] <- Yprop[i]
        }
      }

      # Assign the new samples according to the threshold
      for (i in 1:nrow(AcceptCritLSF)) {
        if (idx[i] && LSFacc[i]<=AcceptCritGi) {
          AcceptCritLSF[i,] <- LSFacc[i,]
        }
      }
      RuntimeLSF <- rbind(RuntimeLSF,AcceptCritLSF)

      # Define the new sample vector
      Xsample <- rbind(Xsample, Xseed) # assign LSF to the mcmcopts to have it somewhere here
      Ysample <- rbind(Ysample, Yseed) # Also store the (calculated) limit-state function evaluations
      if (nrow(Xsample) >= Nsubset - chainnumber) {break;} #Stopping criterion

      loc_pf <- prod(Pfcond[seq_len(ii - 1)])*Pf[length(Pf)] # Probability of failure
      loc_beta <- -stats::qnorm(loc_pf)
      info.print(debug.TAG,debug.level,c("ii","Beta","pf"),c(ii,loc_beta, loc_pf))
    }

    # Select Nsubset - nrow(Q) chain states. Complete Markov chain steps are used as a whole.
    # If the last step is only needed partially, its chains are drawn at random, because the
    # chains are sorted by their seed and taking the first ones would bias the sample.
    n_need <- Nsubset - nrow(Q)
    n_full <- (n_need %/% chainnumber) * chainnumber
    sel <- seq_len(n_full)
    if (n_need > n_full) {
      sel <- c(sel, n_full + sort(sample.int(chainnumber, n_need - n_full)))
    }
    U <- rbind(Q, Xsample[sel, , drop = FALSE])
    PDFeval <- rbind(QPDF, matrix(ncol=1, Ysample[sel, ]))
    LSF <- RuntimeLSF[c(seq_len(chainnumber), chainnumber + sel), ]
    ii <- ii + 1 #  Increase the counter
  }


  # === Final MCS stage: collect the results ==
  ResPf <- prod(Pfcond) # Probability of failure (product of the conditional probabilities)
  ResBeta <- -stats::qnorm(ResPf) # Reliability index
  debug.print(debug.level,debug.TAG,ResPf,"Results Pf:")
  debug.print(debug.level,debug.TAG,ResBeta,"Results Beta:")


  # Coefficient of Variation
  delta2 <- c((1-Pfcond[1])/Pfcond[1]/Nsubset)

  # Computation of gamma
  for (jj in 1:ii) {
    # rearrange the LSFhistory vector to a matrix in order to have each chain separate
    LSFH = matrix(nrow=chainnumber, data=LSFhistory[[jj]][1:(chainlength * chainnumber)])
    Indicator = pracma::zeros(nrow(LSFH),ncol(LSFH));
    if (jj == ii) {
      Indicator[LSFH < 0] = 1;        # last level: failure event g < 0
    } else {
      Indicator[LSFH <= q[jj]] = 1;   # intermediate level
    }

    # Compute the rho (i.e. the auto-correlation between the Markov Chain samples in the i-th subset)
    # and consider every chain independently
    Rik <- vector() #= [];
    for (k in 1:chainlength - 1) {
      sums <- 0;
      for (l in 1 : ((Nsubset/chainnumber)-k)) {
        sums <- sums + sum(Indicator[,l] * Indicator[,l+k]);
      }
      Rik[k] <- sums * 1/(Nsubset - k * chainnumber) - Pfcond[jj]^2;
    }
    # Compute gamma for each each iteration
    rho <- Rik / (Pfcond[jj]*(1-Pfcond[jj]));
    gamma <- 2 * sum( rho*(1-(1:length(rho))*chainnumber / Nsubset) );

    delta2[jj] = (1-Pfcond[jj])/Pfcond[jj]/Nsubset*(1+gamma)
  }

  debug.print(debug.level,debug.TAG,ii,"Number of subsets")
  debug.print(debug.level,debug.TAG,ii*Nsubset,"Number of model evaluations (nominally)")
  debug.print(debug.level,debug.TAG,LFScounter,"Number of model evaluations (effectively)")


  # Compute the confidence bounds
  CoV <- sqrt(sum(delta2))
  ResPfCI <- ResPf * matrix(c(1 + stats::qnorm(Alpha/2)*CoV, 1 + stats::qnorm(1-Alpha/2)*CoV), nrow=1, byrow=TRUE)

  debug.print(debug.level,debug.TAG,CoV,"CoV: ")
  debug.print(debug.level,debug.TAG,ResPfCI,"ResPfCI: ")
  debug.print(debug.level,debug.TAG,t(apply(-stats::qnorm(ResPfCI),1,rev)),"ResBetaCI: ")


  duration<-proc.time()-tic

  output<-list(
    "method"="MC_SubSam",
    "beta"=ResBeta,
    "pf"=ResPf,
    "betaCI"=t(apply(-stats::qnorm(ResPfCI),1,rev)),
    "pfCI"=ResPfCI,
    "CoV"=CoV,
    "NumOfSubsets"=ii,
    "NumOfEvalLSF_nom"=ii*Nsubset,
    "NumOfEvalLSF_eff"=LFScounter,
    "runtime"=duration[1:5]
  )
  info.print(debug.TAG,debug.level,c("Beta","pf","CoV","NumOfEvalLSF_eff"),c(ResBeta,ResPf,CoV,LFScounter))

  debug.print(debug.level,debug.TAG,c(duration), msg="Monte-Carlo Simulation with SubsetSampling finished in: ")
  return(output)
}

# ====  END OF SUBSET SIMULATION WITH MCS ====
