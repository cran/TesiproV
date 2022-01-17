# Reliability Analysis at Biberach University of applied sciences
# Subset Simulation for reliability analysis
# (c) adapted from: S. Marelli, and B. Sudret, UQLab: A framework for uncertainty quantification in Matlab, Proc. 2nd Int. Conf. on Vulnerability, Risk Analysis and Management (ICVRAM2014), Liverpool (United Kingdom), 2014, 2554-2563.
#' @title MonteCarlo with Subset-Sampling
#'
#' @name MC_SubSam
#'
#' @param lsf limit-state function
#' @param lDistr list of basevariables in input space
#' @param Nsubset number of samples in each simulation level
#' @param p0 level probability or conditional probability
#' @param MaxSubsets maximum number of simulation levels that are used to terminate the simulation procedure to avoid infinite loop when the target domain cannot be reached
#' @param Alpha  confidence level
#' @param variance gaussian, uniform
#' @param debug.level If 0 no additional info if 2 high output during calculation
#'
#' @return The results are provided within a list() of the following elements:
#' @return beta
#' @return pf
#' @return betaCI and pfCI are the corresponding confidence intervals
#' @return CoV COV of the result
#' @return NumOfSubsets Amount of Markov-Chains
#' @return NumOfEvalLSF_nom Markov-Chains times Iterations
#' @return NumOfEvalLSF_eff Internal counter that shows the real evaluations of the lsf
#' @return runtime Duration since start to finish of the function
#'
#' @references AU, S. K. & BECK, J. L. Estimation of small failure probabilities in high dimensions by subset simulation. Probabilistic Engineering Mechanics, 2001, 16.4: 263-277.
#'
#' @author (C) 2021 - K. Nille-Hauf, T. Feiri, M. Ricker - Hochschule Biberach, Institut fuer Konstruktiven Ingenieurbau
#' @import pracma
#' @export

MC_SubSam <- function(lsf, # limit-state function
                        lDistr, # dimension of the input space
                        Nsubset=1e5, # number of samples in each simulation level
                        p0=0.1, # level probability or conditional probability
                        MaxSubsets=10, # maximum number of simulation levels that are used to terminate the simulation procedure to avoid infinite loop when the target domain cannot be reached
                        Alpha=0.05, # confidence level
                        variance="uniform",
                        debug.level=0) # gaussian, uniform
{

  debug.TAG <- "MC_SUS"
  debug.print(debug.level,debug.TAG,c(TRUE), msg="Monte-Carlo Simulation with Subset Sampling started...")
  tic<-Sys.time()
  # === Transformation to the standard Gaussian space U and initial algorithm parameters definition ==

  # Isoprobabilistic transformation to the standard Gaussian space U
  Tlsf <- function(X) {
    XU <- stats::pnorm(X)
    for(i in 1:length(X)){
      X[i] <- lDistr[[i]][[1]]$q(XU[i, drop = FALSE])
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
    Pf[ii] <- sum((LSF[]<=0))/Nsubset   # Estimate the failure probability
    LSFhistory[[ii]] <- matrix(LSF)     # Store the history

    # Terminate the iteration due to convergence or maximum number of subsets
    if ((Pf[ii] > p0) || (ii >= MaxSubsets)) {
      Pfcond[ii] <- Pf[ii]
      q[ii] <- 0
      break
    }

        # Estimate the quantile
    q[ii] <- stats::quantile(LSF, p0, type=5)
    Pfcond[ii] <- p0

    # Estimate the seeds for the next MCMC
    idxsort <- sort(LSF,index.return=TRUE)$ix
    Q <- U[idxsort[1:chainnumber],]
    QPDF <- matrix(PDFeval[idxsort[1:chainnumber]])

    # Update the MCMC options
    AcceptCritGi <- q[ii]
    AcceptCritLSF <- matrix(LSF[idxsort[1:chainnumber]]) # Acceptance criterion for the next subset
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

      a = size(Xseed);
      b = size(Xprop);
      c = size(AcceptCritLSF);

      # Metropolies-Hasting MCMC component-wise

      UU <- Xseed  # repmat(Xseed,a(2),1);
      for (i in 2:a[2]) { UU <- rbind(UU,Xseed) }
      UPROP <- Xprop # repmat(Xprop, a(2), 1);
      for (i in 2:a[2]) { UPROP <- rbind(UPROP,Xprop) }
      jj = pracma::kron(pracma::eye(a[2]), pracma::ones(a[1],1));
      UU[jj==1] = UPROP[jj==1];
      YPROP <- apply(stats::dnorm(UU),1,prod) # prod(normpdf(UU),2);
      YY <- Yseed # repmat(Yseed, a(2), 1);
      for (i in 2:a[2]) { YY <- rbind(YY,Yseed) }
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
        if (idx[i] && LSFacc[i]<AcceptCritGi) {
          Xseed[i,] <- Uprop2[i,]
        }
      }
      for (i in 1:nrow(Yseed)) {
        if (idx[i] && LSFacc[i]<AcceptCritGi) {
          Yseed[i] <- Yprop[i]
        }
      }

      # Assign the new samples according to the threshold
      for (i in 1:nrow(AcceptCritLSF)) {
        if (idx[i] && LSFacc[i]<AcceptCritGi) {
          AcceptCritLSF[i,] <- LSFacc[i,]
        }
      }
      RuntimeLSF <- rbind(RuntimeLSF,AcceptCritLSF)

      # Define the new sample vector
      Xsample <- rbind(Xsample, Xseed) # assign LSF to the mcmcopts to have it somewhere here
      Ysample <- rbind(Ysample, Yseed) # Also store the (calculated) limit-state function evaluations
      if (nrow(Xsample) >= round((1 - p0)*Nsubset)) {break;} #Stopping criterion

      loc_pf <- p0^(ii-1)*Pf[length(Pf)] # Probability of failure
      loc_beta <- -stats::qnorm(loc_pf)
      info.print(debug.TAG,c("ii","Beta","pf"),c(ii,loc_beta, loc_pf))
    }

    U <- rbind(Q,Xsample[1:(Nsubset-nrow(Q)),])
    PDFeval <- rbind(QPDF,matrix(ncol=1,Ysample[1:(Nsubset-nrow(Q)),]))
    LSF <- RuntimeLSF[1:Nsubset,]
    ii <- ii + 1 #  Increase the counter
  }


  # === Final MCS stage: collect the results ==
  ResPf <- p0^(ii-1)*Pf[length(Pf)] # Probability of failure
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
    Indicator[LSFH <= q[jj]] = 1;

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


  cat("\n")
  duration<-Sys.time()-tic

  output<-list(
    "method"="MCS",
    "beta"=ResBeta,
    "pf"=ResPf,
    "betaCI"=t(apply(-stats::qnorm(ResPfCI),1,rev)),
    "pfCI"=ResPfCI,
    "CoV"=CoV,
    "NumOfSubsets"=ii,
    "NumOfEvalLSF_nom"=ii*Nsubset,
    "NumOfEvalLSF_eff"=LFScounter,
    "runtime"=duration
  )
  info.print(debug.TAG,debug.level,c("Beta","pf","CoV","NumOfEvalLSF_eff"),c(ResBeta,ResPf,CoV,LFScounter))

  debug.print(debug.level,debug.TAG,c(duration), msg="Monte-Carlo Simulation with SubsetSampling finished in: ")
  return(output)
}

# ====  END OF SUBSET SIMULATION WITH MCS ====
