/// @file sdePost.h

#ifndef sdePost_h
#define sdePost_h

#include <Rcpp.h>
using namespace Rcpp;
typedef Rcpp::LogicalVector Logical;
typedef Rcpp::NumericVector Numeric;
typedef Rcpp::IntegerVector Integer;
typedef Rcpp::List List;
#include <progress.hpp>
#include <progress_bar.hpp>
#include "sdeMCMC.h"
#include "mcmcUtils.h"
#include "sdeRobj.h"

/// The MCMC algorithm employed is a Metropolis-within-Gibbs sampler on the SDE variables at each time point and each parameter.  The proposal for the the SDE variable at time `t` is the Brownian bridge proposal of Eraker (2002).  The proposal for the parameters (and the unknown SDE variables at time `t = 0`) is an adaptive random walk.
///
/// @param[in] initParams Parameter vector of length `nParams` to initialize the sampler.
/// @param[in] initData SDE vector of length `nDims * nComp` to initialize the sampler.
/// @param[in] dT Vector of `nComp-1` interobservation times.
/// @param[in] nDimsPerObs Integer vector of length `nComp` with elements between 0 and `nDims` specifying how many of the corresponding SDE variables are observed (as opposed to missing).  So, `nDimsPerObs[t] = j` means that the first `j` element of the SDE at time `t = 0` are observed and the other `nDims-j` are unobserved.
/// @param[in] fixedParamsIn Logical vector of length `nParams` indicating which parameters are fixed (i.e., known) and which are unknown (i.e., will be sampled by the MCMC algorithm).
/// @param[in] nSamples Number of MCMC samples.
/// @param[in] burn Number of burn-in samples.
/// @param[in] nParamsOut The size of the parameter output, i.e., `nParams * nSamples`.
/// @param[in] nDataOut The size of the data output, which is determined by `dataOutSmp`, `dataOutComp`, and `dataOutDims`.
/// @param[in] dataOutSmp Integer vector with elements between 0 and `nSamples-1` indicating which MCMC iterations (after burn-in) of the SDE variables to keep.
/// @param[in] dataOutComp Integer vector with elements between 0 and `nComp-1` indicating which SDE time points to keep.
/// @param[in] dataOutDims Integer vector with elements between 0 and `nDims` indicating which SDE variables to keep.
/// @param[in] updateParams Scalar specifying when the parameters get updated.  See `updateComponent` in `mcmcUtils.h`.
/// @param[in] updateData Scalar specifying when the SDE variables get updated.
/// @param[in] priorArgsIn A list of vectors specifying the prior.
/// @param[in] tunePar List specifying the adaptive tuning of the MWG algorithm for the parameter updates.  See `mwgAdapt` in `mcmcUtils.h`.
/// @param[in] updateLogLik Whether or not to return the complete data loglikelihood at every iteration.
/// @param[in] nLogLikOut  The size of the loglikelihood output (`nSamples` or 1).
/// @param[in] updateLastMiss Whether or not to return all samples of the missing SDE variables at the last time point.  This is useful for Bayesian forecasting of the SDE trajectory.
/// @param[in] nLastMissOut The size of the `LastMiss` output (`nSamples * (nDims - nDimsPerObs[nComp-1])` or 1).
/// @param[in] nCores Number of parallel cores to use.  Has no effect if `OpenMP` is not enabled.
/// @param[in] displayProgress Whether or not to display a progress bar.
///
/// @return A list with elements:
/// - `paramsOut`: The vector of parameter MCMC draws.
/// - `dataOut`: The vector of SDE variable MCMC draws.
/// - `paramAccept`: A vector of acceptance rates for the parameters.
/// - `gibbsAccept`: A vector of acceptance rates for the SDE variables.
/// - `logLikOut`: The loglikelihood output.
/// - `lastMissOut`: The SDE missing variables for the last time point.
/// - `lastIter`: The parameter and SDE variables for the last MCMC iteration.  This is useful for checkpointing the MCMC.
/// - `mwgSd`: The adapted jump sizes for the parameter (and initial SDE variables) MWG.
template <class sMod, class sPi>
  inline List sdeRobj<sMod, sPi>::Post(Numeric initParams,
				       Numeric initData,
				       Numeric dT,
				       Integer nDimsPerObs,
				       Logical fixedParamsIn,
				       int nSamples, int burn,
				       int nParamsOut, int nDataOut,
				       Integer dataOutSmp,
				       Integer dataOutComp,
				       Integer dataOutDims,
				       double updateParams,
				       double updateData,
				       List priorArgsIn, List tunePar,
				       int updateLogLik,
				       int nLogLikOut,
				       int updateLastMiss,
				       int nLastMissOut, int nCores,
				       bool displayProgress) {
  RNGScope scope;
  int ii, jj, kk;

  // problem dimensions
  int nDims = sMod::nDims;
  int nParams = sMod::nParams;
  int nComp = initData.length()/nDims;
  int nDimsOut = dataOutDims.length();
  int nCompOut = dataOutComp.length();
  // int nMiss0 = nDims-nDimsPerObs[0]; // unobserved states in first observation
  int nMissN = nDims-nDimsPerObs[nComp-1]; // unobserved states in last observation

  // output variables
  Numeric paramsOut(nParamsOut);
  Numeric dataOut(nDataOut);
  // Integer paramAcceptOut(nParams + nMiss0);
  Integer paramAcceptOut(nParams + nDims);
  Integer gibbsAcceptOut(nComp);
  Numeric logLikOut(nLogLikOut);
  Numeric lastMissOut(nLastMissOut);
  Numeric lastIter(nParams + nComp*nDims);
  Numeric mwgSdOut(nParams + nDims);
  // pointers to acceptance rate counters for internal use
  int *paramAccept = INTEGER(paramAcceptOut);
  int *gibbsAccept = INTEGER(gibbsAcceptOut);
  double *mwgSd = REAL(mwgSdOut);
  // convert LogicalVectors to vector of bools
  bool *fixedParams = new bool[nParams];
  bool *tunePar_adapt = new bool[nParams + nDims];
  convert_Logical(fixedParams, fixedParamsIn);
  convert_Logical(tunePar_adapt, tunePar["adapt"]);

  // MCMC tuning parameters
  for(ii=0; ii<nParams+nDims; ii++) {
    mwgSd[ii] = REAL(tunePar["sd"])[ii];
  }
  mwgAdapt tuneMCMC(REAL(tunePar["max"]), REAL(tunePar["rate"]),
		    tunePar_adapt, nParams+nDims);
		    // INTEGER(tunePar["adapt"]), nParams+nDims);

  // prior specification
  // hyper parameters: actual prior gets constructed inside MCMC object
  PriorArgs priorArgs(priorArgsIn);
  /* int nArgs = priorArgs.length(); */
  /* double **phi = new double*[nArgs]; */
  /* int *nEachArg = new int[nArgs]; */
  /* for(ii=0; ii<nArgs; ii++) { */
  /*   if(Rf_isNull(priorArgs[ii])) { */
  /*     nEachArg[ii] = 0; */
  /*   } else { */
  /*     nEachArg[ii] = as<NumericVector>(priorArgs[ii]).length(); */
  /*     phi[ii] = REAL(priorArgs[ii]); */
  /*   } */
  /* } */

  // initialize MCMC
  // prior gets constructed inside of object -- is this really beneficial?
  sdeMCMC<sMod,sPi> mcmc(nComp, REAL(dT), REAL(initData),
			 REAL(initParams),
			 INTEGER(nDimsPerObs), fixedParams,
			 priorArgs.phi, priorArgs.nArgs, priorArgs.nEachArg,
			 nCores);

  // progress bar
  Progress Progress_Bar(burn + nSamples, displayProgress);

  // main MCMC loop
  jj = 0;
  for(int smp = -burn; smp < nSamples; smp++) {
    // user interrupt
    if(smp % (int) 5e3) {
      Rcpp::checkUserInterrupt();
      Progress_Bar.increment();
    }
    // missing data update
    if(updateComponent(updateData, smp)) {
      mcmc.missGibbsUpdate(mwgSd, gibbsAccept, paramAccept);
    }
    // parameter update
    if(updateComponent(updateParams, smp)) {
      mcmc.paramVanillaUpdate(mwgSd, paramAccept);
    }
    // adaptive MCMC
    tuneMCMC.adapt(mwgSd, paramAccept, burn+smp+1);
    // log-likelihood
    // TODO: keep track of this interally after every MCMC step
    if(smp >= 0) {
      if(updateLogLik) logLikOut[smp] = mcmc.loglik(mcmc.currTheta, mcmc.currX);
    }
    // storage
    if(smp == dataOutSmp[jj]) {
      if(updateData > 0.0) {
	for(ii=0; ii<nCompOut; ii++) {
	  //Rprintf("dataOutComp[%i] = %i\n", ii, dataOutComp[ii]);
	  for(kk=0; kk<nDimsOut; kk++) {
	    // Rprintf("dataOut[%i,%i,%i] = currX[%i,%i] = \n",
	    // 	    kk, ii, jj, dataOutDims[kk], dataOutComp[ii]);
	    dataOut[jj*nDimsOut*nCompOut+ii*nDimsOut+kk] = mcmc.currX[dataOutComp[ii]*nDims+dataOutDims[kk]];
	  }
	}
      }
      jj++;
    }
    if((updateParams > 0.0) && (smp >= 0)) {
      for(ii=0; ii<nParams; ii++) paramsOut[smp*nParams + ii] = mcmc.currTheta[ii];
    }
    if(updateLastMiss && smp >= 0) {
      for(ii=0; ii<nMissN; ii++) {
	lastMissOut[smp*nMissN+ii] = mcmc.currX[(nComp-1)*nDims+nDimsPerObs[nComp-1]+ii];
      }
    }
  }

  // store last iteration, to resume MCMC later if needed
  for(ii=0; ii<nParams; ii++) {
    lastIter[ii] = mcmc.currTheta[ii];
  }
  for(ii=0; ii<nDims*nComp; ii++) {
    lastIter[nParams+ii] = mcmc.currX[ii];
  }

  // delete dynamic variables
  delete [] fixedParams;
  // delete [] tunePar_adapt;

  return List::create(_["paramsOut"] = paramsOut,
		      _["dataOut"] = dataOut,
		      _["paramAccept"] = paramAcceptOut,
		      _["gibbsAccept"] = gibbsAcceptOut,
		      _["logLikOut"] = logLikOut,
		      _["lastMissOut"] = lastMissOut,
		      _["lastIter"] = lastIter,
		      _["mwgSd"] = mwgSdOut);
}

#endif
