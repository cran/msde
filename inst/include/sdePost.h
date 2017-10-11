#ifndef sdePost_h
#define sdePost_h

#include <Rcpp.h>
//using namespace Rcpp;
typedef Rcpp::LogicalVector Logical;
typedef Rcpp::NumericVector Numeric;
typedef Rcpp::IntegerVector Integer;
typedef Rcpp::List List;
//[[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include "sdeMCMC.h"
#include "mcmcUtils.h"
#include "sdeInterface.h"

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
  int nMiss0 = nDims-nDimsPerObs[0]; // unobserved states in first observation
  int nMissN = nDims-nDimsPerObs[nComp-1]; // unobserved states in last observation

  // output variables
  Numeric paramsOut(nParamsOut);
  Numeric dataOut(nDataOut);
  Integer paramAcceptOut(nParams + nMiss0);
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
  delete [] tunePar_adapt;

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
