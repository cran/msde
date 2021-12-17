/// @file sdeSim.h

#ifndef sdeSim_h
#define sdeSim_h

#include <Rcpp.h>
// //[[Rcpp::depends("RcppArmadillo")]]
// #include <RcppArmadillo.h>
typedef Rcpp::LogicalVector Logical;
typedef Rcpp::NumericVector Numeric;
typedef Rcpp::IntegerVector Integer;
typedef Rcpp::List List;
//using namespace Rcpp;
#include "rngUtils.h"
#include "sdeUtils.h"
#include "sdeRobj.h"

/// @param[in] nDataOut Size of final data output.
/// @param[in] N length of each SDE trajectory.
/// @param[in] burn Number of initial time points to discard as burn-in.
/// @param[in] reps Number of trajectories to simulate.
/// @param[in] r Number of observations to discard between the ones that we save.  So for example, if `r=1` we keep all observations, if `r=10` we keep every 10.
/// @param[in] dT Effective interobservation time.
/// @param[in] MAXBAD Maximum number of invalid Euler-Maruyama steps before method terminates.
/// @param[in] initData Vector of length `nDims` or `reps * nDims` giving the starting value of each trajectory.
/// @param[in] params Vector of length `nParams` or `reps * nParams` giving the parameters for each trajectory.
/// @param[in] singleX Whether or not `initData` is of length `nDims`.
/// @param[in] singleTheta Whether or not `params` is of length `nParams`.
///
/// @return An `Rcpp::List` with elements `dataOut`, consisting of `reps * N * nDims` SDE observations, and `nBadDraws`, the number of bad Euler-Maruyama steps (total).  Note that `initData` is not included in the output.
template <class sMod, class sPi>
  inline List sdeRobj<sMod, sPi>::Sim(int nDataOut,
				      int N, int burn, int reps, int r,
				      double dT, int MAXBAD,
				      Numeric initData,
				      Numeric params,
				      bool singleX, bool singleTheta) {
  RNGScope scope;

  int nDims = sMod::nDims;
  int nParams = sMod::nParams;
  double sqrtDT = sqrt(dT);
  int bad = 0;
  // output
  Numeric dataOut(nDataOut);
  int nBadDraws;

  // storage
  sMod *sde = new sMod; // sde model
  double *mean = new double[nDims]; // mean
  double *sd = new double[nDims*nDims]; // cholesky factor
  double *X = new double[nDims]; // current value
  double *tmpX = new double[nDims]; // proposed value
  double *Z = new double[nDims]; // random draw
  double *theta; // pointer to parameter
  int ii,jj,kk;

  for(ii = 0; ii < reps; ii++) {
    // initialize chains
    for(kk = 0; kk < nDims; kk++) {
      X[kk] = initData[ii*(!singleX)*nDims + kk];
    }
    theta = &params[ii*(!singleTheta)*nParams];
    // loop through
    for(jj = -burn*r; jj < N*r; jj++) {
      // repeatedly draw from Euler until proposal is valid
      mvEuler<sMod>(mean, sd, X, dT, sqrtDT, theta, sde);
      do {
	for(kk = 0; kk < nDims; kk++) {
	  Z[kk] = sdeRNG::rnorm();
	}
	xmvn<sMod>(X, Z, mean, sd);
	// validate draw
      } while(!sde->isValidData(X, theta) && bad++ < MAXBAD);
      if (bad == MAXBAD) {
	goto stop;
      }
      // store
      if(jj >= 0 && (jj+1) % r == 0) {
	for(kk = 0; kk < nDims; kk++) {
	  dataOut[ii*N*nDims + (jj/r)*nDims + kk] = X[kk];
	}
      }
    }
  }

 stop:
  nBadDraws = bad;

  delete [] X;
  delete [] tmpX;
  delete [] Z;
  delete [] mean;
  delete [] sd;
  delete sde;

  return List::create(_["dataOut"] = dataOut, _["nBadDraws"] = nBadDraws);
}

#endif
