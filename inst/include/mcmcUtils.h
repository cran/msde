/// @file mcmcUtils.h

#ifndef mcmcUtils_h
#define mcmcUtils_h 1

//#include <Rmath.h>
//using namespace Rcpp;
#include "rngUtils.h"

/// Determine whether or not to update a particular component of a Gibbs sampler.
/// 
/// `x` is either:
/// - `true` or `false`: in which case the component always/never gets updated.
/// - an integer, in which case the component gets updated every `x` iterations (i.e., when `iter` is a multiple of `x`).
/// - A double between 0 and 1, in which case the component is randomly updated with probability `x`.
///
/// @param[in] x Scalar determining how to update the component.
/// @param[in] iter Counter for the current MCMC iteration.
///
/// @return Whether or not to update the component.
inline bool updateComponent(double x, int iter) {
  bool update = false;
  if(x > 0.0) {
    if(x < 1.0) {
      update = sdeRNG::runif() <= x;
    }
    else {
      update = (iter % (int) x) == 0;
    }
  }
  return update;
}

/// Metropolis-within-Gibbs adaptive sampler.
///
/// Given a vector of jump sizes, increase or decrease each element depending on whether the cumulative acceptance rate is above or below 0.44.  The amount of change on log-scale is
/// ```
/// delta = min(adaptMax, 1/nIter^adaptRate)
/// ```
class mwgAdapt {
 private:
  double *adaptMax, *adaptRate;
  bool *doAdapt;
  int nRV;
 public:
  /// Class constructor.
  mwgAdapt(double *amax, double *arate, bool *adapt, int n);
  /// Class destructor.
  ~mwgAdapt();
  /// Update the MWG jump sizes.
  void adapt(double *mwgSd, int *nAccept, int nIter);
};

/// @param[in] amax Vector of `n` values of `adaptMax` in formula.
/// @param[in] arate Vector of `n` values of `adaptRate` in formula.
/// @param[in] adapt Vector indicating which of the `n` components are to be adaptively updated.
/// @param[in] n Total number of MWG components.
inline mwgAdapt::mwgAdapt(double *amax, double *arate,
			  bool *adapt, int n) {
  nRV = n;
  adaptMax = new double[nRV];
  adaptRate = new double[nRV];
  doAdapt = new bool[nRV];
  for(int ii=0; ii<nRV; ii++) {
    adaptMax[ii] = amax[ii];
    adaptRate[ii] = arate[ii];
    doAdapt[ii] = adapt[ii];
  }
}

// mwgAdapt(double *amax, double *arate, int *adapt, int n) {
//   nRV = n;
//   adaptMax = new double[nRV];
//   adaptRate = new double[nRV];
//   doAdapt = new bool[nRV];
//   for(int ii=0; ii<nRV; ii++) {
//     adaptMax[ii] = amax[ii];
//     adaptRate[ii] = arate[ii];
//     doAdapt[ii] = (adapt[ii] != 0);
//   }
// }

inline mwgAdapt::~mwgAdapt() {
  delete [] adaptMax;
  delete [] adaptRate;
  delete [] doAdapt;
}

/// @param[in,out] mwgSd MWG jump sizes to be adapted (update performed in-place).
/// @param[in] nAccept Number of accepted MWG draws so far.
/// @param[in] nIter Total number of MWG draws so far.
inline void mwgAdapt::adapt(double *mwgSd, int *nAccept, int nIter) {
  double acc;
  double lsig;
  double delta;
  const double targAcc = 0.44;
  for(int ii=0; ii<nRV; ii++) {
    if(doAdapt[ii]) {
      acc = (double) nAccept[ii] / (double) nIter;
      delta = pow((double) nIter, -adaptRate[ii]);
      if(delta > adaptMax[ii]) delta = adaptMax[ii];
      lsig = log(mwgSd[ii]);
      lsig += acc < targAcc ? -delta : delta;
      mwgSd[ii] = exp(lsig);
    }
  }
  return;
}


#endif
