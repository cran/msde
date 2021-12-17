/// @file sdeUtils.h

#ifndef sdeUtils_h
#define sdeUtils_h 1

// utilities for sde's

#include "LinAlgUtils.h"
#include "mvnUtils.h"
//#include "sdeModel.h"

// --- xmvn, zmvn, llmvn -------------------------------------------------------

/// Calculate `x = sd * z + mean` where `sd` is specified by `sMod`.
///
/// Wrapper to pick between `xmvn_chol` or `xmvn_diag` depending on the diffusion specification of `sMod`.  See `mvnUtils.h`.
template <class sMod>
inline void xmvn(double *x, double *z,
		 double *mean, double *sd, int iFirst, int iLast) {
  if(sMod::diagDiff) {
    xmvn_diag(x, z, mean, sd, sMod::nDims, iFirst, iLast);
  } else {
    xmvn_chol(x, z, mean, sd, sMod::nDims, iFirst, iLast);
  }
  return;
}

/// Calculate `x = sd * z + mean` where `sd` is specified by `sMod`.
///
/// Simplified version with `iFirst = 0` and `iLast = sMod::nDims`.  Wrapper to pick between `xmvn_chol` or `xmvn_diag` depending on the diffusion specification of `sMod`.  See `mvnUtils.h`.
template <class sMod>
inline void xmvn(double *x, double *z, double *mean, double *sd) {
  if(sMod::diagDiff) {
    xmvn_diag(x, z, mean, sd, sMod::nDims);
  } else {
    xmvn_chol(x, z, mean, sd, sMod::nDims);
  }
  return;
}

/// Calculate `z = sd^{-1} * (x - mean)` where `sd` is specificed by `sMod`.
///
/// Wrapper to pick between `zmvn_chol` or `zmvn_diag` depending on the diffusion specification of `sMod`.  See `mvnUtils.h`.
template <class sMod>
inline void zmvn(double *z, double *x,
		      double *mean, double *sd, int iLast) {
  if(sMod::diagDiff) {
    zmvn_diag(z, x, mean, sd, iLast);
  } else {
    zmvn_chol(z, x, mean, sd, sMod::nDims, iLast);
  }
  return;
}

/// Calculate `z = sd^{-1} * (x - mean)` where `sd` is specificed by `sMod`.
///
/// Simplified version with `iLast = sMod::nDims`. Wrapper to pick between `zmvn_chol` or `zmvn_diag` depending on the diffusion specification of `sMod`.  See `mvnUtils.h`.
template <class sMod>
inline void zmvn(double *z, double *x,
		      double *mean, double *sd) {
  if(sMod::diagDiff) {
    zmvn_diag(z, x, mean, sd, sMod::nDims);
  } else {
    zmvn_chol(z, x, mean, sd, sMod::nDims);
  }
  return;
}

/// Calculate the log-density of `x ~ N(mean, sd * sd')` where `sd` is specified by `sMod`.
///
/// Wrapper to pick between `lmvn_chol` or `lmvn_diag` depending on the diffusion specification of `sMod`.  See `mvnUtils.h`.
template <class sMod>
inline double lmvn(double *x, double *z,
		   double *mean, double *sd, int iLast) {
  double ld;
  if(sMod::diagDiff) {
    ld = lmvn_diag(x, mean, sd, iLast);
  } else {
    ld = lmvn_chol(x, z, mean, sd, sMod::nDims, iLast);
  }
  return ld;
}

/// Calculate the log-density of `x ~ N(mean, sd * sd')` where `sd` is specified by `sMod`.
///
/// Simplified version assuming that `iLast = sMod::nDims`.  Wrapper to pick between `lmvn_chol` or `lmvn_diag` depending on the diffusion specification of `sMod`.  See `mvnUtils.h`.
template <class sMod>
inline double lmvn(double *x, double *z,
		   double *mean, double *sd) {
  double ld;
  if(sMod::diagDiff) {
    ld = lmvn_diag(x, mean, sd, sMod::nDims);
  } else {
    ld = lmvn_chol(x, z, mean, sd, sMod::nDims);
  }
  return ld;
}


/// Scale a diffusion matrix by the interobservation time.
///
/// If `sMod::sdDiff = true`, multiply `df` by `sqrtDT`.  Otherwise, multiply `df` by `sqrtDT^2`.
///
/// @param[in,out] df SDE diffusion matrix on the scale specified by `sMod::sdDiff` and `sMod::diagDiff`.
/// @param[in] sqrtDT Square-root of the interobservation time. 
template <class sMod>
inline void scaleDiff(double *df, double sqrtDT) {
  int ii;
  if(sMod::diagDiff) {
    // diagonal diffusion
    if(sMod::sdDiff) {
      // sd scale
      v_mult(df, sqrtDT, sMod::nDims);
    } else {
      // var scale
      for(ii=0; ii<sMod::nDims; ii++) {
	df[ii] = sqrt(df[ii])*sqrtDT;
      }
    }
  } else {
    // matrix diffusion
    if(!sMod::sdDiff) {
      // var scale
      chol_decomp(df, df, sMod::nDims);
    }
    U_mult(df, sqrtDT, sMod::nDims);
  }
  return;
}

/// Calculate the mean and variance of the Euler approximation.
///
/// @param[out] mean Array of length `n = sMod::nDims` giving the mean of the Euler approximation.
/// @param[out] sd Array of length `n^2`giving the (upper) Cholesky factor of the variance when `sMod::diagDiff = false`, or a vector of length `n` giving the standard deviations of each component when `sMod::diagDiff = true`.
/// @param[in] x Array of length `n` at which to evaluate the mean and variance.
/// @param[in] dT Interobservation time.
/// @param[in] sqrtDT Square-root of the interobservation time.
/// @param[in] theta Array of parameter values.
/// @param[in] sde Pointer to an `sMod` object which has methods to evaluate the SDE drift and diffusion functions.
template <class sMod>
inline void mvEuler(double *mean, double *sd,
		    double *x, double dT, double sqrtDT,
		    double *theta, sMod *sde) {
  // mean = x + drift(x,t,theta)*dT
  sde->sdeDr(mean, x, theta);
  v_mult(mean, dT, sMod::nDims);
  for(int ii = 0; ii < sMod::nDims; ii++) {
    mean[ii] += x[ii];
  }
  // sd = diff(x,t,theta)*sqrt(dT)
  sde->sdeDf(sd, x, theta);
  scaleDiff<sMod>(sd, sqrtDT);
  return;
}

#endif
