#ifndef sdeUtils_h
#define sdeUtils_h 1

// utilities for sde's

#include "LinAlgUtils.h"
#include "mvnUtils.h"
//#include "sdeModel.h"

// xmvn, zmvn, llmvn
template <class sMod>
inline void xmvn(double *x, double *z,
		 double *mean, double *sd, int n) {
  if(sMod::diagDiff) {
    xmvn_diag(x, z, mean, sd, n);
  } else {
    xmvn_chol(x, z, mean, sd, n);
  }
  return;
}

template <class sMod>
inline void zmvn(double *z, double *x,
		      double *mean, double *sd, int n, int nMax) {
  if(sMod::diagDiff) {
    zmvn_diag(z, x, mean, sd, n, nMax);
  } else {
    zmvn_chol(z, x, mean, sd, n, nMax);
  }
  return;
}

template <class sMod>
inline double lmvn(double *x, double *z,
		   double *mean, double *sd, int n) {
  double ld;
  if(sMod::diagDiff) {
    ld = lmvn_diag(x, z, mean, sd, n);
  } else {
    ld = lmvn_chol(x, z, mean, sd, n);
  }
  return ld;
}


// rescale diffusion depending on whether it's on the sd or variance scale
// and diagonal or not.
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

// Euler approximation mean and variance
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
