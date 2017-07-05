#ifndef sdeUtils_h
#define sdeUtils_h 1

//[[Rcpp::depends("msde")]]

// utilities for sde's

#include <LinAlgUtils.h>
#include <mvnUtils.h>
#include "sdeModel.h"

// xmvn, zmvn, llmvn
inline void xmvn(double *x, double *z,
		 double *mean, double *sd, int n) {
  if(sdeModel::diagDiff) {
    xmvn_diag(x, z, mean, sd, n);
  } else {
    xmvn_chol(x, z, mean, sd, n);
  }
  return;
}

inline void zmvn(double *z, double *x,
		      double *mean, double *sd, int n, int nMax) {
  if(sdeModel::diagDiff) {
    zmvn_diag(z, x, mean, sd, n, nMax);
  } else {
    zmvn_chol(z, x, mean, sd, n, nMax);
  }
  return;
}

inline double lmvn(double *x, double *z,
		   double *mean, double *sd, int n) {
  double ld;
  if(sdeModel::diagDiff) {
    ld = lmvn_diag(x, z, mean, sd, n);
  } else {
    ld = lmvn_chol(x, z, mean, sd, n);
  }
  return ld;
}


// rescale diffusion depending on whether it's on the sd or variance scale
// and diagonal or not.
inline void scaleDiff(double *df, double sqrtDT) {
  int ii;
  if(sdeModel::diagDiff) {
    // diagonal diffusion
    if(sdeModel::sdDiff) {
      // sd scale
      v_mult(df, sqrtDT, sdeModel::nDims);
    } else {
      // var scale
      for(ii=0; ii<sdeModel::nDims; ii++) {
	df[ii] = sqrt(df[ii])*sqrtDT;
      }
    }
  } else {
    // matrix diffusion
    if(!sdeModel::sdDiff) {
      // var scale
      chol_decomp(df, df, sdeModel::nDims);
    }
    U_mult(df, sqrtDT, sdeModel::nDims);
  }
  return;
}

// Euler approximation mean and variance
inline void mvEuler(double *mean, double *sd,
		    double *x, double dT, double sqrtDT,
		    double *theta, sdeModel *sde) {
  // mean = x + drift(x,t,theta)*dT
  sde->sdeDr(mean, x, theta);
  v_mult(mean, dT, sdeModel::nDims);
  for(int ii = 0; ii < sdeModel::nDims; ii++) {
    mean[ii] += x[ii];
  }
  // sd = diff(x,t,theta)*sqrt(dT)
  sde->sdeDf(sd, x, theta);
  scaleDiff(sd, sqrtDT);
  return;
}

#endif
