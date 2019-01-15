#ifndef mvnUtils_h
#define mvnUtils_h 1

/*
// .5 * log(2*pi)
#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI	0.918938533204672741780329736406	
#endif
*/

// main functions are: xmvn, zmvn, lmvn.  have one for cholSd's and sd's
// which are just diagonal.

// in all functions below, sd = lowTri(cholSD'), i.e., cholesky factor, and
// only upper triangular entries of cholSD are used.

// x = sd * z + mean
// only apply to indices between iFirst and iLast-1, inclusively.
inline void xmvn_chol(double *x, double *z,
		      double *mean, double *cholSd,
		      int n, int iFirst, int iLast) {
  int ii, jj, colI;
  for(ii=iFirst; ii<iLast; ii++) {
    colI = n*ii;
    x[ii] = 0;
    for(jj = 0; jj <= ii; jj++) x[ii] += cholSd[colI + jj] * z[jj];
    x[ii] += mean[ii];
  }
  return;
}

// simplified version
inline void xmvn_chol(double *x, double *z,
		      double *mean, double *cholSd, int n) {
  xmvn_chol(x, z, mean, cholSd, n, 0, n);
  return;
}

// z = sd^{-1} * (x - mean).
// only calculates first iLast values of z.
inline void zmvn_chol(double *z, double *x,
		      double *mean, double *cholSd,
		      int n, int iLast) {
  int ii, jj, colI;
  double tmpSum;
  for(ii=0; ii<iLast; ii++) z[ii] = x[ii] - mean[ii];
  // forward substitution
  for(ii=0; ii<iLast; ii++) {
    colI = n*ii;
    tmpSum = 0.0;
    for(jj=0; jj<ii; jj++) tmpSum += cholSd[colI + jj] * z[jj];
    z[ii] = (z[ii] - tmpSum)/cholSd[colI + ii];
  }
  return;
}

// simplified version
inline void zmvn_chol(double *z, double *x,
		      double *mean, double *cholSd, int n) {
  zmvn_chol(z, x, mean, cholSd, n, n);
  return;
}


// log-normal density evaluation.
// z[] is required as temporary storage of residuals,
// i.e., z = sd^{-1} * (x - mean)
// only evaluate marginal PDF of first iLast observations.
// TODO: include pi factor
inline double lmvn_chol(double *x, double *z,
			double *mean, double *cholSd, int n, int iLast) {
  double ssq = 0.0; // sum(z^2)
  double ldC = 0.0; // log(det(cholSd))
  double resi, tmpSum, val;
  int ii, colI, jj;
  // forward substitution
  colI = 0;
  for(ii=0; ii<iLast; ii++) {
    resi = x[ii] - mean[ii];
    tmpSum = 0.0;
    for(jj = 0; jj < ii; jj++) tmpSum += cholSd[colI + jj] * z[jj];
    val = (resi - tmpSum) / cholSd[colI + ii];
    ldC += log(cholSd[colI + ii]);
    z[ii] = val;
    ssq += (val * val);
    colI += n;
  }
  return(-(.5*ssq + ldC));
}

// simplified version
inline double lmvn_chol(double *x, double *z,
			double *mean, double *cholSd, int n) {
  return lmvn_chol(x, z, mean, cholSd, n, n);
}


// --- diagonal versions -------------------------------------------------------

// don't need iLast versions for zmvn and lmvn

inline void xmvn_diag(double *x, double *z,
		      double *mean, double *diagSd,
		      int n, int iFirst, int iLast) {
  for(int ii=iFirst; ii<iLast; ii++) {
    x[ii] = diagSd[ii]*z[ii] + mean[ii];
  }
  return;
}
inline void xmvn_diag(double *x, double *z,
		      double *mean, double *diagSd,
		      int n) {
  xmvn_diag(x, z, mean, diagSd, n, 0, n);
  return;
}

inline void zmvn_diag(double *z, double *x,
		      double *mean, double *diagSd, int n) {
  for(int ii=0; ii<n; ii++) {
    z[ii] = (x[ii] - mean[ii])/diagSd[ii];
  }
  return;
}
/* inline void zmvn_diag(double *z, double *x, */
/* 		      double *mean, double *diagSd, int n) { */
/*   zmvn_diag(z, x, mean, diagSd, n, n); */
/*   return; */
/* } */

inline double lmvn_diag(double *x, double *z,
			double *mean, double *diagSd, int n) {
  double ssq = 0.0;
  double ldC = 0.0;
  for(int ii=0; ii<n; ii++) {
    z[ii] = (x[ii] - mean[ii])/diagSd[ii];
    ssq += z[ii]*z[ii];
    ldC += log(diagSd[ii]);
  }
  return(-(.5*ssq + ldC));
}

#endif
