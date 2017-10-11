#ifndef mvnUtils_h
#define mvnUtils_h 1

// main functions are: xmvn, zmvn, lmvn.  have one for cholSd's and sd's
// which are just diagonal.

// x = sd * z + mean.
// NOTE: sd = lowTri(cholSD'), i.e. only upper triangular entries of cholSD are used
inline void xmvn_chol(double *x, double *z,
		      double *mean, double *cholSd, int n) {
  int ii, jj, colI;
  for(ii=0; ii<n; ii++) {
    colI = n*ii;
    x[ii] = 0;
    for(jj = 0; jj <= ii; jj++) x[ii] += cholSd[colI + jj] * z[jj];
    x[ii] += mean[ii];
  }
  return;
}

// z = sd^{-1} * (x - mean).  only calculates first nMax values of z.
// NOTE: sd = lowTri(cholSD'), i.e. only upper triangular entries of cholSD are used
inline void zmvn_chol(double *z, double *x,
		      double *mean, double *cholSd, int n, int nMax) {
  int ii, jj, colI;
  double tmpSum;
  for(ii=0; ii<nMax; ii++) z[ii] = x[ii] - mean[ii];
  // forward substitution
  for(ii=0; ii<nMax; ii++) {
    colI = n*ii;
    tmpSum = 0.0;
    for(jj=0; jj<ii; jj++) tmpSum += cholSd[colI + jj] * z[jj];
    z[ii] = (z[ii] - tmpSum)/cholSd[colI + ii];
  }
  return;
}

// log-normal density evaluation.  z[] is required as temporary storage of residuals.
// i.e., z = sd^{-1} * (x - mean)
// NOTE: sd = lowTri(cholSD'), i.e. only upper triangular entries of cholSD are used
// TODO: include pi factor
inline double lmvn_chol(double *x, double *z,
			double *mean, double *cholSd, int n) {
  double ssq = 0.0; // sum(z^2)
  double ldC = 0.0; // log(det(cholSd))
  double resi, tmpSum, val;
  int ii, colI, jj;
  // forward substitution
  colI = 0;
  for(ii=0; ii<n; ii++) {
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

// --- diagonal versions -------------------------------------------------------

inline void xmvn_diag(double *x, double *z,
		      double *mean, double *diagSd, int n) {
  for(int ii=0; ii<n; ii++) {
    x[ii] = diagSd[ii]*z[ii] + mean[ii];
  }
}

inline void zmvn_diag(double *z, double *x,
		      double *mean, double *diagSd, int n, int nMax) {
  for(int ii=0; ii<nMax; ii++) {
    z[ii] = (x[ii] - mean[ii])/diagSd[ii];
  }
  return;
}

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
