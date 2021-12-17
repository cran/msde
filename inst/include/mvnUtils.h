/// @file mvnUtils.h
///
/// Utilities for the multivariate normal distribution.
///
/// Suppose we have `x ~ N(mu, Sigma)` and `sd` is the lower Cholesky factor of `Sigma = sd * sd'`.  The utilities provided are:
///
/// - `xmvn`: Calculate `x = sd * z + mu`.
/// - `zmvn`: Calculate `z = sd^{-1} * (x-mu)`.
/// - `lmvn`: Calculate the log-pdf of `x ~ N(mu, Sigma)`.
///
/// Each utility has an implementation for dense `Sigma` and for diagonal `Sigma`.
///
/// In all functions below, `sd = cholSd'`, i.e., `cholSd` is the uppor triangular Cholesky factor of `Sigma = cholSd' * cholSd`.  Also, only the upper triangular entries of `cholSd` are used, i.e., the lower triangular entries can be arbitrary without affecting the output.

#ifndef mvnUtils_h
#define mvnUtils_h 1

/*
// .5 * log(2*pi)
#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI	0.918938533204672741780329736406	
#endif
*/

/// Calculate `x = sd * z + mean` for lower triangular `sd`.
///
/// Can specify that the calculation is only performed to obtain `x[iFirst]` through `x[iLast]`.
/// 
/// @param[out] x Vector of length `n`.
/// @param[in] z Vector of length `n`.
/// @param[in] mean Vector of length `n`.
/// @param[in] cholSd Vector of length `n^2`, which gets rearranged to a `n x n` upper triangular matrix.
/// @param[in] n Size of the problem.
/// @param[in] iFirst The first desired element of the output.
/// @param[in] iLast The last desired element of the output.
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

/// Calculate `x = sd * z + mean` for lower triangular `sd`.
///
/// Simplified version: assumes `iFirst = 0` and `iLast = n`.
inline void xmvn_chol(double *x, double *z,
		      double *mean, double *cholSd, int n) {
  xmvn_chol(x, z, mean, cholSd, n, 0, n);
  return;
}

/// Calculate `z = sd^{-1} * (x - mean)` for lower triangular `sd`.
///
/// Can specify that the calculation is only performed to obtain `z[0], ..., z[iLast-1]`.
///
/// @param[out] z Vector of length `n`.
/// @param[in] x Vector of length `n`.
/// @param[in] mean Vector of length `n`.
/// @param[in] cholSd Vector of length `n^2`, which gets rearranged to a `n x n` upper triangular matrix.
/// @param[in] n Size of the problem.
/// @param[in] iLast The last desired element of the output.
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

/// Calculate `z = sd^{-1} * (x - mean)` for lower triangular `sd`.
///
/// Simplified version: assumes `iLast = n`.
inline void zmvn_chol(double *z, double *x,
		      double *mean, double *cholSd, int n) {
  zmvn_chol(z, x, mean, cholSd, n, n);
  return;
}


/// Calculate the log-density of `x ~ N(mean, sd * sd')` for lower-triangular `sd`.
///
/// Can specify that the log-density is only evaluated for `x[0], ..., x[iLast-1]`.
/// @notes
/// - `z` is required for temporary storage of residuals, i.e., `z = sd^{-1} * (x - mean)`.
/// - **TODO:** Include factor of `pi`.
///
/// @param[in] x Vector of length `n`.
/// @param[in] z Vector of length `n`.
/// @param[in] mean Vector of length `n`.
/// @param[in] cholSd Vector of length `n^2`, which gets rearranged to a `n x n` upper triangular matrix.
/// @param[in] n Size of the problem.
/// @param[in] iLast Specifies which marginal distribution of `x` to return.
///
/// @return The log-density of `x[0:iLast] ~ N(mean[0:iLast], (sd * sd')[...])`.
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

/// Calculate the log-density of `x ~ N(mean, sd * sd')` for lower-triangular `sd`.
///
/// Simplified version: assumes that `iLast = n`.
inline double lmvn_chol(double *x, double *z,
			double *mean, double *cholSd, int n) {
  return lmvn_chol(x, z, mean, cholSd, n, n);
}


// --- diagonal versions -------------------------------------------------------

// don't need iLast versions for zmvn and lmvn

/// Calculate `x = sd * z + mean` for diagonal `sd`.
///
/// Can specify that the calculation is only performed to obtain `x[iFirst]` through `x[iLast]`.
/// 
/// @param[out] x Vector of length `n`.
/// @param[in] z Vector of length `n`.
/// @param[in] mean Vector of length `n`.
/// @param[in] diagSd Vector of length `n` such that `sd = diag(diagSd)`.
/// @param[in] n Size of the problem.
/// @param[in] iFirst The first desired element of the output.
/// @param[in] iLast The last desired element of the output.
inline void xmvn_diag(double *x, double *z,
		      double *mean, double *diagSd,
		      int n, int iFirst, int iLast) {
  for(int ii=iFirst; ii<iLast; ii++) {
    x[ii] = diagSd[ii]*z[ii] + mean[ii];
  }
  return;
}

/// Calculate `x = sd * z + mean` for diagonal `sd`.
///
/// Simplified version: assumes `iFirst = 0` and `iLast = n`.
inline void xmvn_diag(double *x, double *z,
		      double *mean, double *diagSd,
		      int n) {
  xmvn_diag(x, z, mean, diagSd, n, 0, n);
  return;
}


/// Calculate `z = sd^{-1} * (x - mean)` for diagonal `sd`.
///
/// @param[out] z Vector of length `n`.
/// @param[in] x Vector of length `n`.
/// @param[in] mean Vector of length `n`.
/// @param[in] diagSd Vector of length `n` such that `sd = diag(diagSd)`.
/// @param[in] n Size of the problem.
inline void zmvn_diag(double *z, double *x,
		      double *mean, double *diagSd, int n) {
  for(int ii=0; ii<n; ii++) {
    z[ii] = (x[ii] - mean[ii])/diagSd[ii];
  }
  return;
}


/// Calculate the log-density of `x ~ N(mean, sd * sd')` for diagonal `sd`.
///
/// @param[in] x Vector of length `n`.
/// @param[in] mean Vector of length `n`.
/// @param[in] diagSd Vector of length `n` such that `sd = diag(diagSd)`.
/// @param[in] n Size of the problem.
///
/// @return The log-density of `x ~ N(mean, sd * sd')`.
inline double lmvn_diag(double *x,
			double *mean, double *diagSd, int n) {
  double z;
  double ssq = 0.0;
  double ldC = 0.0;
  for(int ii=0; ii<n; ii++) {
    z = (x[ii] - mean[ii])/diagSd[ii];
    ssq += z*z;
    ldC += log(diagSd[ii]);
  }
  return(-(.5*ssq + ldC));
}

#endif
