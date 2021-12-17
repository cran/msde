/// @file LinAlgUtils.h
///
/// Various linear algebra utilities.

#ifndef LinAlgUtils_h
#define LinAlgUtils_h

/// Multiply vector by a scalar in place.
///
/// @param[in,out] v Vector of length `n`.
/// @param[in] a Scalar by which to multiply `v`.
/// @param[in] n Length of `v`.
inline void v_mult(double *v, double a, int n) {
  for(int ii=0; ii<n; ii++) {
    v[ii] *= a;
  }
  return;
}

/// Multiply upper triangular matrix by a scalar in place.
///
/// @param[in,out] U upper triangular matrix as an array of size `n^2`.
/// @param[in] Scalar by which to multiply `U`.
/// @param[in] n Size of `U`.
inline void U_mult(double *U, double a, int n) {
  int ii,jj,colI;
  for(ii=0; ii<n; ii++) {
    colI = ii*n;
    for(jj=0; jj<=ii; jj++) {
      U[colI+jj] *= a;
    }
  }
  return;
}

/// Cholesky decomposition of a symmetric, positive definite matrix.
///
/// Computes the upper triangular cholesky factor, leaving the other elements of the input unchanged.
/// In other words, `U' %*% U \= A`, but `lowTri(U') %*% upTri(U) = A`.
/// Both `U` and `A` are stacked by column.
/// Can be performed IN-PLACE, i.e., with U and A refering to same memory location.
///
/// @param[out] U Array of size `n^2` into which to compute the Cholesky decomposition.
/// @param[in] A Array of size `n^2` containing the positive definite matrix.  Only uses the upper triangular elements of `A`.  Changing any of the other elements won't affect the result.
/// @param[in] n Size of `A`.
inline void chol_decomp(double *U, double *A, int n) {
  int ii, jj, kk, colI, colJ;
  double tmpSum, tmpInv;
  for(ii = 0; ii < n; ii++) {
    colI = ii*n;
    tmpSum = 0.0;
    for(kk = 0; kk < ii; kk++) {
      tmpSum += U[colI + kk] * U[colI + kk];
    }
    tmpInv = sqrt(A[colI + ii] - tmpSum);
    U[colI + ii] = tmpInv;
    tmpInv = 1.0/tmpInv;
    for(jj = ii+1; jj < n; jj++) {
      colJ = jj*n;
      tmpSum = 0.0;
      for(kk = 0; kk < ii; kk++) tmpSum += U[colJ + kk] * U[colI + kk];
      U[colJ + ii] = tmpInv * (A[colJ + ii] - tmpSum);
    }
  }
  return;
}

#endif
