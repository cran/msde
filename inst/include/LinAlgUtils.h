#ifndef LinAlgUtils_h
#define LinAlgUtils_h

/* // class for proposal mean and variances (or rather sd's) */
/* // also create a dummy space of size mean. */
/* class propMV { */
/* public: */
/*   int nDims; */
/*   double *mean, *sd, *z; */
/*   propMV(int); */
/*   ~propMV(); */
/* }; */

/* inline propMV::propMV(int d) { */
/*   nDims = d; */
/*   mean = new double[nDims]; */
/*   sd = new double[nDims*nDims]; */
/*   z = new double[nDims]; */
/*   // initialize */
/*   int ii; */
/*   for(ii = 0; ii < nDims; ii++) { */
/*     mean[ii] = 0.0; */
/*     z[ii] = 0.0; */
/*   } */
/*   for(ii = 0; ii < nDims*nDims; ii++) { */
/*     sd[ii] = 0.0; */
/*   } */
/* } */

/* inline propMV::~propMV() { */
/*   delete [] mean; */
/*   delete [] sd; */
/*   delete [] z; */
/* } */

// multiply vector by scalar a
inline void v_mult(double *v, double a, int n) {
  for(int ii=0; ii<n; ii++) {
    v[ii] *= a;
  }
  return;
}
// multiply an upper triangular matrix by scalar a
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

// cholesky decomposition of a symmetric, positive definite matrix.
// returns a vector of the *Upper Triangular* cholesy factor, leaving the other elements of the array unchanged.
// in other words, U' %*% U \neq A, but lowTri(U') %*% upTri(U) = A.
// both U and A are stacked by COLUMN
// can be performed IN-PLACE, i.e., with U and A refering to same memory location.
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
