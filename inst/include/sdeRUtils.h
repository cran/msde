/// @file sdeRUtils.h

#ifndef sdeRUtils_h
#define sdeRUtils_h

#include <Rcpp.h>
typedef Rcpp::LogicalVector Logical;
typedef Rcpp::NumericVector Numeric;
typedef Rcpp::IntegerVector Integer;
typedef Rcpp::List List;
//using namespace Rcpp;
#include "sdeUtils.h"
#include "sdeLogLik.h"
#include "sdeRobj.h"

/// @return Number of SDE dimensions.
template <class sMod, class sPi>
  inline int sdeRobj<sMod, sPi>::get_nDims() {
  return sMod::nDims;
}

/// @return Number of SDE parameters.
template <class sMod, class sPi>
  inline int sdeRobj<sMod, sPi>::get_nParams() {
  return sMod::nParams;
}

/// @param[in] xIn Data vector of length `nDims` or `nDims * nReps`.
/// @param[in] thetaIn Parameter vector of length `nParams` or `nParams * nReps`.
/// @param[in] singleX Whether `xIn` is of length `nDims`.
/// @param[in] singleTheta Whether `thetaIn` is of length `nParams`.
/// @param[in] nReps Number of data/parameter pairs to vectorize over.
///
/// @return Logical vector of length `nReps`.
template <class sMod, class sPi>
  inline Logical sdeRobj<sMod, sPi>::isData(Numeric xIn, Numeric thetaIn,
					    bool singleX, bool singleTheta,
					    int nReps) {
  int nDims = sMod::nDims;
  int nParams = sMod::nParams;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  Logical validOut(nReps);
  sMod sde;
  for(int ii = 0; ii < nReps; ii++) {
    validOut[ii] = sde.isValidData(&x[ii*(!singleX)*nDims],
				   &theta[ii*(!singleTheta)*nParams]);
  }
  return validOut;
}

/// @param[in] thetaIn Parameter vector of length `nParams * nReps`.
/// @param[in] nReps Number of parameter vectors to vectorize over.
///
/// @return Logical vector of length `nReps`.
template <class sMod, class sPi>
  inline Logical sdeRobj<sMod, sPi>::isParams(Numeric thetaIn, int nReps) {
  int nParams = sMod::nParams;
  double *theta = REAL(thetaIn);
  Logical validOut(nReps);
  sMod sde;
  for(int ii = 0; ii < nReps; ii++) {
    validOut[ii] = sde.isValidParams(&theta[ii*nParams]);
  }
  return validOut;
}

/// @param[in] xIn Data vector of length `nDims` or `nDims * nReps`.
/// @param[in] thetaIn Parameter vector of length `nParams` or `nParams * nReps`.
/// @param[in] singleX Whether `xIn` is of length `nDims`.
/// @param[in] singleTheta Whether `thetaIn` is of length `nParams`.
/// @param[in] nReps Number of data/parameter pairs to vectorize over.
///
/// @return Vector of length `nDims * nReps`, of which the drift for a given data/parameter pair are contiguous elements.
///
/// @note Best to parallelize these from within R: easier and probably faster since memory can be allocated in parallel there.
template <class sMod, class sPi>
  inline Numeric sdeRobj<sMod, sPi>::Drift(Numeric xIn, Numeric thetaIn,
					   bool singleX, bool singleTheta,
					   int nReps) {
  int nDims = sMod::nDims;
  int nParams = sMod::nParams;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  Numeric drOut(nReps*nDims);
  double *dr = REAL(drOut);
  sMod sde;
  for(int ii = 0; ii < nReps; ii++) {
    sde.sdeDr(&dr[ii*nDims],
	      &x[ii*(!singleX)*nDims], &theta[ii*(!singleTheta)*nParams]);
  }
  return drOut;
}

/// @param[in] xIn Data vector of length `nDims` or `nDims * nReps`.
/// @param[in] thetaIn Parameter vector of length `nParams` or `nParams * nReps`.
/// @param[in] singleX Whether `xIn` is of length `nDims`.
/// @param[in] singleTheta Whether `thetaIn` is of length `nParams`.
/// @param[in] nReps Number of data/parameter pairs to vectorize over.
///
/// @return Vector of length `nDims^2 * nReps`, of which the diffusion matrix on the (upper) cholesky scale for a given data/parameter pair are contiguous elements.
template <class sMod, class sPi>
  inline Numeric sdeRobj<sMod, sPi>::Diff(Numeric xIn, Numeric thetaIn,
					  bool singleX, bool singleTheta,
					  int nReps) {
  int nDims = sMod::nDims;
  int nDims2 = nDims*nDims;
  int nParams = sMod::nParams;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  Numeric dfOut(nReps*nDims2);
  double *df = REAL(dfOut);
  sMod sde;
  int ii,jj;
  for(ii=0; ii<nReps; ii++) {
    sde.sdeDf(&df[ii*nDims2],
	      &x[ii*(!singleX)*nDims], &theta[ii*(!singleTheta)*nParams]);
    if(sMod::diagDiff) {
      if(sMod::sdDiff) {
	for(jj=1; jj<nDims; jj++) {
	  df[ii*nDims2 + jj*nDims + jj] = df[ii*nDims2+jj];
	}
      } else {
	for(jj=1; jj<nDims; jj++) {
	  df[ii*nDims2 + jj*nDims + jj] = sqrt(df[ii*nDims2+jj]);
	}
      }
    } else {
      if(!sMod::sdDiff) {
	chol_decomp(&df[ii*nDims2], &df[ii*nDims2], nDims);
      }
    }
  }
  return dfOut;
}

// SDE log-likelihood evaluation.
template <class sMod, class sPi>
  inline Numeric sdeRobj<sMod, sPi>::LogLik(Numeric xIn, Numeric dTIn,
					    Numeric thetaIn, int nComp,
					    int nReps, bool singleX,
					    bool singleTheta, int nCores) {
  int nDims = sMod::nDims;
  int nParams = sMod::nParams;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  Numeric llOut(nReps);
  double *ll = REAL(llOut);
  sdeLogLik<sMod> sdeLL(nComp, REAL(dTIn), nCores);
  for(int ii=0; ii<nReps; ii++) {
    ll[ii] = sdeLL.loglik(&theta[ii*(!singleTheta)*nParams],
			  &x[ii*(!singleX)*nDims*nComp]);
  }
  return llOut;
}

template <class sMod, class sPi>
  inline Numeric sdeRobj<sMod, sPi>::Prior(Numeric thetaIn, Numeric xIn,
					   bool singleTheta, bool singleX,
					   int nReps, List phiIn) {
  int nDims = sMod::nDims;
  int nParams = sMod::nParams;
  int ii;
  double *x = REAL(xIn);
  double *theta = REAL(thetaIn);
  PriorArgs priorArgs(phiIn);
  /* int nArgs = phiIn.length(); */
  /* double **phi = new double*[nArgs]; */
  /* int *nEachArg = new int[nArgs]; */
  /* for(ii=0; ii<nArgs; ii++) { */
  /*   if(Rf_isNull(phiIn[ii])) { */
  /*     nEachArg[ii] = 0; */
  /*   } else { */
  /*     nEachArg[ii] = as<NumericVector>(phiIn[ii]).length(); */
  /*     phi[ii] = REAL(phiIn[ii]); */
  /*   } */
  /* } */
  sPi prior(priorArgs.phi, priorArgs.nArgs, priorArgs.nEachArg);
  Numeric lpOut(nReps);
  double *lp = REAL(lpOut);
  // NOTE: this can't be parallelized because private storage is common
  // to parallelize need array of sdePrior objects
  for(ii=0; ii<nReps; ii++) {
    lp[ii] = prior.logPrior(&theta[ii*(!singleTheta)*nParams],
			    &x[ii*(!singleX)*nDims]);
  }
  /* delete [] phi; */
  /* delete [] nEachArg; */
  return lpOut;
}

#endif
