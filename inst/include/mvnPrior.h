/// @file mvnPrior.h

#ifndef mvnPrior_h
#define mvnPrior_h 1

#include <mvnUtils.h>
//#include "sdeModel.h"


/// Multivariate normal prior distribution.
///
/// The arguments to the prior consist of:
///
/// 1. The mean of the prior random variables, i.e,. the non-fixed parameters and initial missing data points (params and data concatenated together).
/// 2. The variance (or rather `cholSd`) of these random variables. 
/// 3. The index of each non-fixed parameter.
/// 4. The index of each non-fixed initial missing data point.
class sdePrior {
 private:
  //static const int nDims = sdeModel::nDims;
  //static const int nParams = sdeModel::nParams;
  int nRV, nParamRV, nDataRV; // number of active variables of each type
  int *paramId, *dataId; // index vectors
  double *mean, *cholSd;
  double *tmpX, *tmpZ;
 public:
  /// Evaluate the prior log-density.
  double logPrior(double *theta, double *x);
  /// Constructor.
  sdePrior(double **phi, int nArgs, int *nEachArg);
  ~sdePrior();
};


/// @param[in] phi Pointer to storage for each argument of the prior.
/// @param[in] nArgs Number of prior arguments, which is 4.
/// @param[in] nEachArg Pointer to length of each argument of the prior.
inline sdePrior::sdePrior(double **phi, int nArgs, int *nEachArg) {
  int ii;
  nRV = nEachArg[0];
  if(nRV > 0) {
    nParamRV = nEachArg[2];
    nDataRV = nEachArg[3];
    mean = new double[nRV];
    cholSd = new double[nRV*nRV];
    tmpX = new double[nRV];
    tmpZ = new double[nRV];
    // allocate mean and variance vectors
    for(ii=0; ii<nRV; ii++) {
      mean[ii] = phi[0][ii];
    }
    for(ii=0; ii<nRV*nRV; ii++) {
      cholSd[ii] = phi[1][ii];
    }
    // index arrays for each type of variable
    if(nParamRV > 0) {
      paramId = new int[nParamRV];
      for(ii=0; ii<nParamRV; ii++) {
	paramId[ii] = (int) phi[2][ii] - 1;
      }
    }
    if(nDataRV > 0) {
      dataId = new int[nDataRV];
      for(ii=0; ii<nDataRV; ii++) {
	dataId[ii] = (int) phi[3][ii] - 1;
      }
    }
  }
}

inline sdePrior::~sdePrior() {
  if(nRV > 0) {
    delete [] mean;
    delete [] cholSd;
    delete [] tmpX;
    delete [] tmpZ;
    if(nParamRV > 0) {
      delete [] paramId;
    }
    if(nDataRV > 0) {
      delete [] dataId;
    }
  }
}

/// @param[in] theta Array of parameter values (all of them, fixed ones will be ignored).
/// @param[in] x Array of data values (all of them, fixed ones will be ignored).
///
/// @return Value of the prior log-density.
inline double sdePrior::logPrior(double *theta, double *x) {
  if(nRV == 0) return(0.0);
  double lp;
  int ii;
  for(ii=0; ii<nParamRV; ii++) {
    tmpX[ii] = theta[paramId[ii]];
  }
  for(ii=0; ii<nDataRV; ii++) {
    tmpX[nParamRV+ii] = x[dataId[ii]];
  }
  lp = lmvn_chol(tmpX, tmpZ, mean, cholSd, nRV);
  return(lp);
}


#endif
