#ifndef sdePrior_h
#define sdePrior_h 1

//[[Rcpp::depends("msde")]]

#include "sdeModel.h"

// this is specific to mvnPrior
#include <mvnUtils.h>

class sdePrior {
 private:
  static const int nDims = sdeModel::nDims;
  static const int nParams = sdeModel::nParams;
  int nRV, nParamRV, nDataRV; // number of active variables of each type
  int *paramId, *dataId; // index vectors
  double *mean, *cholSd;
  double *tmpX, *tmpZ;
 public:
  double logPrior(double *theta, double *x);
  sdePrior(double **phi, int nArgs, int *nEachArg);
  ~sdePrior();
};

// constructor
inline sdePrior::sdePrior(double **phi, int nArgs, int *nEachArg) {
  int ii;
  nRV = nEachArg[0];
  nParamRV = nEachArg[2];
  nDataRV = nEachArg[3];
  if(nRV > 0) {
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

// destructor
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
