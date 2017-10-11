#ifndef sdeMCMC_h
#define sdeMCMC_h 1


// posterior inference for msde's

// for speed considerations, and since it is unlikely that users will
// be building their own MCMC routines, all transitions will be members
// of the same object.
// the object will inherit from sdeLogLik.

#include "mvnUtils.h"
#include "sdeLogLik.h"
//#include "mcmcUtils.h"
//#include "sdePrior.h"

template <class sMod, class sPi>
class sdeMCMC : public sdeLogLik<sMod> {
  //int nDims2, nCores; // inherited from sdeLogLik
  using sdeLogLik<sMod>::nDims2;
  using sdeLogLik<sMod>::nCores;
  int *missInd;
  int nMiss, nMiss0, nMissN;
  sPi *prior;
  void mvEraker(double *mean, double *sd,
		double *x0, double *x2,
		double b, double b2, double *theta,
		sMod *sde);
 public:
  //int nComp, nDims, nParams; // inherited from sdeLogLik
  using sdeLogLik<sMod>::nComp;
  using sdeLogLik<sMod>::nDims;
  using sdeLogLik<sMod>::nParams;
  double *currFull, *propFull, *propAccept;
  double *currX, *propX, *currTheta, *propTheta;
  double *propU; // for acceptance rates (predrawn in parallel version)
  //double *dT, *sqrtDT; // inherited from sdeLogLik
  using sdeLogLik<sMod>::dT;
  using sdeLogLik<sMod>::sqrtDT;
  double *B, *sqrtB;
  int *nObsComp;
  bool *fixedTheta;
  //sMod *sde; // inherited sdeLogLik
  using sdeLogLik<sMod>::sde;
  // double *propMean, *propSd, *propZ; // inherited from sdeLogLik
  using sdeLogLik<sMod>::propMean;
  using sdeLogLik<sMod>::propSd;
  using sdeLogLik<sMod>::propZ;
  // double loglik(double *theta, double *x); // inherited from sdeLogLik
  using sdeLogLik<sMod>::loglik;
  void missGibbsUpdate(double *jumpSd, int *gibbsAccept, int *paramAccept);
  void paramVanillaUpdate(double *jumpSd, int *paramAccept);
  sdeMCMC(int n, double *dt, double *xInit, double *thetaInit,
	  int *xIndex, bool *thetaIndex,
	  double **phi, int nArgs, int *nEachArg, int ncores);
  ~sdeMCMC();
};

template <class sMod, class sPi>
  inline sdeMCMC<sMod, sPi>::sdeMCMC(int n, double *dt,
				     double *xInit, double *thetaInit,
				     int *xIndex, bool *thetaIndex,
				     double **phi,
				     int nArgs, int *nEachArg,
				     int ncores) :
  sdeLogLik<sMod>(n, dt, ncores) {
  int ii, jj;
  int nDims2, nCores; // inherited from sdeLogLik
  int nComp, nDims, nParams; // inherited from sdeLogLik
  double *dT, *sqrtDT; // inherited from sdeLogLik
  // inherited template parameters
  nDims2 = this->nDims2;
  nCores = this->nCores;
  nComp = this->nComp;
  nDims = this->nDims;
  nParams = this->nParams;
  dT = this->dT;
  sqrtDT = this->sqrtDT;
  // memory allocation
  B = new double[nComp];
  sqrtB = new double[nComp];
  for(ii=1; ii<nComp-1; ii++) {
    //sqrtDT[ii] = sqrt(dT[ii]);
    //if(ii > 0) {
    B[ii] = dT[ii]/(dT[ii-1] + dT[ii]);
    sqrtB[ii] = sqrt((1-B[ii]) * dT[ii]);
    //}
  }
  // data
  currFull = new double[nComp*nDims + nParams];
  propFull = new double[nCores*nDims + nParams];
  propAccept = new double[nCores];
  propU = new double[nComp];
  currX = currFull + nParams;
  propX = propFull + nParams;
  nObsComp = new int[nComp];
  // initialize
  for(ii=0; ii<nComp; ii++) {
    nObsComp[ii] = xIndex[ii];
    propU[ii] = 0.0;
    for(jj=0; jj<nDims; jj++) {
      currX[ii*nDims + jj] = xInit[ii*nDims + jj];
    }
  }
  for(ii=0; ii<nCores; ii++) {
    propAccept[ii] = 0.0;
    for(jj=0; jj<nDims; jj++) {
      propX[ii*nDims + jj] = currX[ii*nDims + jj];
    }
  }
  // missing data
  nMiss0 = nDims-nObsComp[0]; // unobserved states in first observation
  nMissN = nDims-nObsComp[nComp-1]; // unobserved states in last observation
  // identify missing _intermediate_ data indices,
  // i.e. at least one component to update and not first or last observation
  nMiss = 0;
  for(ii = 1; ii < nComp-1; ii++) {
    nMiss += (nObsComp[ii] < nDims);
  }
  missInd = new int[nMiss + (nMiss == 0)];
  jj = 0;
  for(ii = 1; ii < nComp-1; ii++) {
    if(nObsComp[ii] < nDims) {
      missInd[jj++] = ii;
    }
  }
  // parameters
  fixedTheta = new bool[nParams];
  currTheta = currFull;
  propTheta = propFull;
  for(ii=0; ii<nParams; ii++) {
    //fixedTheta[ii] = (thetaIndex[ii] != 0); // R logicals stored as C++ ints
    fixedTheta[ii] = thetaIndex[ii]; // R logicals stored as C++ ints
    currTheta[ii] = thetaInit[ii];
    propTheta[ii] = currTheta[ii];
  }
  // prior
  // any advantage to passing pointer?
  prior = new sPi(phi, nArgs, nEachArg);
}

template <class sMod, class sPi>
inline sdeMCMC<sMod, sPi>::~sdeMCMC() {
  delete [] B;
  delete [] sqrtB;
  delete [] currFull;
  delete [] propFull;
  delete [] propAccept;
  delete [] propU;
  delete [] missInd;
  delete [] nObsComp;
  delete [] fixedTheta;
  delete prior;
}

#include <MissGibbsUpdate.h>
#include <ParamVanillaUpdate.h>

#endif
