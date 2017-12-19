#ifndef sdeLogLik_h
#define sdeLogLik_h 1

/* -------------------------------------------------------------------------

this is a class which defines the complete data and whatever it needs
to evaluate the log-density.
the log-density of an sde is essentially the sum of conditional normals,
each having its own mean and variance for every point.
i.e.,
logDens(x | theta) = sum dnorm(x_n+1 | mu(x_n, theta), sd(x_n, theta)).

so the constructor should create an array of sdeModels, each with enough
storage to compute log-densities.

public members:
* nComp, tSeq, dT, sqrtDT
* Drift(x, t, theta, i), Diff(x, t, theta, i)
* EulerMV(x, theta, i)
* loglik(x, theta)
*/

#include "LinAlgUtils.h"
#include "mvnUtils.h"
#include "sdeUtils.h"

// parallel implementation
#ifdef _OPENMP
#include <omp.h>
#else
inline int omp_get_thread_num(void) {return 0;}
#endif

template <class sMod>
class sdeLogLik {
 protected:
  int nDims2;
  // >= 1: no sense in "disabling" omp at runtime (instead of at compile time)
  int nCores;
  //void init(int n); // initialize temporary storage
 public:
  double *propMean, *propSd, *propZ; // for storing normal calculations
  sMod *sde; // for storing drift and diffusion calculations
  int nDims, nParams; // internal representations taken from sdeModel
  int nComp; // number of observations INCLUDING first one
  double *dT, *sqrtDT; // times
  // log-density
  double loglik(double *theta, double *x);
  //double loglikPar(double *theta, double *x);
  // constructor and destructor
  //sdeLogLik(int n, double *dt);
  sdeLogLik(int n, double *dt, int ncores);
  ~sdeLogLik();
};


// constructors ------------------------------------------------------------

/* inline sdeLogLik::sdeLogLik(int n, double *dt) { */
/*   nComp = n; */
/*   nDims = sdeModel::nDims; */
/*   nDims2 = sdeModel::diagDiff ? nDims : nDims*nDims; */
/*   nParams = sdeModel::nParams; */
/*   // create storage space */
/*   sde = new sdeModel[nComp]; */
/*   propMean = new double[nComp*nDims]; */
/*   propSd = new double[nComp*nDims2]; */
/*   propZ = new double[nComp*nDims]; */
/*   dT = new double[nComp]; */
/*   sqrtDT = new double[nComp]; */
/*   // timing */
/*   for(int ii=0; ii<nComp-1; ii++) { */
/*     dT[ii] = dt[ii]; */
/*     sqrtDT[ii] = sqrt(dT[ii]); */
/*   } */
/* } */

template <class sMod>
inline sdeLogLik<sMod>::sdeLogLik(int n, double *dt, int ncores) {
  nComp = n;
  nCores = ncores;
  nDims = sMod::nDims;
  nDims2 = sMod::diagDiff ? nDims : nDims*nDims;
  nParams = sMod::nParams;
  // create storage space
  sde = new sMod[nCores];
  propMean = new double[nCores*nDims];
  propSd = new double[nCores*nDims2];
  propZ = new double[nComp*nDims]; // RV draws can't be parallelized
  dT = new double[nComp];
  sqrtDT = new double[nComp];
  // timing
  for(int ii=0; ii<nComp-1; ii++) {
    dT[ii] = dt[ii];
    sqrtDT[ii] = sqrt(dT[ii]);
  }
}

// Destructor
template <class sMod>
inline sdeLogLik<sMod>::~sdeLogLik() {
  delete [] sde;
  delete [] propMean;
  delete [] propSd;
  delete [] propZ;
  delete [] dT;
  delete [] sqrtDT;
}

/* // full log-likelihood evaluation */
/* inline double sdeLogLik::loglik(double *theta, double *x) { */
/*   double ll = 0; */
/*   // *** PARALLELIZABLE FOR-LOOP *** */
/*   for(int ii = 0; ii < nComp-1; ii++) { */
/*     mvEuler(&propMean[ii*nDims], &propSd[ii*nDims2], */
/* 	    &x[ii*nDims], dT[ii], sqrtDT[ii], theta, &sde[ii]); */
/*     ll += lmvn(&x[(ii+1)*nDims], &propZ[ii*nDims], */
/* 	       &propMean[ii*nDims], &propSd[ii*nDims2], nDims); */
/*   } */
/*   return(ll); */
/* } */

// full log-likelihood evaluation
template <class sMod>
inline double sdeLogLik<sMod>::loglik(double *theta, double *x) {
  double ll = 0;
  // *** PARALLELIZABLE FOR-LOOP ***
  #ifdef _OPENMP
#pragma omp parallel for num_threads(nCores) reduction(+: ll) if(nCores > 1)
  #endif
  for(int ii = 0; ii < nComp-1; ii++) {
    int iCore = omp_get_thread_num();
    mvEuler<sMod>(&propMean[iCore*nDims], &propSd[iCore*nDims2],
	    &x[ii*nDims], dT[ii], sqrtDT[ii], theta, &sde[iCore]);
    ll += lmvn<sMod>(&x[(ii+1)*nDims], &propZ[ii*nDims],
		     &propMean[iCore*nDims], &propSd[iCore*nDims2]);
  }
  return(ll);
}

#endif
