/// @file sdeLogLik.h

#ifndef sdeLogLik_h
#define sdeLogLik_h 1

#include "LinAlgUtils.h"
#include "mvnUtils.h"
#include "sdeUtils.h"
#include "sdeData.h"

// parallel implementation
#ifdef _OPENMP
#include <omp.h>
#else
inline int omp_get_thread_num(void) {return 0;}
#endif

/// Evaluates the loglikelihood for an SDE model.
template <class sMod>
class sdeLogLik : public sdeData<sMod> {
 protected:
  using sdeData<sMod>::nDims2;
  // >= 1: no sense in "disabling" omp at runtime (instead of at compile time)
  int nCores;
  //void init(int n); // initialize temporary storage
 public:
  // inherited members from sdeData
  using sdeData<sMod>::nObsComp;
  using sdeData<sMod>::propMean;
  using sdeData<sMod>::propSd;
  using sdeData<sMod>::propZ;
  using sdeData<sMod>::sde;
  using sdeData<sMod>::nDims;
  using sdeData<sMod>::nParams;
  using sdeData<sMod>::nComp;
  using sdeData<sMod>::dT;
  using sdeData<sMod>::sqrtDT;
  /// Calculate the log-density of the SDE model.
  double loglik(double *theta, double *x);
  /// Constructor.
  sdeLogLik(int ncomp, double *dt, int ncores);
  /// Constructor with missing data information.
  sdeLogLik(int ncomp, double *dt, int *nobscomp, int ncores);
  // ~sdeLogLik();
};


// constructors ------------------------------------------------------------


/// @param[in] ncomp Number of complete data observations.
/// @param[in] dt Array of interobservation times.
/// @param[in] ncores Number of cores for parallel processing.
template <class sMod>
inline sdeLogLik<sMod>::sdeLogLik(int ncomp, double *dt, int ncores) :
sdeData<sMod>(ncomp, dt, ncores, ncomp, ncores) {nCores = ncores;}


/// @param[in] ncomp Number of complete data observations.
/// @param[in] dt array of interobservation times
/// @param[in] ncores Number of cores for parallel processing.
template <class sMod>
inline sdeLogLik<sMod>::sdeLogLik(int ncomp, double *dt,
				  int *nobscomp, int ncores) :
sdeData<sMod>(ncomp, dt, nobscomp, ncores, ncomp, ncores) {nCores = ncores;}


/// @param[in] theta Array of parameter values.
/// @param[in] x Array of `nDims * nComp` SDE observations.
///
/// @return Value of the Euler-Maruyama log-density, up to a factor of `pi`.
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
