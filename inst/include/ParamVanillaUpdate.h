#ifndef ParamVanillaUpdate_h
#define ParamVanillaUpdate_h 1

#include "rngUtils.h"
#include "sdeMCMC.h"

// componentwise vanilla MH parameter updates
template <class sMod, class sPi>
  inline void sdeMCMC<sMod, sPi>::paramVanillaUpdate(double *jumpSd,
						     int *paramAccept) {
  double acc, currLoglik, propLoglik;
  int ii;
  // initialize
  for(ii = 0; ii < nParams; ii++) {
    propTheta[ii] = currTheta[ii];
  }
  currLoglik = loglik(currTheta, currX);
  // random walk metropolis updates
  for(ii = 0; ii < nParams; ii++) {
    if(!fixedTheta[ii]) {
      // proposal
      propTheta[ii] = currTheta[ii] + jumpSd[ii] * sdeRNG::rnorm();
      // only calculate acceptance if valid
      if(sde[0].isValidParams(propTheta)) {
	// likelihood
	propLoglik = loglik(propTheta, currX);
	acc = propLoglik - currLoglik;
	// prior
	acc += prior->logPrior(propTheta, currX);
	acc -= prior->logPrior(currTheta, currX);
	// acceptance rate
	if(exp(acc) >= sdeRNG::runif()) {
	  currTheta[ii] = propTheta[ii];
	  currLoglik = propLoglik;
	  paramAccept[ii]++;
	}
      }
      // propTheta and currTheta should only differ by one element
      propTheta[ii] = currTheta[ii];
    }
  }
  return;
}


#endif
