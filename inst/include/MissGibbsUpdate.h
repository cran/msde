#ifndef MissGibbsUpdate_h
#define MissGibbsUpdate_h 1

//[[Rcpp::depends("msde")]]
#include <mvnUtils.h>
#include <sdeMCMC.h>

// eraker proposal mean and standard deviatiation
// NOTE: sde = upper triangular cholesky factor
inline void sdeMCMC::mvEraker(double *mean, double *sd,
	      double *x0, double *x2,
	      double b, double b2, double *theta,
	      sdeModel *sde) {
  double b1 = 1.0-b;
  for(int ii=0; ii<sdeModel::nDims; ii++) {
    mean[ii] = x0[ii] * b + x2[ii] * b1;
  }
  sde->sdeDf(sd, x0, theta);
  scaleDiff(sd, b2);
  /* if(!sdeModel::sdDiff) { */
  /*   chol_decomp(sd, sd, sdeModel::nDims); */
  /* } */
  /* U_mult(sd, b2, sdeModel::nDims); */
  return;
}

inline void sdeMCMC::missGibbsUpdate(double *jumpSd, int *gibbsAccept,
				     int *paramAccept) {
  int ii, II, jj, JJ;
  int iCore;
  //int startII, endII;
  double *mean, *sd, *Z;
  // only elements in missInd are updated.
  // first and last components are handled separately.
  // Markov chain elements are conditionally independent,
  // so every other can be updated in parallel
  // pre-draw missing data (don't know how to parallelize this yet)
  propU[0] = unif_rand(); // also pre-draw acceptance Uniforms
  for(II=0; II<nMiss; II++) {
    ii = missInd[II];
    propU[ii] = unif_rand();
    for(jj=nObsComp[ii]; jj<nDims; jj++) {
      propZ[ii*nDims + jj] = norm_rand();
    }
  }
  ii = nComp-1;
  for(jj=nObsComp[ii]; jj<nDims; jj++) {
    propZ[ii*nDims + jj] = norm_rand();
  }
  // intermediate data points
  for(JJ = 0; JJ < 2; JJ++) {
    // *** PARALLELIZABLE FOR-LOOP ***
    #ifdef _OPENMP
#pragma omp parallel for num_threads(nCores) private(ii, jj, II, mean, sd, Z, iCore) if(nCores > 1)
    #endif
    for(II = JJ; II < nMiss; II = II+2) {
      iCore = omp_get_thread_num();
      ii = missInd[II];
      mean = &propMean[iCore*nDims];
      sd = &propSd[iCore*nDims2];
      Z = &propZ[ii*nDims];
      // intermediate data points
      mvEraker(mean, sd,
	       &currX[(ii-1)*nDims], &currX[(ii+1)*nDims],
	       B[ii], sqrtB[ii], currTheta,
	       &sde[iCore]);
      // partial observations
      if(nObsComp[ii] > 0) {
	zmvn(Z, &currX[ii*nDims], mean, sd, nDims, nObsComp[ii]);
      }
      // proposals
      xmvn(&propX[iCore*nDims], Z, mean, sd, nDims);
      // only calculate acceptance rate if proposal is valid
      if(sde[iCore].isValidData(&propX[iCore*nDims], currTheta)) {
	// acceptance rate
	// proposal
	propAccept[iCore] = lmvn(&currX[ii*nDims], Z, mean, sd, nDims);
	propAccept[iCore] -= lmvn(&propX[iCore*nDims], Z, mean, sd, nDims);
	// target 1
	mvEuler(mean, sd, &currX[(ii-1)*nDims],
		dT[ii-1], sqrtDT[ii-1], currTheta, &sde[iCore]);
	propAccept[iCore] += lmvn(&propX[iCore*nDims], Z, mean, sd, nDims);
	propAccept[iCore] -= lmvn(&currX[ii*nDims], Z, mean, sd, nDims);
	// target 2
	mvEuler(mean, sd, &propX[iCore*nDims],
		dT[ii], sqrtDT[ii], currTheta, &sde[iCore]);
	propAccept[iCore] += lmvn(&currX[(ii+1)*nDims], Z, mean, sd, nDims);
	mvEuler(mean, sd, &currX[ii*nDims],
		dT[ii], sqrtDT[ii], currTheta, &sde[iCore]);
	propAccept[iCore] -= lmvn(&currX[(ii+1)*nDims], Z, mean, sd, nDims);
	// evaluate mh ratio
	if(exp(propAccept[iCore]) >= propU[ii]) {
	  for(jj = 0; jj < nDims; jj++) {
	    currX[ii*nDims + jj] = propX[iCore*nDims + jj];
	  }
	  gibbsAccept[ii]++;
	}
      }
    }
  }
  // last datapoint
  if(nMissN > 0) {
    ii = nComp-1;
    iCore = 0;
    mean = &propMean[iCore*nDims];
    sd = &propSd[iCore*nDims2];
    Z = &propZ[ii*nDims];
    mvEuler(mean, sd, &currX[(ii-1)*nDims],
	    dT[ii-1], sqrtDT[ii-1], currTheta, &sde[iCore]);
    // partial observations
    if(nObsComp[ii] > 0) {
      zmvn(Z, &currX[ii*nDims], mean, sd, nDims, nObsComp[ii]);
    }
    // proposals
    xmvn(&propX[iCore*nDims], Z, mean, sd, nDims);
    // acceptance is 100% as long as the proposal is valid
    if(sde[iCore].isValidData(&propX[iCore*nDims], currTheta)) {
      for(jj = 0; jj < nDims; jj++) {
	currX[ii*nDims + jj] = propX[iCore*nDims + jj];
      }
      gibbsAccept[ii]++;
    }
  }
  // first datapoint
  if(nMiss0 > 0) {
    ii = 0;
    iCore = 0;
    mean = &propMean[iCore*nDims];
    sd = &propSd[iCore*nDims2];
    Z = &propZ[iCore*nDims];
    // initialize
    for(jj = 0; jj < nDims; jj++) {
      propX[jj] = currX[jj];
    }
    // random walk metropolis
    for(jj = 0; jj < nMiss0; jj++) {
      // proposal
      propX[nObsComp[0]+jj] = currX[nObsComp[0]+jj] + jumpSd[nParams+jj] * norm_rand();
      if(sde[iCore].isValidData(&propX[iCore*nDims], currTheta)) {
	// acceptance rate.
	// target 1
	propAccept[iCore] = prior->logPrior(currTheta, propX);
	propAccept[iCore] -= prior->logPrior(currTheta, currX);
	// target 2
	mvEuler(mean, sd, &propX[iCore*nDims],
		dT[ii], sqrtDT[ii], currTheta, &sde[iCore]);
	propAccept[iCore] += lmvn(&currX[(ii+1)*nDims], Z, mean, sd, nDims);
	mvEuler(mean, sd, &currX[ii*nDims],
		dT[ii], sqrtDT[ii], currTheta, &sde[iCore]);
	propAccept[iCore] -= lmvn(&currX[(ii+1)*nDims], Z, mean, sd, nDims);
	// evaluate mh ratio
	if(exp(propAccept[iCore]) >= propU[ii]) {
	  currX[nObsComp[0]+jj] = propX[nObsComp[0]+jj];
	  paramAccept[nParams + jj]++;
	}
	else {
	  propX[nObsComp[0]+jj] = currX[nObsComp[0]+jj];
	}
      }
    }
  }
  return;
}


#endif
