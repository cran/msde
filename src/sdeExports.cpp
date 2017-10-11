// generic R wrappers to exported members of instantiated Cobj

#include <Rcpp.h>
using namespace Rcpp;
#include "sdeInterface.h"

//[[Rcpp::export(".sde_nParams")]]
double sde_nParams(SEXP sdeptr) {
  XPtr<sdeCobj> sde(sdeptr);
  return sde->get_nParams();
}

//[[Rcpp::export(".sde_nDims")]]
double sde_nDims(SEXP sdeptr) {
  XPtr<sdeCobj> sde(sdeptr);
  return sde->get_nDims();
}

//[[Rcpp::export(".sde_isData")]]
LogicalVector isData(SEXP sdeptr, NumericVector xIn, NumericVector thetaIn,
		     bool singleX, bool singleTheta, int nReps) {
  XPtr<sdeCobj> sde(sdeptr);
  return sde->isData(xIn, thetaIn, singleX, singleTheta, nReps);
}

//[[Rcpp::export(".sde_isParams")]]
LogicalVector isParams(SEXP sdeptr, NumericVector thetaIn, int nReps) {
  XPtr<sdeCobj> sde(sdeptr);
  return sde->isParams(thetaIn, nReps);
}

//[[Rcpp::export(".sde_Drift")]]
NumericVector Drift(SEXP sdeptr, NumericVector xIn, NumericVector thetaIn,
		    bool singleX, bool singleTheta, int nReps) {
  XPtr<sdeCobj> sde(sdeptr);
  return sde->Drift(xIn, thetaIn, singleX, singleTheta, nReps);
}

//[[Rcpp::export(".sde_Diff")]]
NumericVector Diff(SEXP sdeptr, NumericVector xIn, NumericVector thetaIn,
		   bool singleX, bool singleTheta, int nReps) {
  XPtr<sdeCobj> sde(sdeptr);
  return sde->Diff(xIn, thetaIn, singleX, singleTheta, nReps);
}

//[[Rcpp::export(".sde_Loglik")]]
NumericVector LogLik(SEXP sdeptr, NumericVector xIn, NumericVector dTIn,
		     NumericVector thetaIn,
		     int nComp, int nReps,
		     bool singleX, bool singleTheta, int nCores) {
  XPtr<sdeCobj> sde(sdeptr);
  return sde->LogLik(xIn, dTIn, thetaIn, nComp, nReps,
		     singleX, singleTheta, nCores);
}

//[[Rcpp::export(".sde_Prior")]]
NumericVector Prior(SEXP sdeptr, NumericVector thetaIn, NumericVector xIn,
		    bool singleTheta, bool singleX,
		    int nReps, List phiIn) {
  XPtr<sdeCobj> sde(sdeptr);
  return sde->Prior(thetaIn, xIn, singleTheta, singleX, nReps, phiIn);
}

//[[Rcpp::export(".sde_Sim")]]
List Sim(SEXP sdeptr, int nDataOut, int N, int burn, int reps,
	 int r, double dT, int MAXBAD, NumericVector initData,
	 NumericVector params, bool singleX, bool singleTheta) {
  XPtr<sdeCobj> sde(sdeptr);
  return sde->Sim(nDataOut, N, burn, reps, r, dT, MAXBAD,
		  initData, params, singleX, singleTheta);
}

//[[Rcpp::export(".sde_Post")]]
List Post(SEXP sdeptr, NumericVector initParams, NumericVector initData,
	  NumericVector dT, IntegerVector nDimsPerObs,
	  LogicalVector fixedParams, int nSamples, int burn,
	  int nParamsOut, int nDataOut, IntegerVector dataOutSmp,
	  IntegerVector dataOutComp, IntegerVector dataOutDims,
	  double updateParams, double updateData, List priorArgs,
	  List tunePar, int updateLogLik, int nLogLikOut,
	  int updateLastMiss, int nLastMissOut,
	  int nCores, bool displayProgress) {
  XPtr<sdeCobj> sde(sdeptr);
  return sde->Post(initParams, initData, dT, nDimsPerObs, fixedParams,
		   nSamples, burn, nParamsOut, nDataOut, dataOutSmp,
		   dataOutComp, dataOutDims, updateParams, updateData,
		   priorArgs, tunePar, updateLogLik, nLogLikOut,
		   updateLastMiss, nLastMissOut, nCores, displayProgress);
}
