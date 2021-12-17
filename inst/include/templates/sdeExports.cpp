/// @file {{ExportFile}}
/// 
/// Exported methods for {{sdeModel}}.

using namespace Rcpp;
//[[Rcpp::depends(msde)]]
#include "sdeInterface.h"
#include "{{ModelFile}}"
#include "{{PriorFile}}"

//[[Rcpp::export("{{RClassName}}_ctor")]]
SEXP {{cppClassName}}_ctor() {
  sdeRobj *sde = new sdeRobj<{{sdeModel}}, {{sdePrior}}>;
  XPtr<sdeRobj> sdeptr(sde, true);
  return sdeptr;
}



//[[Rcpp::export("{{RClassName}}_nParams")]]
double {{cppClassName}}_nParams(SEXP sdeptr) {
  XPtr<sdeRobj<{{sdeModel}}, {{sdePrior}}> > sde(sdeptr);
  return sde->get_nParams();
}

//[[Rcpp::export("{{RClassName}}_nDims")]]
double {{cppClassName}}_nDims(SEXP sdeptr) {
  XPtr<sdeRobj<{{sdeModel}}, {{sdePrior}}> > sde(sdeptr);
  return sde->get_nDims();
}

//[[Rcpp::export("{{RClassName}}_isData")]]
LogicalVector {{cppClassName}}_isData(SEXP sdeptr, NumericVector xIn, NumericVector thetaIn,
		     bool singleX, bool singleTheta, int nReps) {
  XPtr<sdeRobj<{{sdeModel}}, {{sdePrior}}> > sde(sdeptr);
  return sde->isData(xIn, thetaIn, singleX, singleTheta, nReps);
}

//[[Rcpp::export("{{RClassName}}_isParams")]]
LogicalVector {{cppClassName}}_isParams(SEXP sdeptr, NumericVector thetaIn, int nReps) {
  XPtr<sdeRobj<{{sdeModel}}, {{sdePrior}}> > sde(sdeptr);
  return sde->isParams(thetaIn, nReps);
}

//[[Rcpp::export("{{RClassName}}_Drift")]]
NumericVector {{cppClassName}}_Drift(SEXP sdeptr, NumericVector xIn, NumericVector thetaIn,
		    bool singleX, bool singleTheta, int nReps) {
  XPtr<sdeRobj<{{sdeModel}}, {{sdePrior}}> > sde(sdeptr);
  return sde->Drift(xIn, thetaIn, singleX, singleTheta, nReps);
}

//[[Rcpp::export("{{RClassName}}_Diff")]]
NumericVector {{cppClassName}}_Diff(SEXP sdeptr, NumericVector xIn, NumericVector thetaIn,
		   bool singleX, bool singleTheta, int nReps) {
  XPtr<sdeRobj<{{sdeModel}}, {{sdePrior}}> > sde(sdeptr);
  return sde->Diff(xIn, thetaIn, singleX, singleTheta, nReps);
}

//[[Rcpp::export("{{RClassName}}_Loglik")]]
NumericVector {{cppClassName}}_LogLik(SEXP sdeptr, NumericVector xIn, NumericVector dTIn,
		     NumericVector thetaIn,
		     int nComp, int nReps,
		     bool singleX, bool singleTheta, int nCores) {
  XPtr<sdeRobj<{{sdeModel}}, {{sdePrior}}> > sde(sdeptr);
  return sde->LogLik(xIn, dTIn, thetaIn, nComp, nReps,
		     singleX, singleTheta, nCores);
}

//[[Rcpp::export("{{RClassName}}_Prior")]]
NumericVector {{cppClassName}}_Prior(SEXP sdeptr, NumericVector thetaIn, NumericVector xIn,
		    bool singleTheta, bool singleX,
		    int nReps, List phiIn) {
  XPtr<sdeRobj<{{sdeModel}}, {{sdePrior}}> > sde(sdeptr);
  return sde->Prior(thetaIn, xIn, singleTheta, singleX, nReps, phiIn);
}

//[[Rcpp::export("{{RClassName}}_Sim")]]
List {{cppClassName}}_Sim(SEXP sdeptr, int nDataOut, int N, int burn, int reps,
	 int r, double dT, int MAXBAD, NumericVector initData,
	 NumericVector params, bool singleX, bool singleTheta) {
  XPtr<sdeRobj<{{sdeModel}}, {{sdePrior}}> > sde(sdeptr);
  return sde->Sim(nDataOut, N, burn, reps, r, dT, MAXBAD,
		  initData, params, singleX, singleTheta);
}

//[[Rcpp::export("{{RClassName}}_Post")]]
List {{cppClassName}}_Post(SEXP sdeptr, NumericVector initParams, NumericVector initData,
	  NumericVector dT, IntegerVector nDimsPerObs,
	  LogicalVector fixedParams, int nSamples, int burn,
	  int nParamsOut, int nDataOut, IntegerVector dataOutSmp,
	  IntegerVector dataOutComp, IntegerVector dataOutDims,
	  double updateParams, double updateData, List priorArgs,
	  List tunePar, int updateLogLik, int nLogLikOut,
	  int updateLastMiss, int nLastMissOut,
	  int nCores, bool displayProgress) {
  XPtr<sdeRobj<{{sdeModel}}, {{sdePrior}}> > sde(sdeptr);
  return sde->Post(initParams, initData, dT, nDimsPerObs, fixedParams,
		   nSamples, burn, nParamsOut, nDataOut, dataOutSmp,
		   dataOutComp, dataOutDims, updateParams, updateData,
		   priorArgs, tunePar, updateLogLik, nLogLikOut,
		   updateLastMiss, nLastMissOut, nCores, displayProgress);
}
