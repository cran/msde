// on-the-fly compilation of sde models
// specific to R/C++ interface
// (rest of include directory is not, except rngUtils.h)

//[[Rcpp::depends("msde")]]
//[[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::depends("RcppProgress")]]
#include <RcppArmadillo.h>
using namespace Rcpp;
#include <sdeInterface.h>
#include "sdeModel.h"
#include "sdePrior.h"

//[[Rcpp::export(".sde_MakeModel")]]
SEXP construct() {
  sdeCobj *sde = new sdeRobj<sdeModel, sdePrior>;
  XPtr<sdeCobj> sdeptr(sde, true);
  return sdeptr;
}
