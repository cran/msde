#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::depends(msde)]]
//[[Rcpp::depends(RcppProgress)]]
#include <sdeInterface.h>
#include "sdeModel.h"
#include "sdePrior.h"

//[[Rcpp::export(".sde_MakeModel")]]
SEXP construct() {
  sdeCobj *sde = new sdeRobj<sdeModel, sdePrior>;
  XPtr<sdeCobj> sdeptr(sde, true);
  return sdeptr;
}
