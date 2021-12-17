/// @file {{ModuleFile}}
/// 
/// Rcpp module for template class sdeRobj<{{sdeModel}}, {{sdePrior}}>.

// Rcpp include paths from other libraries
//[[Rcpp::depends(msde)]]
//[[Rcpp::depends(RcppProgress)]]


#include <Rcpp.h>
using namespace Rcpp;
#include <sdeRobj.h>
#include "{{ModelFile}}"
#include "{{PriorFile}}"

RCPP_MODULE({{ModuleName}}) {


  class_<sdeRobj<{{sdeModel}}, {{sdePrior}}> >("{{RClassName}}")

    .constructor()

    .method("nDims", &sdeRobj<{{sdeModel}}, {{sdePrior}}>::get_nDims)
    .method("nParams", &sdeRobj<{{sdeModel}}, {{sdePrior}}>::get_nParams)
    .method("isData", &sdeRobj<{{sdeModel}}, {{sdePrior}}>::isData)
    .method("isParams", &sdeRobj<{{sdeModel}}, {{sdePrior}}>::isParams)
    .method("Drift", &sdeRobj<{{sdeModel}}, {{sdePrior}}>::Drift)
    .method("Diff", &sdeRobj<{{sdeModel}}, {{sdePrior}}>::Diff)
    .method("Loglik", &sdeRobj<{{sdeModel}}, {{sdePrior}}>::LogLik)
    .method("Prior", &sdeRobj<{{sdeModel}}, {{sdePrior}}>::Prior)
    .method("Sim", &sdeRobj<{{sdeModel}}, {{sdePrior}}>::Sim)
    .method("Post", &sdeRobj<{{sdeModel}}, {{sdePrior}}>::Post)
    ;
}
