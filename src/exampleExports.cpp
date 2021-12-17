// export sample models

//[[Rcpp::depends(RcppProgress)]]

#include <Rcpp.h>
using namespace Rcpp;
#include "sdeRobj.h"
namespace mvn {
#include "mvnPrior.h"
}
namespace hest {
#include "hestModel.h"
}
#undef sdeModel_h
namespace lotvol {
#include "lotvolModel.h"
}
#undef sdeModel_h
namespace biou {
#include "biouModel.h"
}
#undef sdeModel_h
namespace pgnet {
#include "pgnetModel.h"
}
#undef sdeModel_h
namespace eou {
#include "eouModel.h"
}
#undef sdeModel_h

RCPP_MODULE(class_msde_hestModel) {

    class_<sdeRobj<hest::sdeModel, mvn::sdePrior> >("msde_hestModel")

    .constructor()

    .method("nDims", &sdeRobj<hest::sdeModel, mvn::sdePrior>::get_nDims)
    .method("nParams", &sdeRobj<hest::sdeModel, mvn::sdePrior>::get_nParams)
    .method("isData", &sdeRobj<hest::sdeModel, mvn::sdePrior>::isData)
    .method("isParams", &sdeRobj<hest::sdeModel, mvn::sdePrior>::isParams)
    .method("Drift", &sdeRobj<hest::sdeModel, mvn::sdePrior>::Drift)
    .method("Diff", &sdeRobj<hest::sdeModel, mvn::sdePrior>::Diff)
    .method("Loglik", &sdeRobj<hest::sdeModel, mvn::sdePrior>::LogLik)
    .method("Prior", &sdeRobj<hest::sdeModel, mvn::sdePrior>::Prior)
    .method("Sim", &sdeRobj<hest::sdeModel, mvn::sdePrior>::Sim)
    .method("Post", &sdeRobj<hest::sdeModel, mvn::sdePrior>::Post)
    ;
}

RCPP_MODULE(class_msde_lotvolModel) {

    class_<sdeRobj<lotvol::sdeModel, mvn::sdePrior> >("msde_lotvolModel")

    .constructor()

    .method("nDims", &sdeRobj<lotvol::sdeModel, mvn::sdePrior>::get_nDims)
    .method("nParams", &sdeRobj<lotvol::sdeModel, mvn::sdePrior>::get_nParams)
    .method("isData", &sdeRobj<lotvol::sdeModel, mvn::sdePrior>::isData)
    .method("isParams", &sdeRobj<lotvol::sdeModel, mvn::sdePrior>::isParams)
    .method("Drift", &sdeRobj<lotvol::sdeModel, mvn::sdePrior>::Drift)
    .method("Diff", &sdeRobj<lotvol::sdeModel, mvn::sdePrior>::Diff)
    .method("Loglik", &sdeRobj<lotvol::sdeModel, mvn::sdePrior>::LogLik)
    .method("Prior", &sdeRobj<lotvol::sdeModel, mvn::sdePrior>::Prior)
    .method("Sim", &sdeRobj<lotvol::sdeModel, mvn::sdePrior>::Sim)
    .method("Post", &sdeRobj<lotvol::sdeModel, mvn::sdePrior>::Post)
    ;
}

RCPP_MODULE(class_msde_biouModel) {

    class_<sdeRobj<biou::sdeModel, mvn::sdePrior> >("msde_biouModel")

    .constructor()

    .method("nDims", &sdeRobj<biou::sdeModel, mvn::sdePrior>::get_nDims)
    .method("nParams", &sdeRobj<biou::sdeModel, mvn::sdePrior>::get_nParams)
    .method("isData", &sdeRobj<biou::sdeModel, mvn::sdePrior>::isData)
    .method("isParams", &sdeRobj<biou::sdeModel, mvn::sdePrior>::isParams)
    .method("Drift", &sdeRobj<biou::sdeModel, mvn::sdePrior>::Drift)
    .method("Diff", &sdeRobj<biou::sdeModel, mvn::sdePrior>::Diff)
    .method("Loglik", &sdeRobj<biou::sdeModel, mvn::sdePrior>::LogLik)
    .method("Prior", &sdeRobj<biou::sdeModel, mvn::sdePrior>::Prior)
    .method("Sim", &sdeRobj<biou::sdeModel, mvn::sdePrior>::Sim)
    .method("Post", &sdeRobj<biou::sdeModel, mvn::sdePrior>::Post)
    ;
}

RCPP_MODULE(class_msde_pgnetModel) {

    class_<sdeRobj<pgnet::sdeModel, mvn::sdePrior> >("msde_pgnetModel")

    .constructor()

    .method("nDims", &sdeRobj<pgnet::sdeModel, mvn::sdePrior>::get_nDims)
    .method("nParams", &sdeRobj<pgnet::sdeModel, mvn::sdePrior>::get_nParams)
    .method("isData", &sdeRobj<pgnet::sdeModel, mvn::sdePrior>::isData)
    .method("isParams", &sdeRobj<pgnet::sdeModel, mvn::sdePrior>::isParams)
    .method("Drift", &sdeRobj<pgnet::sdeModel, mvn::sdePrior>::Drift)
    .method("Diff", &sdeRobj<pgnet::sdeModel, mvn::sdePrior>::Diff)
    .method("Loglik", &sdeRobj<pgnet::sdeModel, mvn::sdePrior>::LogLik)
    .method("Prior", &sdeRobj<pgnet::sdeModel, mvn::sdePrior>::Prior)
    .method("Sim", &sdeRobj<pgnet::sdeModel, mvn::sdePrior>::Sim)
    .method("Post", &sdeRobj<pgnet::sdeModel, mvn::sdePrior>::Post)
    ;
}

RCPP_MODULE(class_msde_eouModel) {

    class_<sdeRobj<eou::sdeModel, mvn::sdePrior> >("msde_eouModel")

    .constructor()

    .method("nDims", &sdeRobj<eou::sdeModel, mvn::sdePrior>::get_nDims)
    .method("nParams", &sdeRobj<eou::sdeModel, mvn::sdePrior>::get_nParams)
    .method("isData", &sdeRobj<eou::sdeModel, mvn::sdePrior>::isData)
    .method("isParams", &sdeRobj<eou::sdeModel, mvn::sdePrior>::isParams)
    .method("Drift", &sdeRobj<eou::sdeModel, mvn::sdePrior>::Drift)
    .method("Diff", &sdeRobj<eou::sdeModel, mvn::sdePrior>::Diff)
    .method("Loglik", &sdeRobj<eou::sdeModel, mvn::sdePrior>::LogLik)
    .method("Prior", &sdeRobj<eou::sdeModel, mvn::sdePrior>::Prior)
    .method("Sim", &sdeRobj<eou::sdeModel, mvn::sdePrior>::Sim)
    .method("Post", &sdeRobj<eou::sdeModel, mvn::sdePrior>::Post)
    ;
}

// //[[Rcpp::export(".hest_MakeModel")]]
// SEXP hestMakeModel() {
//   sdeCobj *sde = new sdeRobj<hest::sdeModel, mvn::sdePrior>;
//   XPtr<sdeCobj> sdeptr(sde, true);
//   return sdeptr;
// }

// //[[Rcpp::export(".biou_MakeModel")]]
// SEXP biouMakeModel() {
//   sdeCobj *sde = new sdeRobj<biou::sdeModel, mvn::sdePrior>;
//   XPtr<sdeCobj> sdeptr(sde, true);
//   return sdeptr;
// }

// //[[Rcpp::export(".pgnet_MakeModel")]]
// SEXP pgnetMakeModel() {
//   sdeCobj *sde = new sdeRobj<pgnet::sdeModel, mvn::sdePrior>;
//   XPtr<sdeCobj> sdeptr(sde, true);
//   return sdeptr;
// }

// //[[Rcpp::export(".lotvol_MakeModel")]]
// SEXP lotvolMakeModel() {
//   sdeCobj *sde = new sdeRobj<lotvol::sdeModel, mvn::sdePrior>;
//   XPtr<sdeCobj> sdeptr(sde, true);
//   return sdeptr;
// }

// //[[Rcpp::export(".eou_MakeModel")]]
// SEXP eouMakeModel() {
//   sdeCobj *sde = new sdeRobj<eou::sdeModel, mvn::sdePrior>;
//   XPtr<sdeCobj> sdeptr(sde, true);
//   return sdeptr;
// }
