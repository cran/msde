// export sample models

#include <Rcpp.h>
using namespace Rcpp;
#include "sdeInterface.h"
namespace mvn {
#include "mvnPrior.h"
#undef sdePrior_h
}
namespace hest {
#include "hestModel.h"
#undef sdeModel_h
}
namespace lotvol {
#include "lotvolModel.h"
#undef sdeModel_h
}
namespace biou {
#include "biouModel.h"
#undef sdeModel_h
}
namespace pgnet {
#include "pgnetModel.h"
#undef sdeModel_h
}

//[[Rcpp::export(".hest_MakeModel")]]
SEXP hestMakeModel() {
  sdeCobj *sde = new sdeRobj<hest::sdeModel, mvn::sdePrior>;
  XPtr<sdeCobj> sdeptr(sde, true);
  return sdeptr;
}

//[[Rcpp::export(".biou_MakeModel")]]
SEXP biouMakeModel() {
  sdeCobj *sde = new sdeRobj<biou::sdeModel, mvn::sdePrior>;
  XPtr<sdeCobj> sdeptr(sde, true);
  return sdeptr;
}

//[[Rcpp::export(".pgnet_MakeModel")]]
SEXP pgnetMakeModel() {
  sdeCobj *sde = new sdeRobj<pgnet::sdeModel, mvn::sdePrior>;
  XPtr<sdeCobj> sdeptr(sde, true);
  return sdeptr;
}

//[[Rcpp::export(".lotvol_MakeModel")]]
SEXP lotvolMakeModel() {
  sdeCobj *sde = new sdeRobj<lotvol::sdeModel, mvn::sdePrior>;
  XPtr<sdeCobj> sdeptr(sde, true);
  return sdeptr;
}
