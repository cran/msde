/*

OK interface considerations for the implementation.

1. sde.make.model should create an Xptr instead of returning R/C++ entrypoints
   for all functions.  Therefore, it should trigger the construction of an
   sdeInterface object of which the members interface between R/C++ for the
   necessary functions: drift, diff, sim, post, etc.
2. must be able to hold different models in one R session.  This means that Xptr
   should be a templated class.  That is, each class/function making use of sdeModel and sdePrior should accept templates of these.
3. should have pre-compiled versions of some models.  For this, the src of the package should have pointers to different models in its shared object...

OK:

Abstract base, template derived class
sdeCobj : sdeRobj<sdeModel, sdePrior>

R constructor:
sdeCobj *sde = new sdeRobj<sdeModel, sdePrior>;
XPtr<sdeRobj> sdeptr(sde, true);
return sdeptr;

R methods:
sde.drift/diff
sde.loglik
sde.prior
sde.valid.params/data
sde.sim
sde.post

*/

#ifndef sdeInterface_h
#define sdeInterface_h

#include <Rcpp.h>
typedef Rcpp::LogicalVector Logical;
typedef Rcpp::NumericVector Numeric;
typedef Rcpp::IntegerVector Integer;
using Rcpp::List;
//typedef Rcpp::List List;
//using namespace Rcpp;
//#include "sdeLogLik.h"
//#include "sdeMCMC.h"
//#include "mcmcUtils.h"

// --- helper functions --------------------------------------------------------

// R -> C++ logical vector conversion.
inline void convert_Logical(bool *out, Logical in) {
  for(int ii=0; ii<in.length(); ii++) {
    if(Logical::is_na(in[ii])) {
      Rprintf("convert_Logical: NA detected.\n");
    }
    out[ii] = (in[ii] != 0);
  }
  return;
}

// parse prior arguments from R list
class PriorArgs {
 public:
  int nArgs;
  double **phi;
  int *nEachArg;
  PriorArgs(List phiIn);
  ~PriorArgs() {
    delete [] nEachArg;
    delete [] phi;
  }
};

inline PriorArgs::PriorArgs(List phiIn) {
  nArgs = phiIn.length(); // at least 1, since hyper at least list(NULL)
  phi = new double*[nArgs];
  nEachArg = new int[nArgs];
  for(int ii=0; ii<nArgs; ii++) {
    if(Rf_isNull(phiIn[ii])) {
      nEachArg[ii] = 0;
    } else {
      nEachArg[ii] = as<Numeric>(phiIn[ii]).length();
      phi[ii] = REAL(phiIn[ii]);
    }
  }
}

// --- C -> R object conversion ------------------------------------------------

// for some reason sdeCobj needs a virtual destructor with explicit default.
class sdeCobj {
 public:
  virtual int get_nParams(void) = 0;
  virtual int get_nDims(void) = 0;
  virtual Logical isData(Numeric xIn, Numeric thetaIn,
			       bool singleX, bool singleTheta, int nReps) = 0;
  virtual Logical isParams(Numeric thetaIn, int nReps) = 0;
  virtual Numeric Drift(Numeric xIn, Numeric thetaIn,
			      bool singleX, bool singleTheta, int nReps) = 0;
  virtual Numeric Diff(Numeric xIn, Numeric thetaIn,
			     bool singleX, bool singleTheta, int nReps) = 0;
  virtual Numeric LogLik(Numeric xIn, Numeric dTIn,
			       Numeric thetaIn,
			       int nComp, int nReps,
			       bool singleX, bool singleTheta, int nCores) = 0;
  virtual Numeric Prior(Numeric thetaIn, Numeric xIn,
			      bool singleTheta, bool singleX,
			      int nReps, List phiIn) = 0;
  virtual List Sim(int nDataOut, int N, int burn, int reps, int r, double dT,
		   int MAXBAD, Numeric initData, Numeric params,
		   bool singleX, bool singleTheta) = 0;
  virtual List Post(Numeric initParams, Numeric initData,
		    Numeric dT, Integer nDimsPerObs,
		    Logical fixedParams, int nSamples, int burn,
		    int nParamsOut, int nDataOut, Integer dataOutSmp,
		    Integer dataOutComp, Integer dataOutDims,
		    double updateParams, double updateData, List priorArgs,
		    List tunePar, int updateLogLik, int nLogLikOut,
		    int updateLastMiss, int nLastMissOut,
		    int nCores, bool displayProgress) = 0;
  virtual ~sdeCobj() = 0;
};

// default destructor
inline sdeCobj::~sdeCobj() {
  //Rprintf("sdeCobj destroyed.\n");
}

template <class sMod, class sPi>
class sdeRobj : public sdeCobj {
 public:
  virtual int get_nParams(void);
  virtual int get_nDims(void);
  virtual Logical isData(Numeric xIn, Numeric thetaIn,
			       bool singleX, bool singleTheta, int nReps);
  virtual Logical isParams(Numeric thetaIn, int nReps);
  virtual Numeric Drift(Numeric xIn, Numeric thetaIn,
			      bool singleX, bool singleTheta, int nReps);
  virtual Numeric Diff(Numeric xIn, Numeric thetaIn,
			     bool singleX, bool singleTheta, int nReps);
  virtual Numeric LogLik(Numeric xIn, Numeric dTIn,
			       Numeric thetaIn,
			       int nComp, int nReps,
			       bool singleX, bool singleTheta, int nCores);
  virtual Numeric Prior(Numeric thetaIn, Numeric xIn,
			      bool singleTheta, bool singleX,
			      int nReps, List phiIn);
  virtual List Sim(int nDataOut, int N, int burn, int reps, int r, double dT,
		   int MAXBAD, Numeric initData, Numeric params,
		   bool singleX, bool singleTheta);
  virtual List Post(Numeric initParams, Numeric initData,
		    Numeric dT, Integer nDimsPerObs,
		    Logical fixedParams, int nSamples, int burn,
		    int nParamsOut, int nDataOut, Integer dataOutSmp,
		    Integer dataOutComp, Integer dataOutDims,
		    double updateParams, double updateData, List priorArgs,
		    List tunePar, int updateLogLik, int nLogLikOut,
		    int updateLastMiss, int nLastMissOut,
		    int nCores, bool displayProgress);
  virtual ~sdeRobj() {
    //Rprintf("sdeRobj destroyed.\n");
  };
};

#include "sdeRUtils.h"
#include "sdeSim.h"
#include "sdePost.h"

#endif
