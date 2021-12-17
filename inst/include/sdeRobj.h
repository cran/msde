/// @file sdeRobj.h
///
/// R/C++ class interface.
///
/// `sdeRobj` is a C++ class of which all methods have Rcpp inputs and outputs.  It can then be wrapped on the R side using Rcpp modules.

#ifndef sdeRobj_h
#define sdeRobj_h

#include <Rcpp.h>
typedef Rcpp::LogicalVector Logical;
typedef Rcpp::NumericVector Numeric;
typedef Rcpp::IntegerVector Integer;
typedef Rcpp::NumericMatrix NumericMatrix;
typedef Rcpp::List List;
//using namespace Rcpp;
//#include "sdeLogLik.h"
//#include "sdeMCMC.h"
//#include "mcmcUtils.h"

// --- helper functions --------------------------------------------------------

/// Convert `Rcpp::LogicalVector` to a `bool*`.
///
/// `Rcpp::LogicalVector` is actually an `int*` on the C++ side, in order to handle `true/false/NA`.
///
/// @param[out] out Pointer to C++ logical vector.
/// @param[in] in Rcpp logical vector.
inline void convert_Logical(bool *out, Logical in) {
  for(int ii=0; ii<in.length(); ii++) {
    if(Logical::is_na(in[ii])) {
      Rprintf("convert_Logical: NA detected.\n");
    }
    out[ii] = (in[ii] != 0);
  }
  return;
}

/// Class to parse prior arguments from an `Rcpp::List`.
class PriorArgs {
 public:
  int nArgs; ///< Number of prior arguments.
  double **phi; ///< Pointer to storage for each argument.
  int *nEachArg; ///< Pointer to length of each argument.
  /// Class constructor.
  PriorArgs(List phiIn);
  /// Class destructor.
  ~PriorArgs() {
    delete [] nEachArg;
    delete [] phi;
  }
};

/// Each element of the input list is assumed to be an Rcpp::NumericVector`.  The constructor copies these into its `phi` member as a vector of `double*` (i.e., double pointer).  The other (public) members are `nArgs`, the number of arguments and `nEachArg` an integer vector of argument lengths.
///
/// @param[in] phiIn A list of prior arguments. 
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

// --- sdeRobj -----------------------------------------------------------------

/// C++ class of which methods provide Rcpp entry/exit points for given SDE model methods.
template <class sMod, class sPi>
class sdeRobj {
public:
  /// Get number of parameters of SDE model.
  int get_nParams(void);
  /// Get number of dimensions of SDE model.
  int get_nDims(void);
  /// Determine whether SDE data is valid.
  Logical isData(Numeric xIn, Numeric thetaIn,
		 bool singleX, bool singleTheta, int nReps);
  /// Determine whethere SDE parameters are valid.
  Logical isParams(Numeric thetaIn, int nReps);
  /// SDE drift function.
  Numeric Drift(Numeric xIn, Numeric thetaIn,
		bool singleX, bool singleTheta, int nReps);
  /// SDE diffusion function.
  Numeric Diff(Numeric xIn, Numeric thetaIn,
	       bool singleX, bool singleTheta, int nReps);
  /// Calculate the SDE model loglikelihood.
  Numeric LogLik(Numeric xIn, Numeric dTIn,
		 Numeric thetaIn,
		 int nComp, int nReps,
		 bool singleX, bool singleTheta, int nCores);
  /// Calculate the log prior.
  Numeric Prior(Numeric thetaIn, Numeric xIn,
		bool singleTheta, bool singleX,
		int nReps, List phiIn);
  /// Simulate data from the SDE model.
  List Sim(int nDataOut, int N, int burn, int reps, int r, double dT,
	   int MAXBAD, Numeric initData, Numeric params,
	   bool singleX, bool singleTheta);
  /// MCMC sampling from the SDE model posterior distribution.
  List Post(Numeric initParams, Numeric initData,
	    Numeric dT, Integer nDimsPerObs,
	    Logical fixedParams, int nSamples, int burn,
	    int nParamsOut, int nDataOut, Integer dataOutSmp,
	    Integer dataOutComp, Integer dataOutDims,
	    double updateParams, double updateData, List priorArgs,
	    List tunePar, int updateLogLik, int nLogLikOut,
	    int updateLastMiss, int nLastMissOut,
	    int nCores, bool displayProgress);
};

#include "sdeRUtils.h"
#include "sdeSim.h"
#include "sdePost.h"

#endif
