/// @file sdeData.h

#ifndef sdeData_h
#define sdeData_h 1

/// Hold data for an SDE model.
///
/// Creates space for the data itself, timepoints, and storage space for mean/variance and normal likelihood evaluations.
///
/// @note A lot of public members here should in fact be protected/private.  Also should use standard setter/getter for extraction.
template <class sMod>
class sdeData {
 private:
  /// Allocate dynamic memory (used by constructors).
  void ctor_init(int ncomp, int nmv, int nz, int nsde);
 protected:
  int nDims2; ///< Number of variance dimensions.
 public:
  int nDims; ///< Number of dimensions (get from sMod).
  int nParams; ///< Number of parameters (get from sMod).
  int nComp; ///< Number of (complete data) observations.
  double *dT; ///< Array of interobservation times (length nComp-1).
  double *sqrtDT; ///< Array of sqrt(dT).
  // double *XComp; ///< Complete data.
  int *nObsComp; ///< Number of observed dimensions per observation.
  double *propMean; ///< Storage for SDE means.
  double *propSd; ///< Storage for SDE variances (in sd form).
  sMod *sde; ///< Storage for drift/diffusion calculations.
  double *propZ; ///< Storage for normal likelihoods.
  /// Constructor.
  sdeData(int ncomp, double *dt, int nmv, int nz, int nsde);
  /// Constructor with missing data information.
  sdeData(int ncomp, double *dt, int *nobscomp, int nmv, int nz, int nsde);
  /// Default constructor.
  sdeData();
  /// Destructor.
  ~sdeData();
};


/// @param[in] ncomp Number of complete data observations.
/// @param[in] nmv Number of mean/variance units to store.
/// @param[in] nz Number of additional `z` vectors to store.
/// @param[in] nsde Number of `sMod` objects to store.
template <class sMod>
inline void sdeData<sMod>::ctor_init(int ncomp, int nmv, int nz, int nsde) {
  nDims = sMod::nDims;
  nDims2 = sMod::diagDiff ? nDims : nDims*nDims;
  nParams = sMod::nParams;
  nComp = ncomp;
  // create storage space
  dT = new double[nComp];
  sqrtDT = new double[nComp];
  // XComp = new double[nComp*nDims];
  propMean = new double[nmv*nDims];
  propSd = new double[nmv*nDims2];
  propZ = new double[nz*nDims];
  sde = new sMod[nsde];
  nObsComp = new int[nComp];
  return;
}

/// @param[in] ncomp Number of complete data observations.
/// @param[in] dt Array of interobservation times.
/// @param[in] nmv Number of mean/variance units to store.
/// @param[in] nz Number of additional z vectors to store.
/// @param[in] nsde Number of sdeModel objects to store.
template <class sMod>
inline sdeData<sMod>::sdeData(int ncomp, double *dt,
			      int nmv, int nz, int nsde) {
  ctor_init(ncomp, nmv, nz, nsde);
  // assignments
  for(int ii=0; ii<nComp-1; ii++) {
    dT[ii] = dt[ii];
    sqrtDT[ii] = sqrt(dT[ii]);
  }
}

/// @param[in] ncomp Number of complete data observations.
/// @param[in] dt Array of interobservation times.
/// @param[in] nobscomp Number of observed sde components per time point.
/// @param[in] nmv Number of mean/variance units to store.
/// @param[in] nz Number of additional `z` vectors to store.
/// @param[in] nsde Number of `sMod` objects to store.
template <class sMod>
inline sdeData<sMod>::sdeData(int ncomp, double *dt, int *nobscomp, 
			      int nmv, int nz, int nsde) {
  ctor_init(ncomp, nmv, nz, nsde);
  // assignments
  for(int ii=0; ii<nComp-1; ii++) {
    dT[ii] = dt[ii];
    sqrtDT[ii] = sqrt(dT[ii]);
    nObsComp[ii] = nobscomp[ii];
  }
  nObsComp[nComp-1] = nobscomp[nComp-1];
  // // assign partial observations
  // for(int ii=0; ii<nComp; ii++) {
  //   nObsComp[ii] = nobscomp[ii];
  // }
}

/// The default constructor allocates a single unit to each memory variable. 
template <class sMod>
inline sdeData<sMod>::sdeData() {
  ctor_init(1, 1, 1, 1);
}

template <class sMod>
inline sdeData<sMod>::~sdeData() {
  delete [] nObsComp;
  delete [] sde;
  delete [] propMean;
  delete [] propSd;
  delete [] propZ;
  delete [] dT;
  delete [] sqrtDT;
}

#endif
