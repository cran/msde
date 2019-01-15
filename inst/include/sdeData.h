#ifndef sdeData_h
#define sdeData_h 1

/// Hold data for an SDE model.
///
/// Creates space for the data itself, timepoints,
/// and storage space for mean/variance and normal likelihood evaluations.
///
/// \note A lot of public members here should in fact be protected/private.  Also should use standard setter/getter for extraction.

template <class sMod>
class sdeData {
 private:
  /// allocate dynamic memory (used by constructors)
  void ctor_init(int ncomp, int nmv, int nz, int nsde);
 protected:
  int nDims2; ///< number of variance dimensions
 public:
  int nDims; ///< number of dimensions (get from sMod)
  int nParams; ///< number of parameters (get from sMod)
  int nComp; ///< number of (complete data) observations
  double *dT; ///< array of interobservation times (length nComp-1)
  double *sqrtDT; ///< array of sqrt(dT)
  // double *XComp; ///< complete data
  int *nObsComp; ///< number of observed dimensions per observation
  double *propMean; ///< storage for sde means
  double *propSd; ///< storage for sde variances (in sd form)
  sMod *sde; ///< storage for drift/diffusion calculations
  double *propZ; ///< storage for normal likelihoods
  /// constructor
  sdeData(int ncomp, double *dt, int nmv, int nz, int nsde);
  sdeData(int ncomp, double *dt, int *nobscomp, int nmv, int nz, int nsde);
  /// default constructor
  sdeData();
  /// destructor
  ~sdeData();
};


/// \param nc number of (complete data observations)
/// \param dt array of interobservation times
/// \param nmv number of mean/variance units to store
/// \param nz number of additional z vectors to store
/// \param nsde number of sdeModel objects to store
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

/// \param nc number of (complete data observations)
/// \param dt array of interobservation times
/// \param nmv number of mean/variance units to store
/// \param nz number of additional z vectors to store
/// \param nsde number of sdeModel objects to store
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

/// \param nc number of (complete data observations)
/// \param dt array of interobservation times
/// \param nobscomp number of observed sde components per time point
/// \param nmv number of mean/variance units to store
/// \param nz number of additional z vectors to store
/// \param nsde number of sdeModel objects to store
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

/// Default constructor allocates a single unit to each memory variable. 
template <class sMod>
inline sdeData<sMod>::sdeData() {
  ctor_init(1, 1, 1, 1);
}

template <class sMod>
inline sdeData<sMod>::~sdeData() {
  // delete [] XComp;
  delete [] nObsComp;
  delete [] sde;
  delete [] propMean;
  delete [] propSd;
  delete [] propZ;
  delete [] dT;
  delete [] sqrtDT;
}

#endif
