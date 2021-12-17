/// @file biouModel.h

#ifndef biouModel_h
#define biouModel_h 1

/// SDE model class for the bivariate Ornstein-Uhlenbeck model.
///
/// The model is given by
/// ```
/// dYt = (Gamma * Yt + Lambda) * dt + Psi dBt
/// ```
/// where `Yt = (Y1t, Y2t)`, `Gamma` is a `2x2` matrix, `Lambda` is a `2x1` vector, and `Psi` is a `2x2` lower triangular matrix with positive diagonal elements.
/// The data vector is `x = Yt` and the parameter vector is `theta = (Gamma, Lambda, lowertri(Psi)`).

// sde model object
class sdeModel {
 public:
  static const int nParams = 9; ///< Number of SDE parameters.  We have `|Gamma| = 4`, `|Lambda = 2|`, and |Psi = 3|, such that `nParams = 9`.
  static const int nDims = 2; ///< Numbers of SDE dimensions.
  static const bool sdDiff = true; ///< Diffusion is on the standard deviation scale.
  static const bool diagDiff = false; ///< Diffusion is not diagonal.
  /// SDE drift function.
  void sdeDr(double *dr, double *x, double *theta);
  /// SDE diffusion function.
  void sdeDf(double *df, double *x, double *theta);
  /// SDE data validator.
  bool isValidData(double *x, double *theta);
  /// SDE parameter validator.
  bool isValidParams(double *theta);
};

/// @param[out] dr Array into which to store the calculated drift.
/// @param[in] x Array of SDE components at a given time point.
/// @param[in] theta Array of SDE parameters.
inline void sdeModel::sdeDr(double *dr, double *x, double *theta) {
  dr[0] = theta[0]*x[0] + theta[2]*x[1] + theta[4];
  dr[1] = theta[1]*x[0] + theta[3]*x[1] + theta[5];
  return;
}

/// @param[out] df Array into which to store the calculated diffusion matrix.
/// @param[in] x Array of SDE components at a given time point.
/// @param[in] theta Array of SDE parameters.
inline void sdeModel::sdeDf(double *df, double *x, double *theta) {
  df[0] = theta[6];
  df[1] = 0.0;
  df[2] = theta[7];
  df[3] = theta[8];
  return;
}

/// @param[in] x Array of SDE components at a given time point.
/// @param[in] theta Array of SDE parameters.
///
/// @return Whether or not the SDE data `x` is valid.  In this case all data is valid.
inline bool sdeModel::isValidData(double *x, double *theta) {
  return true;
}

/// @param[in] theta Array of SDE parameters.
///
/// @return Whether or not the SDE parameters `theta` are valid.  In this case we must have `Psi[0,0], Psi[0,1] > 0`, i.e., `theta[6], theta[8] > 0`.
inline bool sdeModel::isValidParams(double *theta) {
  return (theta[6] > 0.0) && (theta[8] > 0.0);
}

#endif
