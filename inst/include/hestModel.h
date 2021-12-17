/// @file hestModel.h

#ifndef hestModel_h
#define hestModel_h 1

/// SDE model class for Heston's stochastic volatility model.
///
/// The model is given by
/// ```
/// dXt = (alpha - .125 * Zt^2)dt + .5 * Zt dB_Xt
/// dZt = (beta/Zt - .5*gamma * Zt)dt + sigma * dB_Zt
/// cor(B_Xt, B_Zt) = rho
/// ```
///
/// The data vector is `x = (X, Z)` and the parameter vector is `theta = (alpha, gamma, beta, sigma, rho)`.
class sdeModel {
 public:
  static const int nParams = 5; ///< Number of model parameters.
  static const int nDims = 2; ///< Number of SDE dimensions.
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
  dr[0] = (theta[0] - .125 * x[1]*x[1]); // x
  dr[1] = (theta[2]/x[1] - .5 * theta[1]*x[1]); // z
  return;
}

/// @param[out] df Array into which to store the calculated diffusion matrix.
/// @param[in] x Array of SDE components at a given time point.
/// @param[in] theta Array of SDE parameters.
inline void sdeModel::sdeDf(double *df, double *x, double *theta) {
  df[0] = .5 * x[1];
  df[2] = theta[3];
  df[3] = sqrt(1.0-theta[4]*theta[4]) * df[2];
  df[2] *= theta[4];
  return;
}

/// @param[in] x Array of SDE components at a given time point.
/// @param[in] theta Array of SDE parameters.
///
/// @return Whether or not the SDE data `x` is valid.  In this case we must have `Zt > 0`.
inline bool sdeModel::isValidData(double *x, double *theta) {
  return(x[1] > 0.0);
}

/// @param[in] theta Array of SDE parameters.
///
/// @return Whether or not the SDE parameters `theta` are valid.  In this case we must have `gamma, sigma > 0`, `beta > sigma^2/2`, and `|rho| < 1`.
inline bool sdeModel::isValidParams(double *theta) {
  bool isValid;
  isValid = (theta[1] > 0) && (theta[3] > 0);
  isValid = isValid && (-1.0 < theta[4]) && (1.0 > theta[4]);
  isValid = isValid && (theta[2] > 0.5 * theta[3] * theta[3]);
  return(isValid);
}

#endif
