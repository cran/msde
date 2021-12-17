/// @file lotvolModel.h

#ifndef lotvolModel_h
#define lotvolModel_h 1

/// SDE model class for the Lotka-Volterra predator-prey model.
///
/// The model is given by
/// \f[
/// \begin{bmatrix} \mathrm{d} H_t \\ \mathrm{d} L_t \end{bmatrix} = \begin{bmatrix} \alpha H_t - \beta H_tL_t \\ \beta H_tL_t - \gamma L_t \end{bmatrix}\, \mathrm{d} t + \begin{bmatrix} \alpha H_t + \beta H_tL_t & -\beta H_tL_t \\ -\beta H_tL_t & \beta H_tL_t + \gamma L_t\end{bmatrix}^{1/2} \begin{bmatrix} \mathrm{d} B_{1t} \\ \mathrm{d} B_{2t} \end{bmatrix}.
/// \f]
/// The data vector is `x = (H, L)` and the parameter vector is `theta = (alpha, beta, gamma)`.
class sdeModel {
 public:
  static const int nParams = 3; ///< Number of model parameters.
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
  dr[0] = theta[0]*x[0] - theta[1]*x[0]*x[1]; // alpha * H - beta * H*L
  dr[1] = theta[1]*x[0]*x[1] - theta[2]*x[1]; // beta * H*L - gamma * L
  return;
}

/// @param[out] df Array into which to store the calculated diffusion matrix.
/// @param[in] x Array of SDE components at a given time point.
/// @param[in] theta Array of SDE parameters.
inline void sdeModel::sdeDf(double *df, double *x, double *theta) {
  double bHL = theta[1]*x[0]*x[1]; // beta * H*L
  df[0] = sqrt(theta[0]*x[0] + bHL); // sqrt(alpha * H + bHL)
  df[2] = -bHL/df[0];
  df[3] = sqrt(bHL + theta[2]*x[1] - df[2]*df[2]);
  return;
}

/// @param[in] x Array of SDE components at a given time point.
/// @param[in] theta Array of SDE parameters.
///
/// @return Whether or not the SDE data `x` is valid.  In this case we must have `Ht > 0` and `Lt > 0`.
inline bool sdeModel::isValidData(double *x, double *theta) {
  return (x[0] > 0.0) && (x[1] > 0.0);
}


/// @param[in] theta Array of SDE parameters.
///
/// @return Whether or not the SDE parameters `theta` are valid.  In this case we must have `alpha, beta, gamma > 0`.
inline bool sdeModel::isValidParams(double *theta) {
  bool val = theta[0] > 0.0;
  val = val && theta[1] > 0.0;
  val = val && theta[2] > 0.0;
  return val;
}


#endif
