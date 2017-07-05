/////////////////////////////////////////

#ifndef sdeModel_h
#define sdeModel_h 1

/////////////////////////////////////////

// Heston's model:
// dXt = (alpha - .125 * Zt^2)dt + .5 * Zt dB_Xt
// dZt = (beta/Zt - .5*gamma * Zt)dt + sigma * dB_Zt
// cor(B_Xt, B_Zt) = rho

// sde model object
class sdeModel {
  // put private storage variables here
  // end private storage
 public:
  static const int nParams = 5;
  static const int nDims = 2;
  static const bool sdDiff = true;
  static const bool diagDiff = false;
  void sdeDr(double *dr, double *x, double *theta);
  void sdeDf(double *df, double *x, double *theta);
  bool isValidData(double *x, double *theta);
  bool isValidParams(double *theta);
};

// drift function
inline void sdeModel::sdeDr(double *dr, double *x, double *theta) {
  dr[0] = (theta[0] - .125 * x[1]*x[1]); // x
  dr[1] = (theta[2]/x[1] - .5 * theta[1]*x[1]); // z
  return;
}

// diffusion function
inline void sdeModel::sdeDf(double *df, double *x, double *theta) {
  df[0] = .5 * x[1];
  df[2] = theta[3];
  df[3] = sqrt(1.0-theta[4]*theta[4]) * df[2];
  df[2] *= theta[4];
  return;
}

// data validator
inline bool sdeModel::isValidData(double *x, double *theta) {
  return(x[1] > 0.0);
}

// parameter validator
inline bool sdeModel::isValidParams(double *theta) {
  bool isValid;
  isValid = (theta[1] > 0) && (theta[3] > 0);
  isValid = isValid && (-1.0 < theta[4]) && (1.0 > theta[4]);
  isValid = isValid && (theta[2] > 0.5 * theta[3] * theta[3]);
  return(isValid);
}

#endif
