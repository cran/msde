/////////////////////////////////////////

#ifndef sdeModel_h
#define sdeModel_h 1

/////////////////////////////////////////

// bivariate Ornstein-Uhlenbeck model:
//
// dYt = (Gamma * Yt + Lambda) * dt + Psi dBt,
//
// where Yt = c(Y1_t, Y2_t)
// Gamma is a 2x2 matrix
// Lambda is a 2x1 vector
// Psi is a 2x2 lower triangular matrix with positive diagonal elements.
// the parameter vector is theta = (Gamma, Lambda, lowertri(Psi)).

// sde model object
class sdeModel {
 public:
  static const int nParams = 9; // Gamma = 4, Lambda = 2, Psi = 3
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
  dr[0] = theta[0]*x[0] + theta[2]*x[1] + theta[4];
  dr[1] = theta[1]*x[0] + theta[3]*x[1] + theta[5];
  return;
}

// diffusion function
inline void sdeModel::sdeDf(double *df, double *x, double *theta) {
  df[0] = theta[6];
  df[1] = 0.0;
  df[2] = theta[7];
  df[3] = theta[8];
  return;
}

// data validator
inline bool sdeModel::isValidData(double *x, double *theta) {
  return true;
}

// parameter validator
inline bool sdeModel::isValidParams(double *theta) {
  return (theta[6] > 0.0) && (theta[8] > 0.0);
}

#endif
