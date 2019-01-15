/////////////////////////////////////////

#ifndef sdeModel_h
#define sdeModel_h 1

/////////////////////////////////////////

// Prokaryotic Autoregulatory Gene Network model
// of Golightly-Wilkinson 2005

// sde model object
class sdeModel {
private:
  double K, eps; // numerical constants (double so can't be static const)
 public:
  static const int nParams = 8;
  static const int nDims = 4;
  static const bool sdDiff = true;
  static const bool diagDiff = false;
  void sdeDr(double *dr, double *x, double *theta);
  void sdeDf(double *df, double *x, double *theta);
  bool isValidData(double *x, double *theta);
  bool isValidParams(double *theta);
  sdeModel() {K = 10.0; eps = 0.05;}
};

// drift function
inline void sdeModel::sdeDr(double *dr, double *x, double *theta) {
  // double K = 10.0;
  dr[3] = exp(theta[1]) * (K - x[3]) - exp(theta[0]) * x[3] * x[2];
  dr[1] = exp(theta[4]) * x[1] * (x[1]-1.0);
  dr[2] = dr[3] + 0.5 * dr[1];
  dr[0] = exp(theta[5]) * x[2];
  dr[2] = dr[2] - dr[0];
  dr[1] = 2 * dr[0] - dr[1] + exp(theta[3]) * x[0] - exp(theta[7]) * x[1];
  dr[0] = exp(theta[2]) * x[3] - exp(theta[6]) * x[0];
  return;
}

// diffusion function
inline void sdeModel::sdeDf(double *df, double *x, double *theta) {
  // double K = 10.0;
  df[0] = sqrt(exp(theta[2]) * x[3] + exp(theta[6]) * x[0]);
  df[1] = exp(theta[4]) * x[1] * (x[1]-1.0);
  df[2] = exp(theta[5]);
  df[5] = exp(theta[7]) * x[1] + 4*df[2] * x[2] + exp(theta[3]) * x[0] + 2.0*df[1];
  df[14] = exp(theta[0]) * x[3] * x[2] + exp(theta[1]) * (K - x[3]);
  df[9] = -2.0*df[2] * x[2] - df[1];
  df[15] = df[2] * x[2] + df[14] + 0.5*df[1];
  df[10] = df[15] - df[9] * df[9]/ df[5];
  df[15] = sqrt(df[14] - df[14] * df[14] / df[10]);
  df[10] = sqrt(df[10]);
  df[14] = df[14] / df[10];
  df[5] = sqrt(df[5]);
  df[9] = df[9] / df[5];
  df[1] = 0.0;
  df[2] = 0.0;
  df[4] = 0.0;
  df[8] = 0.0;
  df[12] = 0.0;
  df[13] = 0.0;
  return;
}

// data validator
inline bool sdeModel::isValidData(double *x, double *theta) {
  bool isValid;
  // double K = 10.0;
  // double eps = 0.05;
  isValid = x[0] > 1.0+eps;
  isValid = isValid && x[1] > 1.0+eps;
  isValid = isValid && x[2] > 1.0+eps;
  isValid = isValid && ((x[3] > 1.0+eps) && (x[3] < K-eps));
  return(isValid);
}

// parameter validator
inline bool sdeModel::isValidParams(double *theta) {
  return(true);
}

#endif
