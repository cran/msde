#ifndef sdeModel_h
#define sdeModel_h 1

// Lotka-Volterra Predator-Prey model

// class definition
class sdeModel {
 public:
  static const int nParams = 3; // number of model parameters
  static const int nDims = 2; // number of sde dimensions
  static const bool diagDiff = false; // whether diffusion function is diagonal
  static const bool sdDiff = true; // whether diffusion is on sd or var scale
  void sdeDr(double *dr, double *x, double *theta); // drift function
  void sdeDf(double *df, double *x, double *theta); // diffusion function
  bool isValidParams(double *theta); // parameter validator
  bool isValidData(double *x, double *theta); // data validator
};

// drift function
inline void sdeModel::sdeDr(double *dr, double *x, double *theta) {
  dr[0] = theta[0]*x[0] - theta[1]*x[0]*x[1]; // alpha * H - beta * H*L
  dr[1] = theta[1]*x[0]*x[1] - theta[2]*x[1]; // beta * H*L - gamma * L
  return;
}

// diffusion function (sd scale)
inline void sdeModel::sdeDf(double *df, double *x, double *theta) {
  double bHL = theta[1]*x[0]*x[1]; // beta * H*L
  df[0] = sqrt(theta[0]*x[0] + bHL); // sqrt(alpha * H + bHL)
  df[2] = -bHL/df[0];
  df[3] = sqrt(bHL + theta[2]*x[1] - df[2]*df[2]);
  return;
}

// parameter validator
inline bool sdeModel::isValidParams(double *theta) {
  bool val = theta[0] > 0.0;
  val = val && theta[1] > 0.0;
  val = val && theta[2] > 0.0;
  return val;
}

// data validator
inline bool sdeModel::isValidData(double *x, double *theta) {
  return (x[0] > 0.0) && (x[1] > 0.0);
}

#endif
