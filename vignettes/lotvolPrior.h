#ifndef sdePrior_h
#define sdePrior_h 1

#include <Rcpp.h> // contains R's dlnorm function

// Prior for Lotka-Volterra Model

class sdePrior {
  private:
  static const int nHyper = 4; // (alpha, beta, gamma, L)
  double *mean, *sd; // log-normal mean and standard deviation vectors
 public:
  double logPrior(double *theta, double *x); // log-prior function
  sdePrior(double **phi, int nArgs, int *nEachArg); // constructor
  ~sdePrior(); // destructor
};

// constructor
inline sdePrior::sdePrior(double **phi, int nArgs, int *nEachArg) {
  // allocate memory for hyperparameters
  mean = new double[nHyper];
  sd = new double[nHyper];
  // hard-copy hyperparameters into Prior object
  for(int ii=0; ii<nHyper; ii++) {
    mean[ii] = phi[0][ii];
    sd[ii] = phi[1][ii];
  }
}

// destructor
inline sdePrior::~sdePrior() {
  // deallocate memory to avoid memory leaks
  delete [] mean;
  delete [] sd;
}

// log-prior function itself:
// independent log-normal densities for (alpha,beta,gamma,L)
inline double sdePrior::logPrior(double *theta, double *x) {
  double lpi = 0.0;
  // alpha,beta,gamma
  for(int ii=0; ii < 3; ii++) {
    lpi += R::dlnorm(theta[ii], mean[ii], sd[ii], 1);
  }
  // L
  lpi += R::dlnorm(x[1], mean[3], sd[3], 1);
  return(lpi);
}

#endif
