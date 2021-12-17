/// @file rngUtils.h
///
/// Entrypoint for specifying random number generators.
/// Primarily useful for porting C++ library to something other than R.


#ifndef rngUtils_h
#define rngUtils_h


#include <Rmath.h>

namespace sdeRNG {
  /// Generate a uniform random number.
  inline double runif(void) {
    return unif_rand();
  }
  /// Generate a standard normal random number.
  inline double rnorm(void) {
    return norm_rand();
  }
}

/// Fast draw from a (standard) normal distribution.
///
/// Implements a fast version of the Box-Muller algoritm.
///
/// Not currently used, because Box-Muller "counter" issue can't be passed to `RNGkind` on the R side.  However, this makes R-side debugging with `rnorm_fast` extremely difficult.
inline double rnorm_fast(void) {
  static bool doCalc = true;
  static double y1, y2;
  double x1, x2, w, retval;
  if(doCalc) {
    do {
      x1 = 2.0 * unif_rand() - 1.0;
      x2 = 2.0 * unif_rand() - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w;
    retval = y1;
    doCalc = false;
  }
  else {
    retval = y2;
    doCalc = true;
  }
  return retval;
}

#endif
