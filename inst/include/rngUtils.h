#ifndef rngUtils_h
#define rngUtils_h

// entrypoint for specifying random number generators.
// primarily useful for porting C++ library to something other than R.

#include <Rmath.h>

namespace sdeRNG {
  inline double runif(void) {
    return unif_rand();
  }

  inline double rnorm(void) {
    return norm_rand();
  }
}

// --- fast normal draws -------------------------------------------------------

// not currently used, because box-muller "counter" issue can't be passed to
// RNGkind, which would be ideal for debuggin...
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
