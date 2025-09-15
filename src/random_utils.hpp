#ifndef RANDOM_UTILS_H
#define RANDOM_UTILS_H

#include <Rcpp.h>

// Use R's random number generators directly - avoid namespace conflicts
inline double runif(double min, double max) {
  return ::Rf_runif(min, max);
}

inline double rnorm(double mean, double sd) {
  return ::Rf_rnorm(mean, sd);
}

inline double rexp(double rate) {
  return ::Rf_rexp(rate);
}

#endif // RANDOM_UTILS_H