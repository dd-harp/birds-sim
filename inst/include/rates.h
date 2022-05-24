/*
 * rates.h
 * Sean L. Wu (slwood89@gmail.com)
 * May 2022
 */

#ifndef RATES_H
#define RATES_H

#include <RcppArmadillo.h>

inline double mating(const double t, const Rcpp::List& pars) {
  double kB = pars["kB"];
  double lambdaB = pars["lambdaB"];
  double phiB = pars["phiB"];
  
  return 0.5 * kB * exp(-lambdaB * pow(cos((M_PI * t) + phiB / 365.0), 2.0));
}

#endif