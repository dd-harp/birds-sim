/*
 * bird.h
 * Sean L. Wu (slwood89@gmail.com)
 * May 2022
*/

#ifndef BIRD_H
#define BIRD_H

#include <Rmath.h>
#include <RcppArmadillo.h>

#include <algorithm> // for max_element
#include <vector>

#include "rates.h"

// --------------------------------------------------------------------------------
//   model struct
// --------------------------------------------------------------------------------

template <typename T>
struct bird {
  
  arma::Mat<T> Egg;
  
  arma::Mat<T> Fledgling_S;
  arma::Mat<T> Fledgling_I;
  arma::Mat<T> Fledgling_R;
  
  arma::Row<T> Juvenile_S;
  arma::Row<T> Juvenile_I;
  arma::Row<T> Juvenile_R;
  
  arma::Row<T> Adult_S;
  arma::Row<T> Adult_I;
  arma::Row<T> Adult_R;
  
  arma::Mat<double> psi; // mating dispersal
  arma::Mat<double> theta; // typical home range
  
  arma::SpMat<int> shiftEgg;
  arma::SpMat<int> shiftFledgling;
  
  arma::Mat<double> K; // carrying capacity, patch x time
  double mu_egg;
  double mu_fledge;
  double mu_juvenile;
  double mu_adult;
  
  int step;
  int p;
  double dt;
  
  bird(
    const arma::Mat<double> psi_, arma::Mat<double> theta_,
    const double mu_egg_, const double mu_fledge_,
    const double mu_juvenile_, const double mu_adult_,
    const int delay_egg, const int delay_fledge,
    const int p_, const double dt_
  );
  ~bird() = default;
  
  void update(const Rcpp::List& parameters);
  
};

// constructor
template <typename T>
inline bird<T>::bird(const int p_, const double dt_, const arma::Mat<double>& psi_) :
  psi(psi_), step(0), p(p_), dt(dt_), 
  f(p, arma::fill::zeros), q(n_species, p_, arma::fill::zeros), kappa(n_species, p_, arma::fill::zeros),
  tau_E(tau_E_), tau_L(tau_L_), tau_P(tau_P_), tau_EIP(tau_EIP_)
{
  int maxE = *max_element(tau_E.begin(), tau_E.end());
  int maxL = *max_element(tau_L.begin(), tau_L.end());
  int maxP = *max_element(tau_P.begin(), tau_P.end());
  int maxEIP = *max_element(tau_EIP.begin(), tau_EIP.end());
  
  // state
  E = arma::Mat<T>(maxE, p, arma::fill::zeros);
  L = arma::Mat<T>(maxL, p, arma::fill::zeros);
  P = arma::Mat<T>(maxP, p, arma::fill::zeros);
  
  E_I = arma::Mat<T>(maxE, p, arma::fill::zeros);
  L_I = arma::Mat<T>(maxL, p, arma::fill::zeros);
  P_I = arma::Mat<T>(maxP, p, arma::fill::zeros);
  
  A_S = arma::Row<T>(p, arma::fill::zeros);
  A_E = arma::Mat<T>(maxEIP, p, arma::fill::zeros);
  A_I = arma::Row<T>(p, arma::fill::zeros);
  
  // shift matrices (multiply on left)
  arma::umat fillE(2, maxE);
  for (auto i = 0u; i < (maxE - 1); ++i) {
    fillE(0, i) = i;
    fillE(1, i) = i+1;
  }
  arma::Col<int> fillEvals(maxE, arma::fill::ones);
  shiftE = arma::SpMat<int>(fillE, fillEvals, maxE, maxE);
  
  arma::umat fillL(2, maxL);
  for (auto i = 0u; i < (maxL - 1); ++i) {
    fillL(0, i) = i;
    fillL(1, i) = i+1;
  }
  arma::Col<int> fillLvals(maxL, arma::fill::ones);
  shiftL = arma::SpMat<int>(fillL, fillLvals, maxL, maxL);
  
  arma::umat fillP(2, maxP);
  for (auto i = 0u; i < (maxP - 1); ++i) {
    fillP(0, i) = i;
    fillP(1, i) = i+1;
  }
  arma::Col<int> fillPvals(maxP, arma::fill::ones);
  shiftP = arma::SpMat<int>(fillP, fillPvals, maxP, maxP);
  
  arma::umat fillEIP(2, maxEIP);
  for (auto i = 0u; i < (maxEIP - 1); ++i) {
    fillEIP(0, i) = i;
    fillEIP(1, i) = i+1;
  }
  arma::Col<int> fillEIPvals(maxEIP, arma::fill::ones);
  shiftEIP = arma::SpMat<int>(fillEIP, fillEIPvals, maxEIP, maxEIP);
  
};


#endif