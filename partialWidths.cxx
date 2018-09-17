// Copyright (C) 2018 - Fabian A.J. Thiele, <fabian.thiele@cern.ch>

#include "Lepton.h"
#include "Meson.h"
#include "HNL.h"
#include "Config.h"
#include <iostream>
#include <map>
#include <vector>
#include "auxfunctions.h"
#include "partialWidths.h"
#include "Logger.h"

Double_t pw_nualpha_lbeta_lbeta(std::shared_ptr<Config> cfg, const Lepton &alpha, const Lepton &beta, const HNL &N) {
  if(not N.mixesWith(alpha)) return 0;
  if(N.getMass() < 2.*beta.getMass()) return 0; // this means we don't have enough mass in the HNL to produce decay product on-shell

  mpfr_t fermiC, fermiCsq, pi, VUDsq, SOL, HBAR;
  unsigned int BITS = cfg->getBITS();
  mpfr_init2(fermiC, BITS);
  mpfr_init2(fermiCsq, BITS);
  mpfr_init2(pi, BITS);
  mpfr_init2(VUDsq, BITS);
  mpfr_init2(SOL, BITS);
  mpfr_init2(HBAR, BITS);


  cfg->getFermiCsq(fermiCsq);
  cfg->getFermiC(fermiC);
  cfg->getPi(pi);
  cfg->getVUDsq(VUDsq);
  cfg->getSOL(SOL);
  cfg->getHBAR(HBAR);

  mpfr_t denominator, temp, x, factor, betamass, HNLmass, NZ, angle;

  mpfr_init2(temp, BITS);
  mpfr_init2(x, BITS);
  mpfr_init2(factor, BITS);
  mpfr_init2(betamass, BITS);
  mpfr_init2(HNLmass, BITS);
  mpfr_init2(NZ, BITS);
  mpfr_init2(angle, BITS);
  mpfr_init2(denominator, BITS);

  mpfr_set_d(NZ, 1, MPFR_RNDD);
  mpfr_set_d(betamass, beta.getMass(), MPFR_RNDD);
  mpfr_set_d(HNLmass, N.getMass(), MPFR_RNDD);
  mpfr_set_d(angle, N.getAngle(), MPFR_RNDD);
  mpfr_div(x, betamass, HNLmass, MPFR_RNDD);

  // calculate first factors NZ * GF^2 * M^5 * U^2/(192*pi^3)
  mpfr_mul(factor, NZ, fermiCsq, MPFR_RNDD);
  mpfr_pow_ui(temp, HNLmass, 5, MPFR_RNDD);
  mpfr_mul(factor, factor, temp, MPFR_RNDD);
  mpfr_pow_ui(denominator, pi, 3, MPFR_RNDD);
  mpfr_mul_ui(denominator, denominator, 192, MPFR_RNDD);
  mpfr_div(factor, factor, denominator, MPFR_RNDD);
  mpfr_mul(factor, factor, angle, MPFR_RNDD);
  // end calculation first factors

  // pre-calculate square root of 1-4*x^2
  mpfr_t sq1m4x2;
  mpfr_init2(sq1m4x2, BITS);
  mpfr_pow_ui(sq1m4x2, x, 2, MPFR_RNDD);
  mpfr_mul_ui(sq1m4x2, sq1m4x2, 4, MPFR_RNDD);
  mpfr_ui_sub(sq1m4x2, 1, sq1m4x2, MPFR_RNDD);
  mpfr_sqrt(sq1m4x2, sq1m4x2, MPFR_RNDD);
  // end pre-calculation

  // calculate L(x) (logarithm function)
  mpfr_t L;
  mpfr_init2(L, BITS);
  mpfr_pow_ui(L, x, 2, MPFR_RNDD);
  mpfr_mul_ui(L, L, 3, MPFR_RNDD);
  mpfr_ui_sub(L, 1, L, MPFR_RNDD);
  mpfr_pow_ui(temp, x, 2, MPFR_RNDD);
  mpfr_ui_sub(temp, 1, temp, MPFR_RNDD);
  mpfr_mul(temp, temp, sq1m4x2, MPFR_RNDD);
  mpfr_sub(L, L, temp, MPFR_RNDD);
  mpfr_pow_ui(denominator, x, 2, MPFR_RNDD);
  mpfr_add_ui(temp, sq1m4x2, 1, MPFR_RNDD);
  mpfr_mul(denominator, denominator, temp, MPFR_RNDD);
  mpfr_div(L, L, denominator, MPFR_RNDD);
  mpfr_log(L, L, MPFR_RNDD);
  // end calculation L(x)

  mpfr_t weinberg;
  mpfr_init2(weinberg, BITS);
  mpfr_set_d(weinberg, 0.2223, MPFR_RNDD);

  mpfr_t c1, c2, c1_factor, c2_factor;
  mpfr_init2(c1, BITS);
  mpfr_init2(c2, BITS);
  mpfr_init2(c1_factor, BITS);
  mpfr_init2(c2_factor, BITS);

  if(alpha==beta) {
    mpfr_pow_ui(c1, weinberg, 2, MPFR_RNDD);
    mpfr_mul_ui(c1, c1, 8, MPFR_RNDD);
    mpfr_mul_ui(temp, weinberg, 4, MPFR_RNDD);
    mpfr_add(c1, c1, temp, MPFR_RNDD);
    mpfr_add_ui(c1, c1, 1, MPFR_RNDD);
    mpfr_div_ui(c1, c1, 4, MPFR_RNDD);

    mpfr_mul_ui(c2, weinberg, 2, MPFR_RNDD);
    mpfr_add_ui(c2, c2, 1, MPFR_RNDD);
    mpfr_mul(c2, weinberg, c2, MPFR_RNDD);
    mpfr_div_ui(c2, c2, 2, MPFR_RNDD);
  } else {
    mpfr_pow_ui(c1, weinberg, 2, MPFR_RNDD);
    mpfr_mul_ui(c1, c1, 8, MPFR_RNDD);
    mpfr_mul_ui(temp, weinberg, 4, MPFR_RNDD);
    mpfr_ui_sub(temp, 1, temp, MPFR_RNDD);
    mpfr_add(c1, c1, temp, MPFR_RNDD);
    mpfr_div_ui(c1, c1, 4, MPFR_RNDD);

    mpfr_mul_ui(c2, weinberg, 2, MPFR_RNDD);
    mpfr_sub_ui(c2, c2, 1, MPFR_RNDD);
    mpfr_mul(c2, weinberg, c2, MPFR_RNDD);
    mpfr_div_ui(c2, c2, 2, MPFR_RNDD);
  }

  mpfr_t temp2;
  mpfr_init2(temp2, BITS);

  mpfr_pow_ui(c1_factor, x, 2, MPFR_RNDD);
  mpfr_mul_ui(c1_factor, c1_factor, 14, MPFR_RNDD);
  mpfr_ui_sub(c1_factor, 1, c1_factor, MPFR_RNDD);
  mpfr_pow_ui(temp, x, 4, MPFR_RNDD);
  mpfr_mul_ui(temp, temp, 2, MPFR_RNDD);
  mpfr_sub(c1_factor, c1_factor, temp, MPFR_RNDD);
  mpfr_pow_ui(temp, x, 6, MPFR_RNDD);
  mpfr_mul_ui(temp, temp, 12, MPFR_RNDD);
  mpfr_sub(c1_factor, c1_factor, temp, MPFR_RNDD);
  mpfr_mul(c1_factor, c1_factor, sq1m4x2, MPFR_RNDD);
  mpfr_pow_ui(temp, x, 4, MPFR_RNDD);
  mpfr_mul_ui(temp, temp,12, MPFR_RNDD);

  mpfr_pow_ui(temp2, x, 4, MPFR_RNDD);
  mpfr_sub_ui(temp2, temp2, 1, MPFR_RNDD);
  mpfr_mul(temp, temp, temp2, MPFR_RNDD);
  mpfr_mul(temp, temp, L, MPFR_RNDD);
  mpfr_add(c1_factor, c1_factor, temp, MPFR_RNDD);

  mpfr_pow_ui(c2_factor, x, 4, MPFR_RNDD);
  mpfr_mul_ui(c2_factor, c2_factor, 12, MPFR_RNDD);
  mpfr_pow_ui(temp, x, 2, MPFR_RNDD);
  mpfr_mul_ui(temp, temp, 10, MPFR_RNDD);
  mpfr_sub(c2_factor, temp, c2_factor, MPFR_RNDD);
  mpfr_add_ui(c2_factor, c2_factor, 2, MPFR_RNDD);
  mpfr_pow_ui(temp, x, 2, MPFR_RNDD);
  mpfr_mul(c2_factor, temp, c2_factor, MPFR_RNDD);
  mpfr_mul(c2_factor, c2_factor, sq1m4x2, MPFR_RNDD);
  mpfr_pow_ui(temp2, x, 4, MPFR_RNDD);
  mpfr_mul_ui(temp2, temp2, 2, MPFR_RNDD);
  mpfr_pow_ui(temp, x, 2, MPFR_RNDD);
  mpfr_mul_ui(temp, temp, 2, MPFR_RNDD);
  mpfr_ui_sub(temp, 1, temp, MPFR_RNDD);
  mpfr_add(temp, temp, temp2, MPFR_RNDD);
  mpfr_mul(temp, temp, L, MPFR_RNDD);
  mpfr_pow_ui(temp2, x, 4, MPFR_RNDD);
  mpfr_mul_ui(temp2, temp2, 6, MPFR_RNDD);
  mpfr_mul(temp, temp, temp2, MPFR_RNDD);
  mpfr_add(c2_factor, c2_factor, temp, MPFR_RNDD);

  mpfr_t result;
  mpfr_init2(result, BITS);
  mpfr_mul(temp, c2, c2_factor, MPFR_RNDD);
  mpfr_mul_ui(temp, temp, 4, MPFR_RNDD);
  mpfr_mul(temp2, c1, c1_factor, MPFR_RNDD);
  mpfr_add(result, temp, temp2, MPFR_RNDD);
  mpfr_mul(result, factor, result, MPFR_RNDD);

  Double_t rval = mpfr_get_d(result, MPFR_RNDD);

  mpfr_clears(fermiC, fermiCsq, pi, VUDsq, SOL, HBAR, result, denominator, temp, x, factor, betamass, HNLmass, NZ, angle, sq1m4x2, L, weinberg, c1, c2, c1_factor, c2_factor, temp2, (mpfr_ptr) 0);
  return rval;
}


Double_t pw_lalpha_lbeta_nubeta(std::shared_ptr<Config> cfg, const Lepton &alpha, const Lepton &beta, const HNL &N) {

  if(not N.mixesWith(alpha)) return 0;
  if(alpha == beta) return 0;
  if(N.getMass() < alpha.getMass() + beta.getMass()) return 0; // this means we don't have enough mass in the HNL to produce decay product on-shell

  mpfr_t fermiC, fermiCsq, pi, VUDsq, SOL, HBAR;
  unsigned int BITS = cfg->getBITS();
  mpfr_init2(fermiC, BITS);
  mpfr_init2(fermiCsq, BITS);
  mpfr_init2(pi, BITS);
  mpfr_init2(VUDsq, BITS);
  mpfr_init2(SOL, BITS);
  mpfr_init2(HBAR, BITS);
  cfg->getFermiCsq(fermiCsq);
  cfg->getFermiC(fermiC);
  cfg->getPi(pi);
  cfg->getVUDsq(VUDsq);
  cfg->getSOL(SOL);
  cfg->getHBAR(HBAR);

  mpfr_t NW, angle, HNLmass;
  mpfr_init2(NW, BITS);
  mpfr_init2(angle, BITS);
  mpfr_init2(HNLmass, BITS);
  mpfr_set_d(NW, 1, MPFR_RNDD);
  mpfr_set_d(angle, N.getAngle(), MPFR_RNDD);
  mpfr_set_d(HNLmass, N.getMass(), MPFR_RNDD);

  Double_t xl = alpha.getMass()/N.getMass();
  Double_t xu = beta.getMass()/N.getMass();
  Double_t xd = 0;

  mpfr_t I;
  mpfr_init2(I, BITS);
  mpfr_set_d(I, I_xu_xd_xl(xu, xd, xl), MPFR_RNDD);

  // temporary vars and result value
  mpfr_t result, temp;
  mpfr_init2(result, BITS);
  mpfr_init2(temp, BITS);

  mpfr_pow_ui(temp, pi, 3, MPFR_RNDD);
  mpfr_mul_ui(temp, temp, 192, MPFR_RNDD);
  mpfr_pow_ui(result, HNLmass, 5, MPFR_RNDD);
  mpfr_mul(result, result, fermiCsq, MPFR_RNDD);
  mpfr_mul(result, result, NW, MPFR_RNDD);
  mpfr_div(result, result, temp, MPFR_RNDD);
  mpfr_mul(result, result, angle, MPFR_RNDD);
  mpfr_mul(result, result, I, MPFR_RNDD);

  Double_t rval = mpfr_get_d(result, MPFR_RNDD);
  mpfr_clears(fermiC, fermiCsq, pi, VUDsq, SOL, HBAR, I, NW, angle, HNLmass, result, temp, (mpfr_ptr) 0);

  return rval;
}

Double_t pw_nualpha_nubeta_nubeta(std::shared_ptr<Config> cfg, const Lepton &alpha, const Lepton &beta, const HNL &N) {
  if(not N.mixesWith(alpha)) return 0;

  /*
   * We just need to load configuration for precision bits and
   * nature constants
   */
   mpfr_t fermiC, fermiCsq, pi, VUDsq, SOL, HBAR;
  unsigned int BITS = cfg->getBITS();
  mpfr_init2(fermiC, BITS);
  mpfr_init2(fermiCsq, BITS);
  mpfr_init2(pi, BITS);
  mpfr_init2(VUDsq, BITS);
  mpfr_init2(SOL, BITS);
  mpfr_init2(HBAR, BITS);
  cfg->getFermiCsq(fermiCsq);
  cfg->getFermiC(fermiC);
  cfg->getPi(pi);
  cfg->getVUDsq(VUDsq);
  cfg->getSOL(SOL);
  cfg->getHBAR(HBAR);
  // end of loading configuration

  // initialize and set vars to M_N, U^2
  mpfr_t HNLmass, angle;
  mpfr_init2(HNLmass, BITS);
  mpfr_init2(angle, BITS);
  mpfr_set_d(angle, N.getAngle(), MPFR_RNDD);
  mpfr_set_d(HNLmass, N.getMass(), MPFR_RNDD);

  unsigned int equal = 0;
  if(alpha==beta) {
    equal += 1;
  }
  equal += 1; // increase by one to get the right factor (1+delta_{alpha,beta})

  mpfr_t result, temp;
  mpfr_init2(result, BITS);
  mpfr_init2(temp, BITS);

  mpfr_pow_ui(result, HNLmass, 5, MPFR_RNDD);
  mpfr_mul(result, fermiCsq, result, MPFR_RNDD);
  mpfr_mul_ui(result, result, equal, MPFR_RNDD);
  mpfr_mul(result, result, angle, MPFR_RNDD);
  mpfr_pow_ui(temp, pi, 3, MPFR_RNDD);
  mpfr_mul_ui(temp, temp, 768, MPFR_RNDD);
  mpfr_div(result, result, temp, MPFR_RNDD);

  Double_t rval = mpfr_get_d(result, MPFR_RNDD);

  mpfr_clears(fermiC, fermiCsq, pi, VUDsq, SOL, HBAR, HNLmass, angle, result, temp, (mpfr_ptr) 0);
  return rval;
}

Double_t pw_neutral_pseudoscalar_mesons(std::shared_ptr<Config> cfg, const Lepton &alpha, const Meson &m, const HNL &N) {
  if(not N.mixesWith(alpha)) return 0;
  if(N.getMass() < m.getMass()) return 0; // this means we don't have enough mass in the HNL to produce decay product on-shell

  mpfr_t fermiC, fermiCsq, pi, VUDsq, SOL, HBAR;
  unsigned int BITS = cfg->getBITS();
  mpfr_init2(fermiC, BITS);
  mpfr_init2(fermiCsq, BITS);
  mpfr_init2(pi, BITS);
  mpfr_init2(VUDsq, BITS);
  mpfr_init2(SOL, BITS);
  mpfr_init2(HBAR, BITS);
  cfg->getFermiCsq(fermiCsq);
  cfg->getFermiC(fermiC);
  cfg->getPi(pi);
  cfg->getVUDsq(VUDsq);
  cfg->getSOL(SOL);
  cfg->getHBAR(HBAR);

  // initialize high precision variables
  mpfr_t xhsq, fhsq;
  mpfr_t mesonMass, HNLmass, angle;
  mpfr_init2(xhsq, BITS);
  mpfr_init2(fhsq, BITS);
  mpfr_init2(mesonMass, BITS);
  mpfr_init2(HNLmass, BITS);
  mpfr_init2(angle, BITS);

  // set high precision variables to values
  mpfr_set_d(mesonMass, m.getMass(), MPFR_RNDD);
  mpfr_set_d(HNLmass, N.getMass(), MPFR_RNDD);
  mpfr_set_d(angle, N.getAngle(), MPFR_RNDD);
  mpfr_set_d(fhsq, m.getDecayConstant(), MPFR_RNDD);
  mpfr_pow_ui(fhsq, fhsq, 2, MPFR_RNDD);
  mpfr_div(xhsq, mesonMass, HNLmass, MPFR_RNDD);
  mpfr_pow_ui(xhsq, xhsq, 2, MPFR_RNDD);

  // create result and temp variables
  mpfr_t result, temp;
  mpfr_init2(temp, BITS);
  mpfr_init2(result, BITS);

  mpfr_mul(result, fermiCsq, fhsq, MPFR_RNDD);
  mpfr_mul(result, result, angle, MPFR_RNDD);
  mpfr_pow_ui(temp, HNLmass, 3, MPFR_RNDD);
  mpfr_mul(result, result, temp, MPFR_RNDD);
  mpfr_mul_ui(temp, pi, 32, MPFR_RNDD);
  mpfr_div(result, result, temp, MPFR_RNDD);

  mpfr_ui_sub(temp, 1, xhsq, MPFR_RNDD);
  mpfr_pow_ui(temp, temp, 2, MPFR_RNDD);
  mpfr_mul(result, result, temp, MPFR_RNDD);

  Double_t rval = mpfr_get_d(result, MPFR_RNDD);
  mpfr_clears(fermiC, fermiCsq, pi, VUDsq, SOL, HBAR, xhsq, fhsq, mesonMass, HNLmass, angle, result, temp, (mpfr_ptr) 0);

  return rval;
}

Double_t pw_charged_pseudoscalar_mesons(std::shared_ptr<Config> cfg, const Lepton &alpha, const Meson &m, const HNL &N) {
  if(not N.mixesWith(alpha)) return 0;
  if(N.getMass() < alpha.getMass() + m.getMass()) return 0; // this means we don't have enough mass in the HNL to produce decay product on-shell

  mpfr_t fermiC, fermiCsq, pi, VUDsq, SOL, HBAR;
  unsigned int BITS = cfg->getBITS();
  mpfr_init2(fermiC, BITS);
  mpfr_init2(fermiCsq, BITS);
  mpfr_init2(pi, BITS);
  mpfr_init2(VUDsq, BITS);
  mpfr_init2(SOL, BITS);
  mpfr_init2(HBAR, BITS);
  cfg->getFermiCsq(fermiCsq);
  cfg->getFermiC(fermiC);
  cfg->getPi(pi);
  cfg->getVUDsq(VUDsq);
  cfg->getSOL(SOL);
  cfg->getHBAR(HBAR);

  // initialize high precision variables
  mpfr_t xhsq, xlsq, fh;
  mpfr_t mesonMass, alphaMass, HNLmass, angle;
  mpfr_init2(xhsq, BITS);
  mpfr_init2(xlsq, BITS);
  mpfr_init2(fh, BITS);
  mpfr_init2(mesonMass, BITS);
  mpfr_init2(alphaMass, BITS);
  mpfr_init2(HNLmass, BITS);
  mpfr_init2(angle, BITS);

  // set high precision variables to values
  mpfr_set_d(mesonMass, m.getMass(), MPFR_RNDD);
  mpfr_set_d(alphaMass, alpha.getMass(), MPFR_RNDD);
  mpfr_set_d(HNLmass, N.getMass(), MPFR_RNDD);
  mpfr_set_d(angle, N.getAngle(), MPFR_RNDD);
  mpfr_set_d(fh, m.getDecayConstant(), MPFR_RNDD);
  mpfr_div(xlsq, alphaMass, HNLmass, MPFR_RNDD);
  mpfr_pow_ui(xlsq, xlsq, 2, MPFR_RNDD);
  mpfr_div(xhsq, mesonMass, HNLmass, MPFR_RNDD);
  mpfr_pow_ui(xhsq, xhsq, 2, MPFR_RNDD);

  // create result and temp variables
  mpfr_t result, temp, temp2;
  mpfr_init2(temp, BITS);
  mpfr_init2(temp2, BITS);
  mpfr_init2(result, BITS);

  mpfr_pow_ui(temp, fh, 2, MPFR_RNDD);
  mpfr_mul(result, fermiCsq, temp, MPFR_RNDD);
  mpfr_mul(result, result, VUDsq, MPFR_RNDD);
  mpfr_mul(result, result, angle, MPFR_RNDD);
  mpfr_pow_ui(temp, HNLmass, 3, MPFR_RNDD);
  mpfr_mul(result, result, temp, MPFR_RNDD);
  mpfr_mul_ui(temp, pi, 16, MPFR_RNDD);
  mpfr_div(result, result, temp, MPFR_RNDD);

  mpfr_ui_sub(temp, 1, xlsq, MPFR_RNDD);
  mpfr_pow_ui(temp, temp, 2, MPFR_RNDD);

  mpfr_add_ui(temp2, xlsq, 1, MPFR_RNDD);
  mpfr_mul(temp2, temp2, xhsq, MPFR_RNDD);
  mpfr_sub(temp, temp, temp2, MPFR_RNDD);
  mpfr_mul(result, result, temp, MPFR_RNDD);

  mpfr_t one;
  mpfr_init2(one, BITS);
  mpfr_set_d(one, 1, MPFR_RNDD);

  kaellen(cfg, temp, one, xhsq, xlsq);
  mpfr_sqrt(temp, temp, MPFR_RNDD);

  mpfr_mul(result, result, temp, MPFR_RNDD);

  Double_t rval = mpfr_get_d(result, MPFR_RNDD);
  mpfr_clears(fermiC, fermiCsq, pi, VUDsq, SOL, HBAR, xhsq, xlsq, fh, mesonMass, alphaMass, HNLmass, angle, result, temp, temp2, one, (mpfr_ptr) 0);

  return rval;
}

Double_t pw_charged_vector_mesons(std::shared_ptr<Config> cfg, const Lepton &alpha, const Meson &m, const HNL &N) {
  if(not N.mixesWith(alpha)) return 0;
  if(N.getMass() < alpha.getMass() + m.getMass()) return 0; // this means we don't have enough mass in the HNL to produce decay product on-shell

  mpfr_t fermiC, fermiCsq, pi, VUDsq, SOL, HBAR;
  unsigned int BITS = cfg->getBITS();
  mpfr_init2(fermiC, BITS);
  mpfr_init2(fermiCsq, BITS);
  mpfr_init2(pi, BITS);
  mpfr_init2(VUDsq, BITS);
  mpfr_init2(SOL, BITS);
  mpfr_init2(HBAR, BITS);
  cfg->getFermiCsq(fermiCsq);
  cfg->getFermiC(fermiC);
  cfg->getPi(pi);
  cfg->getVUDsq(VUDsq);
  cfg->getSOL(SOL);
  cfg->getHBAR(HBAR);

  // initialize high precision variables
  mpfr_t xhsq, xlsq, gh;
  mpfr_t mesonMass, alphaMass, HNLmass, angle;
  mpfr_init2(xhsq, BITS);
  mpfr_init2(xlsq, BITS);
  mpfr_init2(gh, BITS);
  mpfr_init2(mesonMass, BITS);
  mpfr_init2(alphaMass, BITS);
  mpfr_init2(HNLmass, BITS);
  mpfr_init2(angle, BITS);

  // set high precision variables to values
  mpfr_set_d(mesonMass, m.getMass(), MPFR_RNDD);
  mpfr_set_d(alphaMass, alpha.getMass(), MPFR_RNDD);
  mpfr_set_d(HNLmass, N.getMass(), MPFR_RNDD);
  mpfr_set_d(angle, N.getAngle(), MPFR_RNDD);
  mpfr_set_d(gh, m.getDecayConstant(), MPFR_RNDD);
  mpfr_div(xlsq, alphaMass, HNLmass, MPFR_RNDD);
  mpfr_pow_ui(xlsq, xlsq, 2, MPFR_RNDD);
  mpfr_div(xhsq, mesonMass, HNLmass, MPFR_RNDD);
  mpfr_pow_ui(xhsq, xhsq, 2, MPFR_RNDD);

  // create result and temp variables
  mpfr_t result, temp, temp2;
  mpfr_init2(temp, BITS);
  mpfr_init2(temp2, BITS);
  mpfr_init2(result, BITS);

  mpfr_pow_ui(temp, gh, 2, MPFR_RNDD);
  mpfr_mul(result, fermiCsq, temp, MPFR_RNDD);
  mpfr_mul(result, result, VUDsq, MPFR_RNDD);
  mpfr_mul(result, result, angle, MPFR_RNDD);
  mpfr_pow_ui(temp, HNLmass, 3, MPFR_RNDD);
  mpfr_mul(result, result, temp, MPFR_RNDD);
  mpfr_mul_ui(temp, pi, 16, MPFR_RNDD);
  mpfr_pow_ui(temp2, mesonMass, 2, MPFR_RNDD);
  mpfr_mul(temp, temp, temp2, MPFR_RNDD);
  mpfr_div(result, result, temp, MPFR_RNDD);

  mpfr_ui_sub(temp, 1, xlsq, MPFR_RNDD);
  mpfr_pow_ui(temp, temp, 2, MPFR_RNDD);

  mpfr_add_ui(temp2, xlsq, 1, MPFR_RNDD);
  mpfr_mul(temp2, temp2, xhsq, MPFR_RNDD);
  mpfr_add(temp, temp, temp2, MPFR_RNDD);
  mpfr_pow_ui(temp2, xhsq, 2, MPFR_RNDD);
  mpfr_mul_ui(temp2, temp2, 2, MPFR_RNDD);
  mpfr_sub(temp2, temp, temp2, MPFR_RNDD);

  mpfr_t one;
  mpfr_init2(one, BITS);
  mpfr_set_d(one, 1, MPFR_RNDD);

  kaellen(cfg, temp, one, xhsq, xlsq);
  mpfr_sqrt(temp, temp, MPFR_RNDD);
  mpfr_mul(temp, temp, temp2, MPFR_RNDD);
  LOG_DEBUG("value of after factor: " << mpfr_get_d(temp, MPFR_RNDD));
  LOG_DEBUG("value of Gamma0: " << mpfr_get_d(result, MPFR_RNDD));

  mpfr_mul(result, result, temp, MPFR_RNDD);

  Double_t rval = mpfr_get_d(result, MPFR_RNDD);
  mpfr_clears(fermiC, fermiCsq, pi, VUDsq, SOL, HBAR, one, result, temp, temp2, xhsq, xlsq, gh, mesonMass, alphaMass, HNLmass, angle, (mpfr_ptr) 0);

  return rval;
}

Double_t pw_neutral_vector_mesons(std::shared_ptr<Config> cfg, const Lepton &alpha, const Meson &m, const HNL &N) {
  if(not N.mixesWith(alpha)) return 0;
  if(not m.hasValue("kh"))
    throw std::runtime_error("Neutral mesons need to have dimensionless kh factor stored as 'kh' value.");
  if(N.getMass() < m.getMass()) return 0; // this means we don't have enough mass in the HNL to produce decay product on-shell

  mpfr_t fermiC, fermiCsq, pi, VUDsq, SOL, HBAR;
  unsigned int BITS = cfg->getBITS();
  mpfr_init2(fermiC, BITS);
  mpfr_init2(fermiCsq, BITS);
  mpfr_init2(pi, BITS);
  mpfr_init2(VUDsq, BITS);
  mpfr_init2(SOL, BITS);
  mpfr_init2(HBAR, BITS);
  cfg->getFermiCsq(fermiCsq);
  cfg->getFermiC(fermiC);
  cfg->getPi(pi);
  cfg->getVUDsq(VUDsq);
  cfg->getSOL(SOL);
  cfg->getHBAR(HBAR);

  // initialize high precision variables
  mpfr_t xh, gh, kh;
  mpfr_t mesonMass, HNLmass, angle;
  mpfr_init2(xh, BITS);
  mpfr_init2(gh, BITS);
  mpfr_init2(kh, BITS);
  mpfr_init2(mesonMass, BITS);
  mpfr_init2(HNLmass, BITS);
  mpfr_init2(angle, BITS);

  // set high precision variables to values
  mpfr_set_d(mesonMass, m.getMass(), MPFR_RNDD);
  mpfr_set_d(HNLmass, N.getMass(), MPFR_RNDD);
  mpfr_set_d(angle, N.getAngle(), MPFR_RNDD);
  mpfr_set_d(gh, m.getDecayConstant(), MPFR_RNDD);
  mpfr_set_d(kh, m.getValue("kh"), MPFR_RNDD);
  mpfr_div(xh, mesonMass, HNLmass, MPFR_RNDD);

  // create result and temp variables
  mpfr_t result, temp, temp2;
  mpfr_init2(temp, BITS);
  mpfr_init2(temp2, BITS);
  mpfr_init2(result, BITS);

  mpfr_pow_ui(temp, gh, 2, MPFR_RNDD);
  mpfr_pow_ui(temp2, kh, 2, MPFR_RNDD);
  mpfr_mul(result, fermiCsq, temp, MPFR_RNDD);
  mpfr_mul(result, result, temp2, MPFR_RNDD);
  mpfr_mul(result, result, VUDsq, MPFR_RNDD);
  mpfr_mul(result, result, angle, MPFR_RNDD);
  mpfr_pow_ui(temp, HNLmass, 3, MPFR_RNDD);
  mpfr_mul(result, result, temp, MPFR_RNDD);
  mpfr_mul_ui(temp, pi, 32, MPFR_RNDD);
  mpfr_pow_ui(temp2, mesonMass, 2, MPFR_RNDD);
  mpfr_mul(temp, temp, temp2, MPFR_RNDD);
  mpfr_div(result, result, temp, MPFR_RNDD);

  mpfr_pow_ui(temp, xh, 2, MPFR_RNDD);
  mpfr_mul_ui(temp, temp, 2, MPFR_RNDD);
  mpfr_add_ui(temp, temp, 1, MPFR_RNDD);

  mpfr_pow_ui(temp2, xh, 2, MPFR_RNDD);
  mpfr_ui_sub(temp2, 1, temp2, MPFR_RNDD);
  mpfr_pow_ui(temp2, temp2, 2, MPFR_RNDD);
  mpfr_mul(result, result, temp, MPFR_RNDD);
  mpfr_mul(result, result, temp2, MPFR_RNDD);

  Double_t rval = mpfr_get_d(result, MPFR_RNDD);
  mpfr_clears(fermiC, fermiCsq, pi, VUDsq, SOL, HBAR, result, temp, temp2, xh, gh, kh,mesonMass, HNLmass, angle, (mpfr_ptr) 0);

  return rval;
}
