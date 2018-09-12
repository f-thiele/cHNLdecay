// Copyright (C) 2018 - Fabian A.J. Thiele, <fabian.thiele@cern.ch>

#ifndef   CONFIG_H
#define   CONFIG_H

#include <gmp.h>
#include <mpfr.h>
#include "TMath.h"

class Config {
 public:
  Config() {
    name = "default";
    BITS = 400;
    initialize();
  }
  Config(TString n, unsigned int b) {
    name = n;
    BITS = b;
    initialize();
  }

  ~Config(void) {
    mpfr_clears(fermiC, fermiCsq, pi, VUDsq, SOL, HBAR, (mpfr_ptr) 0);
  }

  TString getName() const {
    return name;
  }

  unsigned int getBITS() {
    return BITS;
  }
  void getFermiC(mpfr_t result) {
    mpfr_set(result, fermiC, MPFR_RNDD);
  }
  void getFermiCsq(mpfr_t result) {
    mpfr_set(result, fermiCsq, MPFR_RNDD);
  }
  void getPi(mpfr_t result) {
    mpfr_set(result, pi, MPFR_RNDD);
  }
  void getVUDsq(mpfr_t result) {
    mpfr_set(result, VUDsq, MPFR_RNDD);
  }
  void getSOL(mpfr_t result) {
    mpfr_set(result, SOL, MPFR_RNDD);
  }
  void getHBAR(mpfr_t result) {
    mpfr_set(result, HBAR, MPFR_RNDD);
  }

 private:
  TString name;
  unsigned int BITS;
  mpfr_t fermiC, fermiCsq, pi, VUDsq, SOL, HBAR;

  void initialize () {
    mpfr_init2(fermiC, BITS);
    mpfr_set_d(fermiC, 1.1663787e-11, MPFR_RNDD); //  GF/(hbar c)^3 fermi constant in MeV^-2

    mpfr_init2(fermiCsq, BITS);
    mpfr_pow_ui(fermiCsq, fermiC, 2, MPFR_RNDD); // save square of fermiC for easier use

    mpfr_init2(pi, BITS);
    mpfr_set_d(pi, TMath::Pi(), MPFR_RNDD);

    mpfr_init2(VUDsq, BITS);
    mpfr_set_d(VUDsq, 0.97425, MPFR_RNDD);
    mpfr_pow_ui(VUDsq, VUDsq, 2, MPFR_RNDD);

    mpfr_init2(SOL, BITS);
    mpfr_set_d(SOL, 299792458, MPFR_RNDD);

    mpfr_init2(HBAR, BITS);
    mpfr_set_d(HBAR, 6.582119514e-16, MPFR_RNDD);
  }
};

#endif
