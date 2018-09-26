//  cHNLdecay  --  calculate decay widths of Heavy Neutral Leptons
//  Copyright (C) 2018 - Fabian A.J. Thiele, <fabian.thiele@posteo.de>
//
//  This file is part of cHNLdecay.
//
//  cHNLdecay is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  cHNLdecay is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <https://www.gnu.org/licenses/>.

#ifndef CONFIG_H
#define CONFIG_H

#include "TMath.h"
#include "TString.h"
#include "Meson.h"
#include <gmp.h>
#include <mpfr.h>

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
    mpfr_clears(fermiC, fermiCsq, pi, SOL, HBAR, (mpfr_ptr)0);
  }

  TString getName() const { return name; }

  unsigned int getBITS() { return BITS; }
  void getFermiC(mpfr_t result) { mpfr_set(result, fermiC, MPFR_RNDD); }
  void getFermiCsq(mpfr_t result) { mpfr_set(result, fermiCsq, MPFR_RNDD); }
  void getPi(mpfr_t result) { mpfr_set(result, pi, MPFR_RNDD); }
  void getVUDsq(mpfr_t result, const Meson &m) {
    mpfr_set_d(result, ckm.at({m.getU(), m.getD()}), MPFR_RNDD);
    mpfr_pow_ui(result, result, 2, MPFR_RNDD);
  }
  void getSOL(mpfr_t result) { mpfr_set(result, SOL, MPFR_RNDD); }
  void getHBAR(mpfr_t result) { mpfr_set(result, HBAR, MPFR_RNDD); }

private:
  TString name;
  unsigned int BITS;
  mpfr_t fermiC, fermiCsq, pi, SOL, HBAR;

    std::map<std::pair<Quark, Quark>, Double_t> ckm = {
                                                          {{Quark::up, Quark::down}, 0.97427},
                                                          {{Quark::up, Quark::strange}, 0.22534},
                                                          {{Quark::up, Quark::bottom}, 0.00351},
                                                          {{Quark::charm, Quark::down}, 0.22520},
                                                          {{Quark::charm, Quark::strange}, 0.97344},
                                                          {{Quark::charm, Quark::bottom}, 0.0412},
                                                          {{Quark::top, Quark::down}, 0.00867},
                                                          {{Quark::top, Quark::strange}, 0.0404},
                                                          {{Quark::top, Quark::bottom}, 0.999146}
    };


  void initialize() {
    mpfr_init2(fermiC, BITS);
    mpfr_set_d(fermiC, 1.1663787e-11,
               MPFR_RNDD); //  GF/(hbar c)^3 fermi constant in MeV^-2

    mpfr_init2(fermiCsq, BITS);
    mpfr_pow_ui(fermiCsq, fermiC, 2,
                MPFR_RNDD); // save square of fermiC for easier use

    mpfr_init2(pi, BITS);
    mpfr_set_d(pi, TMath::Pi(), MPFR_RNDD);

    mpfr_init2(SOL, BITS);
    mpfr_set_d(SOL, 299792458, MPFR_RNDD);

    mpfr_init2(HBAR, BITS);
    mpfr_set_d(HBAR, 6.582119514e-16, MPFR_RNDD);
  }
};

#endif
