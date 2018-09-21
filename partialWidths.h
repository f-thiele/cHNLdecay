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

#ifndef PARTIALWIDTHS_H
#define PARTIALWIDTHS_H
#include "HNL.h"
#include "Lepton.h"
#include "Meson.h"
#include <vector>

Double_t pw_nualpha_lbeta_lbeta(std::shared_ptr<Config> cfg,
                                const Lepton &alpha, const Lepton &beta,
                                const HNL &N);
Double_t pw_lalpha_lbeta_nubeta(std::shared_ptr<Config> cfg,
                                const Lepton &alpha, const Lepton &beta,
                                const HNL &N);
Double_t pw_nualpha_nubeta_nubeta(std::shared_ptr<Config> cfg,
                                  const Lepton &alpha, const Lepton &beta,
                                  const HNL &N);
Double_t pw_neutral_pseudoscalar_mesons(std::shared_ptr<Config> cfg,
                                        const Lepton &alpha, const Meson &m,
                                        const HNL &N);
Double_t pw_charged_pseudoscalar_mesons(std::shared_ptr<Config> cfg,
                                        const Lepton &alpha, const Meson &m,
                                        const HNL &N);
Double_t pw_charged_vector_mesons(std::shared_ptr<Config> cfg,
                                  const Lepton &alpha, const Meson &m,
                                  const HNL &N);
Double_t pw_neutral_vector_mesons(std::shared_ptr<Config> cfg,
                                  const Lepton &alpha, const Meson &m,
                                  const HNL &N);
#endif
