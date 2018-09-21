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

#ifndef PLOTS_H
#define PLOTS_H
#include "HNL.h"
#include "Lepton.h"
#include "Meson.h"
#include "TString.h"
#include <iostream>

void plot_meson_pw(std::shared_ptr<Config> cfg, Lepton alpha,
                   std::vector<Meson> mesons, HNL N, TString output,
                   Int_t lowMass, Int_t highMass, Int_t stepsize = 10);
void plot_br(std::shared_ptr<Config> cfg, std::vector<Lepton> leptons,
             std::vector<Meson> mesons, HNL N, TString output, Int_t lowMass,
             Int_t highMass, Int_t stepsize = 10);
void plot_qcd_coupling(TString output);
void plot_qcd_correction(TString output);
void plot_br_low(std::shared_ptr<Config> cfg, std::vector<Lepton> leptons,
                 std::vector<Meson> mesons, HNL N, TString output,
                 Int_t lowMass, Int_t highMass, Int_t stepsize = 10);
#endif
