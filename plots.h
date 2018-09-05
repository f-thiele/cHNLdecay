#ifndef   PLOTS_H
#define   PLOTS_H
#include "Lepton.h"
#include "Meson.h"
#include "HNL.h"
#include <iostream>
#include "TString.h"

void plot_meson_pw(std::shared_ptr<Config> cfg, Lepton alpha, std::vector<Meson> mesons, HNL N, TString output, Int_t lowMass, Int_t highMass, Int_t stepsize=10);
#endif
