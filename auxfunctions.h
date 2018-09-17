// Copyright (C) 2018 - Fabian A.J. Thiele, <fabian.thiele@cern.ch>

#ifndef   AUXFUNCTIONS_H
#define   AUXFUNCTIONS_H
#include "TGraph.h"
#include <gmp.h>
#include <mpfr.h>
#include "Config.h"
#include "Lepton.h"
#include "Meson.h"
#include "HNL.h"

void kaellen(std::shared_ptr<Config> cfg, mpfr_t result, mpfr_t a, mpfr_t b, mpfr_t c);
Double_t kaellen(Double_t a, Double_t b, Double_t c);
Double_t int__i(Double_t* x, Double_t* par);
Double_t I_xu_xd_xl(Double_t xu, Double_t xd, Double_t xl);
Double_t gamma2ctau(std::shared_ptr<Config> cfg, Double_t gamma);
Double_t qcd_coupling(Double_t* x, Double_t* par);
Double_t qcd_coupling(Double_t mass);
Double_t qcd_correction(Double_t alpha);
Double_t f_qcd_correction(Double_t* x, Double_t* par);
Double_t ctauToU2(std::shared_ptr<Config> cfg, Double_t target, const std::vector<Lepton> &leptons, const std::vector<Meson> &mesons, HNL &N, Double_t start=1e-3, Double_t tol=1e-6, Double_t stepsize=1-1e-4);

TGraph* create_graph(Double_t xd, Double_t xl, Float_t low=0, Float_t high=0.5, Float_t stepsize=0.001);
void plot_I();
#endif
