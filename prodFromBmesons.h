#ifndef PRODFROMBMESONS_H
#define PRODFROMBMESONS_H
#include "HNL.h"
#include "Lepton.h"
#include "Meson.h"
#include "auxfunctions.h"
#include <vector>

Double_t kallen(Double_t a, Double_t b, Double_t c);
Double_t Gamma(Double_t eta, Double_t yhp, Double_t yl, Double_t yN);
Double_t Gminus(Double_t eta, Double_t yhp, Double_t yl, Double_t yN);
Double_t Gplus (Double_t xi, Double_t yhp, Double_t yl, Double_t yN);
Double_t F(Double_t xi, Double_t yhp) ;


Double_t compute_fP1(std::shared_ptr<Config> cfg, HNL N, Lepton l, Meson meson, Meson mesonp, Double_t eta_);
Double_t compute_fP2(std::shared_ptr<Config> cfg, HNL N, Lepton l, Meson meson, Meson mesonp, Double_t eta_);
Double_t compute_fP3(std::shared_ptr<Config> cfg, HNL N, Lepton l, Meson meson, Meson mesonp, Double_t eta_);

Double_t integral_fP(Double_t a, Double_t b, Double_t nsteps, std::shared_ptr<Config> cfg, HNL N, Lepton l, Meson meson, Meson mesonp, int fct);

Double_t compute_fV(std::shared_ptr<Config> cfg, HNL N, Lepton l, Meson meson, Meson mesonp, Double_t xi_, int part);
Double_t compute_integral(Double_t a, Double_t b, Double_t nsteps, std::shared_ptr<Config>, HNL, Lepton, Meson, Meson);
Double_t MONITORING(Double_t nsteps, std::shared_ptr<Config>, Lepton, Meson, Meson);

Double_t pw_prodFromBmeson_semileptonic(std::shared_ptr<Config> cfg, HNL N, Lepton l, Meson meson, Meson mesonp);
Double_t pw_prodFromBmeson_leptonic(std::shared_ptr<Config> cfg, HNL N, Lepton l, Meson meson);

Double_t compute_dint_Lambdab_proton(std::shared_ptr<Config> cfg, HNL N, Lepton l, Double_t mh, Double_t mhp, Double_t q2);
Double_t compute_dint_Lambdab_Lambdac(std::shared_ptr<Config> cfg, HNL N, Lepton l, Double_t mh, Double_t mhp, Double_t q2);

Double_t BR_prod_Lambdab_proton(std::shared_ptr<Config> cfg, HNL N, Lepton l);
Double_t BR_prod_Lambdab_Lambdac(std::shared_ptr<Config> cfg, HNL N, Lepton l);

void display_fform(HNL N, Lepton l, Meson meson, Meson mesonp);
#endif
