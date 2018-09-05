#ifndef   PARTIALWIDTHS_H
#define   PARTIALWIDTHS_H
#include <vector>
#include "Meson.h"
#include "Lepton.h"
#include "HNL.h"

Double_t pw_nualpha_lbeta_lbeta(std::shared_ptr<Config> cfg, Lepton alpha, Lepton beta, HNL N);
Double_t pw_lalpha_lbeta_nubeta(std::shared_ptr<Config> cfg, Lepton alpha, Lepton beta, HNL N);
Double_t pw_nualpha_nubeta_nubeta(std::shared_ptr<Config> cfg, Lepton alpha, Lepton beta, HNL N);
Double_t pw_neutral_pseudoscalar_mesons(std::shared_ptr<Config> cfg, Lepton alpha, Meson m, HNL N);
Double_t pw_charged_pseudoscalar_mesons(std::shared_ptr<Config> cfg, Lepton alpha, Meson m, HNL N);
Double_t pw_charged_vector_mesons(std::shared_ptr<Config> cfg, Lepton alpha, Meson m, HNL N);
Double_t pw_neutral_vector_mesons(std::shared_ptr<Config> cfg, Lepton alpha, Meson m, HNL N);
#endif
