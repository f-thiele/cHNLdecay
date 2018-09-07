// Copyright (C) 2018 - Fabian A.J. Thiele, <fabian.thiele@cern.ch>

#ifndef   PARTIALWIDTHS_H
#define   PARTIALWIDTHS_H
#include <vector>
#include "Meson.h"
#include "Lepton.h"
#include "HNL.h"

Double_t pw_nualpha_lbeta_lbeta(std::shared_ptr<Config> cfg, const Lepton &alpha, const Lepton &beta, const HNL &N);
Double_t pw_lalpha_lbeta_nubeta(std::shared_ptr<Config> cfg, const Lepton &alpha, const Lepton &beta, const HNL &N);
Double_t pw_nualpha_nubeta_nubeta(std::shared_ptr<Config> cfg, const Lepton &alpha, const Lepton &beta, const HNL &N);
Double_t pw_neutral_pseudoscalar_mesons(std::shared_ptr<Config> cfg, const Lepton &alpha, const Meson &m, const HNL &N);
Double_t pw_charged_pseudoscalar_mesons(std::shared_ptr<Config> cfg, const Lepton &alpha, const Meson &m, const HNL &N);
Double_t pw_charged_vector_mesons(std::shared_ptr<Config> cfg, const Lepton &alpha, const Meson &m, const HNL &N);
Double_t pw_neutral_vector_mesons(std::shared_ptr<Config> cfg, const Lepton &alpha, const Meson &m, const HNL &N);
#endif
