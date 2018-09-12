// Copyright (C) 2018 - Fabian A.J. Thiele, <fabian.thiele@cern.ch>

#include "HNL.h"
#include "partialWidths.h"
#include "auxfunctions.h"
#include "Logger.h"
#include "TF1.h"

Double_t HNL::getPartialWidth(std::shared_ptr<Config> cfg, const Lepton &alpha, const Lepton &beta, bool invisible) {
  Double_t nununu = 0;
  Double_t nua_lb_lb = std::max(0., pw_nualpha_lbeta_lbeta(cfg, alpha, beta, *this));
  Double_t la_lb_nub = std::max(0., pw_lalpha_lbeta_nubeta(cfg, alpha, beta, *this));
  if(invisible) {
    nununu += std::max(0., pw_nualpha_nubeta_nubeta(cfg, alpha, beta, *this));
  }

  this->newDecayChannel("\\nu_" + alpha.getName() + ", " + beta.getName() + ", " + beta.getName(), nua_lb_lb);
  this->newDecayChannel(alpha.getName() + "^-, " + beta.getName() + ", \\nu_" + beta.getName(), la_lb_nub);
  this->newDecayChannel("\\nu_" + alpha.getName() + ", \\nu_" + beta.getName() + ", \\nu_" + beta.getName(), nununu);

  Double_t pw = nununu + nua_lb_lb + la_lb_nub;

  return pw;
}


Double_t HNL::getPartialWidthInv(std::shared_ptr<Config> cfg, const Lepton &alpha, const Lepton &beta) {
  Double_t nununu = std::max(0., pw_nualpha_nubeta_nubeta(cfg, alpha, beta, *this));

  TString channel_name = "\\nu_" + alpha.getName() + ", \\nu_" + beta.getName() + ", \\nu_" + beta.getName();
  if(not existDecayChannel(name)) {
    this->newDecayChannel(name, nununu);
  }

  return nununu;
}


Double_t HNL::getPartialWidth(std::shared_ptr<Config> cfg, const Lepton &alpha, const Meson &m) {
  Double_t dw = 0;
  Double_t temp = 0;

  if(m.getMesonType() == MesonType::pseudoscalar and m.getCharge() == Charge::charged) {
    temp = std::max(0., pw_charged_pseudoscalar_mesons(cfg, alpha, m, *this));
    this->newDecayChannel(alpha.getName() + "^-, " + m.getName(), temp);

    dw += temp;

  } else if (m.getMesonType() == MesonType::pseudoscalar and m.getCharge() == Charge::neutral) {
    temp = std::max(0., pw_neutral_pseudoscalar_mesons(cfg, alpha, m, *this));
    this->newDecayChannel("\\nu_" + alpha.getName() + ", " + m.getName(), temp);

    dw += temp;

  } else if (m.getMesonType() == MesonType::vector and m.getCharge() == Charge::charged) {
    temp = std::max(0., pw_charged_vector_mesons(cfg, alpha, m, *this));
    this->newDecayChannel(alpha.getName() + "^-, " + m.getName(), temp);

    dw += temp;

  } else if (m.getMesonType() == MesonType::vector and m.getCharge() == Charge::neutral) {
    temp = std::max(0., pw_neutral_vector_mesons(cfg, alpha, m, *this));
    this->newDecayChannel("\\nu_" + alpha.getName() + ", " + m.getName(), temp);

    dw += temp;

  } else if (m.getMesonType() == MesonType::Unknown) {
    throw std::runtime_error("Cannot use meson with undefined MesonType");
  }

  return dw;
}

Double_t HNL::getTotalWidth(std::shared_ptr<Config> cfg, const std::vector<Lepton> &leptons, const std::vector<Meson> &mesons) {
  Double_t tw_lept = 0;
  Double_t tw_mes = 0;
  this->clearDecayChannels();

  for(auto l1 : leptons) {
    for(auto l2 : leptons) {
      Double_t t = this->getPartialWidth(cfg, l1, l2);
      if(t>0) LOG_DEBUG(l1.getName() << " " << l2.getName() << ": " << t);
      tw_lept += t;
    }
    for(auto m : mesons) {
      Double_t t = this->getPartialWidth(cfg, l1, m);
      if(t>0) LOG_DEBUG(l1.getName() << " " << m.getName() << ": " << t);
      tw_mes += t;
    }
  }
  TF1* f = new TF1("#Delta_{QCD}", qcd_coupling, 1, 100, 0);
  Double_t qcd_corr = qcd_correction(f->Eval(this->getMass()/1000.));
  LOG_DEBUG("QCD correction: " << qcd_corr);
  delete f;

  Double_t totalWidth = (1+qcd_corr)*tw_mes + tw_lept;
  totalWidth *= 1.; // multiply by two for majorana channels

  return totalWidth;
}
