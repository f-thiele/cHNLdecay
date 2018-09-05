#include "HNL.h"
#include "partialWidths.h"

Double_t HNL::getPartialWidth(std::shared_ptr<Config> cfg, Lepton &alpha, Lepton &beta, bool invisible) const {
  Double_t dw = 0;
  dw += std::max(0., pw_nualpha_lbeta_lbeta(cfg, alpha, beta, *this));
  dw += std::max(0., pw_lalpha_lbeta_nubeta(cfg, alpha, beta, *this));
  if(invisible)
    dw += std::max(0., pw_nualpha_nubeta_nubeta(cfg, alpha, beta, *this));

  return dw;
}


Double_t HNL::getPartialWidthInv(std::shared_ptr<Config> cfg, Lepton &alpha, Lepton &beta) const {
  return std::max(0., pw_nualpha_nubeta_nubeta(cfg, alpha, beta, *this));
}


Double_t HNL::getPartialWidth(std::shared_ptr<Config> cfg, Lepton &alpha, Meson &m) const {
  Double_t dw = 0;

  if(m.getMesonType() == MesonType::pseudoscalar) {
    if (m.getCharge() == Charge::charged) {
      dw += std::max(0., pw_charged_pseudoscalar_mesons(cfg, alpha, m, *this));
    } else if (m.getCharge() == Charge::neutral) {
      dw += std::max(0., pw_neutral_pseudoscalar_mesons(cfg, alpha, m, *this));
    }
  } else if (m.getMesonType() == MesonType::vector) {
    if (m.getCharge() == Charge::charged) {
      dw += std::max(0., pw_charged_vector_mesons(cfg, alpha, m, *this));
    } else if (m.getCharge() == Charge::neutral) {
      dw += std::max(0., pw_neutral_vector_mesons(cfg, alpha, m, *this));
    }
  } else if (m.getMesonType() == MesonType::Unknown) {
    throw std::runtime_error("Cannot use meson with undefined MesonType");
  }
  return dw;
}
