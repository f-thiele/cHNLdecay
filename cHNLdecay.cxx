// Copyright (C) 2018 - Fabian A.J. Thiele, <fabian.thiele@cern.ch>

#include "TGraph.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TMath.h"
#include "TF1.h"

#include "HNL.h"
#include "Lepton.h"
#include "Meson.h"

#include <iostream>
#include <gmp.h>
#include <mpfr.h>
#include "Config.h"
#include "auxfunctions.h"
#include "partialWidths.h"
#include "plots.h"
#include "Logger.h"

Level gLOGLEVEL;

int main() {
  gLOGLEVEL = Level::INFO;

  auto cfg = std::make_shared<Config>(); // set BITS and initialize constants

  Double_t muonMass = 105.6583715; // MeV
  Double_t electronMass = 0.5109989461; // MeV
  Double_t tauMass = 1776.82; // MeV

  Lepton mu = Lepton("\\mu", muonMass);
  Lepton el = Lepton("e", electronMass);
  Lepton tau = Lepton("\\tau", tauMass);
  std::vector<Lepton> all_leptons = {el, mu, tau};

  Meson pi = Meson("\\pi^+", 139.57018, 130.2, MesonType::pseudoscalar, Charge::charged);
  Meson K = Meson("K^+", 493.677, 155.6, MesonType::pseudoscalar, Charge::charged);
  Meson D = Meson("D", 1869.62, 212, MesonType::pseudoscalar, Charge::charged);
  Meson Ds = Meson("D_s", 1968.47, 249, MesonType::pseudoscalar, Charge::charged);
  Meson B = Meson("B", 5279.29, 187, MesonType::pseudoscalar, Charge::charged);
  Meson Bc = Meson("B_c", 6275.1, 434, MesonType::pseudoscalar, Charge::charged);

  Meson pi0 = Meson("\\pi^0", 134.9766, 130.2, MesonType::pseudoscalar, Charge::neutral);
  Meson eta = Meson("\\eta", 547.862, 81.7, MesonType::pseudoscalar, Charge::neutral);
  Meson etaprime = Meson("\\eta'", 957.78, -94.7, MesonType::pseudoscalar, Charge::neutral);
  Meson etac = Meson("\\eta_c", 2983.6, 237, MesonType::pseudoscalar, Charge::neutral);

  Meson rho = Meson("\\rho", 775.11, 162000, MesonType::vector, Charge::charged);
  Meson Dstar = Meson("D^{\\ast+}", 2010.26, 535000, MesonType::vector, Charge::charged);
  Meson Dstars = Meson("D^{\\ast+}_s", 2112.1, 650000, MesonType::vector, Charge::charged);


  Double_t weinberg = 0.2223;
  Meson rho0 = Meson("\\rho^0", 775.26, 162000, MesonType::vector, Charge::neutral);
  rho0.insertValue("kh", 1.-2.*weinberg);
  Meson omega = Meson("\\omega", 782.65, 153000, MesonType::vector, Charge::neutral);
  omega.insertValue("kh", 4./3.*weinberg);
  Meson phi = Meson("\\phi", 1019.461, 234000, MesonType::vector, Charge::neutral);
  phi.insertValue("kh", (4./3.*weinberg)-1.);
  Meson jpsi = Meson("J/\\Psi", 3096.916, 1290000, MesonType::vector, Charge::neutral);
  jpsi.insertValue("kh", 1.-8./3.*weinberg);

  std::vector<Meson> mesons = {pi, K, D, Ds, B, Bc, pi0, eta, etaprime, etac, rho, Dstar, Dstars, rho0, omega, phi, jpsi};

  HNL N = HNL("20G", 20000, 3.55e-6, {el, mu, tau}); // 20 GeV but notated in MeV, allow mixing with all flavours (only possible implemtation currently is 1:1:1)

  Double_t Gamma = N.getTotalWidth(cfg, all_leptons, mesons);
  LOG_INFO("mass=" << N.getMass()/1000. << " GeV, " <<
           "ctau=" << gamma2ctau(cfg, Gamma) << " mm");

  N = HNL("20G", 20000, 3.55e-6, mu); // 20 GeV but notated in MeV, allow mixing with only muon flavour
  Gamma = N.getTotalWidth(cfg, all_leptons, mesons);
  LOG_INFO("mass=" << N.getMass()/1000. << " GeV, " <<
           "ctau=" << gamma2ctau(cfg, Gamma) << " mm");

  auto ch = N.getDecayChannels();
  for (auto it=ch.begin(); it!=ch.end(); ++it) {
    if(it->second > 0)
      std::cout << "\\Gamma(" << it->first << ") = \\SI{" << it->second << "}{}\\\\\n";
  }

  LOG_INFO("Searching now for U2 for 0.1 mm ctau... This might take a while!");
  Double_t foundU2 = ctauToU2(cfg, 0.1, all_leptons, mesons, N);
  LOG_INFO(std::endl << "FOUND U2: " << foundU2);

  std::vector<Meson> vector_charged;
  for(auto m : mesons) {
    if(m.getMesonType() != MesonType::vector) continue;
    if(m.getCharge() != Charge::charged) continue;
    vector_charged.emplace_back(m);
  }

  N.setMass(3000);
  LOG_INFO("pw: " << N.getPartialWidth(cfg, mu, rho));

  plot_meson_pw(cfg, mu, vector_charged, N, "mesons_vector_charged.pdf", 500, 5000);
  plot_I();
  plot_br(cfg, all_leptons, mesons, N, "BR.pdf", 1000, 5000);
  plot_qcd_correction("qcd_corr.pdf");


  return EXIT_SUCCESS;
}
