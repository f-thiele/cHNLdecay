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

int main() {
  auto cfg = std::make_shared<Config>(); // set BITS and initialize constants

  Double_t muonMass = 105.6583715; // MeV
  Double_t electronMass = 0.5109989461; // MeV
  Double_t tauMass = 1776.82; // MeV

  Lepton mu = Lepton("Muon", muonMass);
  Lepton el = Lepton("Electron", electronMass);
  Lepton tau = Lepton("Tau", tauMass);
  std::vector<Lepton> all_leptons = {el, mu, tau};

  Meson pi = Meson("Pi", 139.57018, 130.2, MesonType::pseudoscalar, Charge::charged);
  Meson K = Meson("K", 493.677, 155.6, MesonType::pseudoscalar, Charge::charged);
  Meson D = Meson("D", 1869.62, 212, MesonType::pseudoscalar, Charge::charged);
  Meson Ds = Meson("D_s", 1968.47, 249, MesonType::pseudoscalar, Charge::charged);
  Meson B = Meson("B", 5279.29, 187, MesonType::pseudoscalar, Charge::charged);
  Meson Bc = Meson("B_c", 6275.1, 434, MesonType::pseudoscalar, Charge::charged);

  Meson pi0 = Meson("Pi0", 134.9766, 130.2, MesonType::pseudoscalar, Charge::neutral);
  Meson eta = Meson("eta", 547.862, 81.7, MesonType::pseudoscalar, Charge::neutral);
  Meson etaprime = Meson("etaprime", 957.78, -94.7, MesonType::pseudoscalar, Charge::neutral);
  Meson etac = Meson("etac", 2983.6, 237, MesonType::pseudoscalar, Charge::neutral);

  Meson rho = Meson("rho", 775.11, 162000, MesonType::vector, Charge::charged);
  Meson Dstar = Meson("D star", 2010.26, 535000, MesonType::vector, Charge::charged);
  Meson Dstars = Meson("D star strange", 2112.1, 650000, MesonType::vector, Charge::charged);


  Double_t weinberg = 0.2223;
  Meson rho0 = Meson("rho0", 775.26, 162000, MesonType::vector, Charge::neutral);
  rho0.insertValue("kh", 1.-2.*weinberg);
  Meson omega = Meson("omega", 782.65, 153000, MesonType::vector, Charge::neutral);
  omega.insertValue("kh", 4./3.*weinberg);
  Meson phi = Meson("phi", 1019.461, 234000, MesonType::vector, Charge::neutral);
  phi.insertValue("kh", (4./3.*weinberg)-1.);
  Meson jpsi = Meson("jpsi", 3096.916, 1290000, MesonType::vector, Charge::neutral);
  jpsi.insertValue("kh", 1.-8./3.*weinberg);

  std::vector<Meson> mesons = {pi, K, D, Ds, B, Bc, pi0, eta, etaprime, etac, rho, Dstar, Dstars, rho0, omega, phi, jpsi};

  HNL N = HNL("20G", 20000, 1e-5, mu); // 20 GeV but notated in MeV
  TF1* f = new TF1("#alpha_{s}", qcd_coupling, 1, 100, 0);
  Double_t qcorr = qcd_correction(f->Eval(N.getMass()/1000.));
  std::cout << "QCD correction: " << qcorr << std::endl;


  Double_t totalWidth = 0;
  Double_t tw_lept = 0;
  Double_t tw_mes = 0;

  for(auto l1 : all_leptons) {
    for(auto l2 : all_leptons) {
      tw_lept += N.getPartialWidth(cfg, l1, l2);
    }
    for(auto m : mesons) {
      tw_mes += N.getPartialWidth(cfg, l1, m);
    }
  }

  totalWidth = qcorr*(tw_mes) + tw_lept;
  std::cout << "mass=" << N.getMass()/1000. << " GeV, "
            << "ctau=" << gamma2ctau(cfg, totalWidth) << " mm" << std::endl;

  std::vector<Meson> vector_charged;
  for(auto m : mesons) {
    if(m.getMesonType() != MesonType::vector) continue;
    if(m.getCharge() != Charge::charged) continue;
    vector_charged.emplace_back(m);
  }
  plot_meson_pw(cfg, mu, vector_charged, N, "mesons_vector_charged.pdf", 500, 5000);


  return EXIT_SUCCESS;
}
