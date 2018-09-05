#include "Lepton.h"
#include "Meson.h"
#include "HNL.h"
#include "Config.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include <iostream>
#include <map>
#include <vector>

void plot_meson_pw(std::shared_ptr<Config> cfg, Lepton alpha, std::vector<Meson> mesons, HNL N, TString output, Int_t lowMass, Int_t highMass, Int_t stepsize=10) {
  std::vector<Double_t> res_m;

  std::map<TString, std::vector<Double_t>> meson_pws;

  for(auto m :mesons) {
    meson_pws.insert(std::pair<TString, std::vector<Double_t>>(m.getName(), {}));
  }

  for(Int_t mass = lowMass; mass < highMass; mass+=stepsize) {
    N.setMass(mass);
    for(auto m : mesons) {
      Double_t pw =N.getPartialWidth(cfg, alpha, m);
      mpfr_t fermiC, fermiCsq, pi, VUDsq;
      unsigned int BITS = cfg->getBITS();
      mpfr_init2(fermiC, BITS);
      mpfr_init2(fermiCsq, BITS);
      mpfr_init2(pi, BITS);
      mpfr_init2(VUDsq, BITS);
      cfg->getFermiCsq(fermiCsq);
      cfg->getFermiC(fermiC);
      cfg->getPi(pi);
      cfg->getVUDsq(VUDsq);

      mpfr_t temp,factor, mesonMass, HNLmass, angle, gh;

      mpfr_init2(temp, BITS);
      mpfr_init2(factor, BITS);
      mpfr_init2(mesonMass, BITS);
      mpfr_init2(HNLmass, BITS);
      mpfr_init2(angle, BITS);
      mpfr_init2(gh, BITS);

      mpfr_set_d(mesonMass, m.getMass(), MPFR_RNDD);
      mpfr_set_d(HNLmass, N.getMass(), MPFR_RNDD);
      mpfr_set_d(angle, N.getAngle(), MPFR_RNDD);
      mpfr_set_d(gh, m.getDecayConstant(), MPFR_RNDD);

      // calculate first factors NZ * GF^2 * M^5 * U^2/(192*pi^3)
      mpfr_pow_ui(factor, HNLmass, 3, MPFR_RNDD);
      mpfr_mul(factor, factor, angle, MPFR_RNDD);
      mpfr_mul(factor, factor, VUDsq, MPFR_RNDD);
      mpfr_pow_ui(temp, gh, 2, MPFR_RNDD);
      mpfr_mul(factor, factor, temp, MPFR_RNDD);
      mpfr_mul(factor, factor, fermiCsq, MPFR_RNDD);
      mpfr_mul_ui(temp, mesonMass, 2, MPFR_RNDD);
      mpfr_mul(temp, temp, pi, MPFR_RNDD);
      mpfr_mul_ui(temp, temp, 16, MPFR_RNDD);
      mpfr_div(factor,factor,temp, MPFR_RNDD);
      // end calculation first factors
      Double_t Gamma0 =mpfr_get_d(factor, MPFR_RNDD);
      // std::cout << "mass: " << mass << " -- " << pw/Gamma0 << std::endl;
      meson_pws.at(m.getName()).emplace_back(pw/Gamma0);
    }

    res_m.emplace_back(mass);
  }

  TCanvas* c1 = new TCanvas("c1", "c1", 500, 400);
  TLegend* leg = new TLegend(0.15, 0.15, 0.4, 0.4);
  std::vector<TGraph*> graphs;
  Color_t style = 1;
  for(auto m : mesons) {
    TGraph* g = new TGraph(res_m.size(), &(res_m[0]), &(meson_pws.at(m.getName())[0]));
    g->SetLineColor(style);
    style++;
    leg->AddEntry(g, m.getName());
    graphs.emplace_back(g);
  }
  bool first = true;
  for(auto g : graphs) {
    if(first) {
      g->Draw("AL");
      first = false;
    } else {
      g->Draw("SAME");
    }
  }
  leg->Draw();
  c1->SetLogy();
  c1->SaveAs(output);
}
