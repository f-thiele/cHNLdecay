//  cHNLdecay  --  calculate decay widths of Heavy Neutral Leptons
//  Copyright (C) 2018 - Fabian A.J. Thiele, <fabian.thiele@posteo.de>
//
//  This file is part of cHNLdecay.
//
//  cHNLdecay is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  cHNLdecay is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include "Config.h"
#include "HNL.h"
#include "Lepton.h"
#include "Logger.h"
#include "Meson.h"
#include "ParticleCatalogue.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLegend.h"
#include "auxfunctions.h"
#include <iostream>
#include <map>
#include <vector>

void plot_meson_pw(std::shared_ptr<Config> cfg, Lepton alpha,
                   std::vector<Meson> mesons, HNL N, TString output,
                   Int_t lowMass, Int_t highMass, Int_t stepsize = 10) {
  std::vector<Double_t> res_m;

  std::map<TString, std::vector<Double_t>> meson_pws;

  for (auto m : mesons) {
    meson_pws.insert(
        std::pair<TString, std::vector<Double_t>>(m.getName(), {}));
  }

  for (Int_t mass = lowMass; mass <= highMass; mass += stepsize) {
    N.setMass(mass);
    for (auto m : mesons) {
      LOG_DEBUG("Meson: " << m.getName());
      Double_t pw = N.getPartialWidth(cfg, alpha, m);
      mpfr_t fermiC, fermiCsq, pi, VUDsq;
      unsigned int BITS = cfg->getBITS();
      mpfr_init2(fermiC, BITS);
      mpfr_init2(fermiCsq, BITS);
      mpfr_init2(pi, BITS);
      mpfr_init2(VUDsq, BITS);
      cfg->getFermiCsq(fermiCsq);
      cfg->getFermiC(fermiC);
      cfg->getPi(pi);
      cfg->getVUDsq(VUDsq, m);

      mpfr_t temp, factor, mesonMass, HNLmass, angle, ghsq;

      mpfr_init2(temp, BITS);
      mpfr_init2(factor, BITS);
      mpfr_init2(mesonMass, BITS);
      mpfr_init2(HNLmass, BITS);
      mpfr_init2(angle, BITS);
      mpfr_init2(ghsq, BITS);

      mpfr_set_d(mesonMass, m.getMass(), MPFR_RNDD);
      mpfr_set_d(HNLmass, N.getMass(), MPFR_RNDD);
      mpfr_set_d(angle, N.getAngle(), MPFR_RNDD);
      mpfr_set_d(ghsq, m.getDecayConstant(), MPFR_RNDD);
      mpfr_pow_ui(ghsq, ghsq, 2, MPFR_RNDD);

      // calculate first factors NZ * GF^2 * M^5 * U^2/(192*pi^3)
      mpfr_pow_ui(factor, HNLmass, 3, MPFR_RNDD);
      mpfr_mul(factor, factor, angle, MPFR_RNDD);
      mpfr_mul(factor, factor, VUDsq, MPFR_RNDD);
      mpfr_mul(factor, factor, ghsq, MPFR_RNDD);
      mpfr_mul(factor, factor, fermiCsq, MPFR_RNDD);
      mpfr_pow_ui(temp, mesonMass, 2, MPFR_RNDD);
      mpfr_mul(temp, temp, pi, MPFR_RNDD);
      mpfr_mul_ui(temp, temp, 16, MPFR_RNDD);
      mpfr_div(factor, factor, temp, MPFR_RNDD);

      // end calculation first factors
      Double_t Gamma0 = mpfr_get_d(factor, MPFR_RNDD);

      mpfr_clears(fermiC, fermiCsq, pi, VUDsq, mesonMass, HNLmass, angle,
                  factor, temp, ghsq, (mpfr_ptr)0);
      LOG_DEBUG("Gamma0:" << Gamma0);
      LOG_DEBUG("pw:" << pw);
      LOG_DEBUG("mass: " << mass << " -- " << pw / Gamma0);
      meson_pws.at(m.getName()).emplace_back(pw / Gamma0);
    }

    res_m.emplace_back(mass);
  }

  TCanvas *c1 = new TCanvas("c1", "c1", 500, 400);
  TLegend *leg = new TLegend(0.15, 0.15, 0.4, 0.4);
  std::vector<TGraph *> graphs;
  Color_t style = 1;
  for (auto m : mesons) {
    TGraph *g =
        new TGraph(res_m.size(), &(res_m[0]), &(meson_pws.at(m.getName())[0]));
    g->SetLineColor(style);
    style++;
    leg->AddEntry(g, alpha.getName() + m.getName());
    graphs.emplace_back(g);
  }
  bool first = true;
  for (auto g : graphs) {
    if (first) {
      g->GetYaxis()->SetRangeUser(0.01, 1.1);
      g->Draw("AL");
      first = false;
    } else {
      g->Draw("SAME");
    }
  }
  leg->Draw();
  c1->SetLogy();
  c1->SaveAs(output);
  delete c1;
  delete leg;
  for (auto g : graphs) {
    delete g;
  }
}

void plot_br(std::shared_ptr<Config> cfg, std::vector<Lepton> leptons,
             std::vector<Meson> mesons, HNL N, TString output, Int_t lowMass,
             Int_t highMass, Int_t stepsize) {
  std::vector<Double_t> res_m;
  std::vector<Double_t> res_mes;
  std::vector<Double_t> res_lep;
  std::vector<Double_t> res_inv;

  for (Int_t mass = lowMass; mass < highMass; mass += stepsize) {
    N.setMass(mass);

    Double_t tw_mes = 0;
    Double_t tw_lep = 0;
    Double_t tw_inv = 0;
    for (auto l1 : leptons) {
      for (auto l2 : leptons) {
        tw_lep += N.getPartialWidth(cfg, l1, l2, false);
        tw_inv += N.getPartialWidthInv(cfg, l1, l2);
      }
      for (auto m : mesons) {
        tw_mes += N.getPartialWidth(cfg, l1, m);
      }
    }
    LOG_DEBUG("mass: " << N.getMass() << ", meson pw: " << tw_mes);
    res_m.emplace_back(mass);
    Double_t tot = tw_mes + tw_lep + tw_inv;
    res_mes.emplace_back(tw_mes / tot);
    res_lep.emplace_back(tw_lep / tot);
    res_inv.emplace_back(tw_inv / tot);
  }

  TCanvas *c1 = new TCanvas("c1", "c1", 500, 400);

  TGraph *g_mes = new TGraph(res_m.size(), &(res_m[0]), &(res_mes[0]));
  TGraph *g_lep = new TGraph(res_m.size(), &(res_m[0]), &(res_lep[0]));
  TGraph *g_inv = new TGraph(res_m.size(), &(res_m[0]), &(res_inv[0]));

  g_mes->SetLineColor(kBlue);
  g_mes->SetLineWidth(2);
  g_lep->SetLineColor(kRed);
  g_lep->SetLineWidth(2);
  g_lep->SetLineStyle(9);
  g_inv->SetLineColor(kGreen);
  g_inv->SetLineWidth(2);
  g_inv->SetLineStyle(7);

  TLegend *leg = new TLegend(0.15, 0.15, 0.6, 0.3);
  leg->AddEntry(g_mes, "quarks");
  leg->AddEntry(g_lep, "leptons");
  leg->AddEntry(g_inv, "invisible");
  // c1->SetLogy();
  c1->SetLogx();

  g_mes->Draw("AL");
  g_mes->GetYaxis()->SetRangeUser(0.05, 1);
  g_mes->GetXaxis()->SetRangeUser(1000, 5000);
  g_lep->Draw("SAME");
  g_inv->Draw("SAME");
  leg->Draw();
  c1->SaveAs(output);
  delete leg;
  delete c1;
  delete g_mes;
  delete g_lep;
  delete g_inv;
}

void plot_qcd_coupling(TString output) {
  TF1 *f = new TF1("#alpha_{s}", qcd_coupling, 1, 100, 0);
  f->SetMinimum(0);
  TCanvas *c1 = new TCanvas("c1", "c1", 500, 400);
  f->Draw();
  c1->SaveAs(output);
  delete c1;
  delete f;
}

void plot_qcd_correction(TString output) {
  TF1 *f = new TF1("#Delta_{QCD}", f_qcd_correction, 1, 5, 0);
  f->SetMinimum(0);
  f->SetMaximum(0.35);
  TCanvas *c1 = new TCanvas();
  f->Draw();
  c1->SaveAs(output);
  delete f;
  delete c1;
}

void plot_br_low(std::shared_ptr<Config> cfg, std::vector<Lepton> leptons,
                 std::vector<Meson> mesons, HNL N, TString output,
                 Int_t lowMass, Int_t highMass, Int_t stepsize) {
  ParticleCatalogue pc;

  std::map<TString, std::vector<Meson>> plot_sep;
  std::vector<Meson> v_pi = {
      pc.getMeson(211), // pi+
      pc.getMeson(111)  // pi0
  };
  plot_sep["pi"] = v_pi;

  std::vector<Meson> v_K = {pc.getMeson(321)}; // K+
  plot_sep["K"] = v_K;

  std::vector<Meson> v_eta = {
      pc.getMeson(221), // eta
      pc.getMeson(331), // eta prime
      pc.getMeson(441)  // eta c
  };
  plot_sep["eta"] = v_eta;

  std::vector<Meson> v_rho = {
      pc.getMeson(213), // rho+
      pc.getMeson(113)  // rho0
  };
  plot_sep["rho"] = v_rho;

  std::vector<Double_t> res_m;

  std::map<TString, TGraph *> graphs;
  std::map<TString, std::vector<Double_t>> res_values;
  std::map<TString, Double_t> tws;

  for (auto entry : plot_sep) {
    res_values[entry.first] = {};
    graphs[entry.first] = nullptr;
    tws[entry.first] = 0;
  }

  std::vector<Double_t> res_lep;
  std::vector<Double_t> res_inv;

  for (Int_t mass = lowMass; mass < highMass; mass += stepsize) {
    N.setMass(mass);
    for (auto entry : plot_sep) {
      tws[entry.first] = 0; // zero the accoumulated widths for new mass
    }

    Double_t tw_mes = 0;
    Double_t tw_lep = 0;
    Double_t tw_inv = 0;
    for (auto l1 : leptons) {
      for (auto l2 : leptons) {
        tw_lep += N.getPartialWidth(cfg, l1, l2, false);
        tw_inv += N.getPartialWidthInv(cfg, l1, l2);
      }
      for (auto m : mesons) {
        Double_t t_pw = N.getPartialWidth(cfg, l1, m);
        tw_mes += t_pw;
        for (auto entry : plot_sep) {
          for (auto c_m : entry.second) {
            if (m == c_m)
              tws.at(entry.first) += t_pw;
          }
        }
      }
    }
    res_m.emplace_back(mass);
    Double_t tot = tw_mes + tw_lep + tw_inv;

    for (auto entry : plot_sep)
      res_values.at(entry.first).emplace_back(tws.at(entry.first) / tot);

    res_lep.emplace_back(tw_lep / tot);
    res_inv.emplace_back(tw_inv / tot);
  }

  TCanvas *c1 = new TCanvas("c1", "c1", 500, 400);

  Short_t lc = 1;

  for (auto entry : plot_sep) {
    graphs.at(entry.first) =
        new TGraph(res_m.size(), &(res_m[0]), &(res_values.at(entry.first)[0]));
    graphs.at(entry.first)->SetLineColor(lc);
    graphs.at(entry.first)->SetLineWidth(2);
    lc++;
  }

  TGraph *g_lep = new TGraph(res_m.size(), &(res_m[0]), &(res_lep[0]));
  TGraph *g_inv = new TGraph(res_m.size(), &(res_m[0]), &(res_inv[0]));

  g_lep->SetLineColor(kRed);
  g_lep->SetLineWidth(2);
  g_lep->SetLineStyle(9);
  g_inv->SetLineColor(kGreen);
  g_inv->SetLineWidth(2);
  g_inv->SetLineStyle(7);

  TLegend *leg = new TLegend(0.15, 0.15, 0.6, 0.3);
  leg->AddEntry(g_lep, "leptons");
  leg->AddEntry(g_inv, "invisible");
  for (auto entry : plot_sep) {
    leg->AddEntry(graphs.at(entry.first), entry.first);
  }
  c1->SetLogy();
  c1->SetLogx();

  bool isFirst = true;
  for (auto entry : plot_sep) {
    if (isFirst) {
      graphs.at(entry.first)->Draw("AL");
      graphs.at(entry.first)->GetYaxis()->SetRangeUser(0.0001, 1);
      isFirst = false;
    } else {
      graphs.at(entry.first)->Draw("SAME");
    }
  }
  g_lep->Draw("SAME");
  g_inv->Draw("SAME");
  leg->Draw();
  c1->SaveAs(output);
  delete leg;
  delete c1;
  for (auto entry : plot_sep) {
    delete graphs.at(entry.first);
  }
}
