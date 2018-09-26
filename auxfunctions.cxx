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

#include "auxfunctions.h"
#include "Config.h"
#include "HNL.h"
#include "Lepton.h"
#include "Logger.h"
#include "Math/GaussIntegrator.h"
#include "Math/WrappedTF1.h"
#include "Meson.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLegend.h"
#include <fstream>
#include <sstream>
#include <iostream>

void kaellen(std::shared_ptr<Config> cfg, mpfr_t result, mpfr_t a, mpfr_t b,
             mpfr_t c) {
  unsigned int BITS = cfg->getBITS();

  mpfr_t temp;
  mpfr_init2(temp, BITS);
  mpfr_pow_ui(result, a, 2, MPFR_RNDD);
  mpfr_pow_ui(temp, b, 2, MPFR_RNDD);
  mpfr_add(result, result, temp, MPFR_RNDD);
  mpfr_pow_ui(temp, c, 2, MPFR_RNDD);
  mpfr_add(result, result, temp, MPFR_RNDD);
  mpfr_mul(temp, a, b, MPFR_RNDD);
  mpfr_mul_ui(temp, temp, 2, MPFR_RNDD);
  mpfr_sub(result, result, temp, MPFR_RNDD);
  mpfr_mul(temp, a, c, MPFR_RNDD);
  mpfr_mul_ui(temp, temp, 2, MPFR_RNDD);
  mpfr_sub(result, result, temp, MPFR_RNDD);
  mpfr_mul(temp, b, c, MPFR_RNDD);
  mpfr_mul_ui(temp, temp, 2, MPFR_RNDD);
  mpfr_sub(result, result, temp, MPFR_RNDD);
  mpfr_clear(temp);
}

Double_t kaellen(Double_t a, Double_t b, Double_t c) {
  return a * a + b * b + c * c - 2 * a * b - 2 * a * c - 2 * b * c;
}

Double_t int__i(Double_t *x, Double_t *par) {
  Double_t xx = x[0];
  Double_t xu = par[0];
  Double_t xd = par[1];
  Double_t xl = par[2];
  // Double_t _s = (TMath::Power(xx, 2)+TMath::Power(xl, 4)+TMath::Power(xd,
  // 4)-2*xx*xl*xl-2*xx*xd*xd-2*xl*xl*xd*xd)*(TMath::Power(1,
  // 2)+TMath::Power(xx, 2)+TMath::Power(xu, 4)-2*xx-2*xu*xu-2*xx*xu*xu);
  Double_t f =
      (xx - xl * xl - xd * xd) / xx * (1 + xu * xu - xx) *
      TMath::Sqrt(kaellen(xx, TMath::Power(xl, 2), TMath::Power(xd, 2)) *
                  kaellen(1, xx, TMath::Power(xu, 2)));

  return f;
}
Double_t I_xu_xd_xl(Double_t xu, Double_t xd, Double_t xl) {
  TF1 *f = new TF1("I(x_{d},x_{u},x_{l}) function", int__i, 0, 1, 3);
  f->SetParNames("xu", "xd", "xl");
  f->SetParameters(xu, xd, xl);

  ROOT::Math::WrappedTF1 wf1(*f);

  int status = 1;
  Double_t val = -1;
  while (status != 0) {
    ROOT::Math::GaussIntegrator ig;

    ig.SetFunction(wf1);
    ig.SetRelTolerance(0.0001);
    val = ig.Integral(TMath::Power(xd + xl, 2), TMath::Power(1 - xu, 2));
    status = ig.Status();
  }
  delete f;

  return 12 * val;
}
Double_t gamma2ctau(std::shared_ptr<Config> cfg, Double_t gamma) {
  // gamma is expected in [MeV]
  // output ctau is in [mm]

  mpfr_t fermiC, fermiCsq, pi, VUDsq, SOL, HBAR;
  unsigned int BITS = cfg->getBITS();
  mpfr_init2(fermiC, BITS);
  mpfr_init2(fermiCsq, BITS);
  mpfr_init2(pi, BITS);
  mpfr_init2(VUDsq, BITS);
  mpfr_init2(SOL, BITS);
  mpfr_init2(HBAR, BITS);
  cfg->getFermiCsq(fermiCsq);
  cfg->getFermiC(fermiC);
  cfg->getPi(pi);
  cfg->getVUDsq(VUDsq);
  cfg->getSOL(SOL);
  cfg->getHBAR(HBAR);
  mpfr_t _gamma, result;
  mpfr_init2(_gamma, BITS);
  mpfr_init2(result, BITS);

  mpfr_mul(result, SOL, HBAR, MPFR_RNDD);      // eV s m /s = [eV m]
  mpfr_mul_ui(result, result, 1e3, MPFR_RNDD); // [eV mm]
  mpfr_set_d(_gamma, gamma, MPFR_RNDD);        // given in [MeV]
  mpfr_mul_ui(_gamma, _gamma, 1e6, MPFR_RNDD); // now in [eV]
  mpfr_div(result, result, _gamma, MPFR_RNDD); // [mm]

  Double_t rval = mpfr_get_d(result, MPFR_RNDD);
  mpfr_clears(fermiC, fermiCsq, pi, VUDsq, SOL, HBAR, _gamma, result,
              (mpfr_ptr)0);

  return rval;
}

Double_t qcd_coupling(Double_t *x, Double_t *par) {
  UInt_t nc = 3;
  UInt_t nf = 6;
  Double_t beta0 = (11. * nc - 2. * nf) / (12 * TMath::Pi());

  // Double_t alpha_mu = 0.1184;
  // Double_t mu = 91.1876;
  Double_t alpha_mu = 0.31593;
  Double_t mu = 1.776;
  return alpha_mu /
         (1. + beta0 * alpha_mu * TMath::Log(x[0] * x[0] / (mu * mu)));
}

Double_t qcd_coupling(Double_t mass) {
  TF1 *f = new TF1("#alpha_{s}", qcd_coupling, 1, 100, 0);
  Double_t ret = f->Eval(mass);
  delete f;
  return ret;
}

Double_t qcd_correction(Double_t alpha) {
  return alpha / TMath::Pi() + 5.2 * TMath::Power(alpha / TMath::Pi(), 2) +
         26.4 * TMath::Power(alpha / TMath::Pi(), 3);
}

Double_t f_qcd_correction(Double_t *x, Double_t *par) {
  TF1 *f = new TF1("#alpha_{s}", qcd_coupling, 1, 100, 0);
  Double_t alpha = f->Eval(x[0]);
  delete f;

  return alpha / TMath::Pi() + 5.2 * TMath::Power(alpha / TMath::Pi(), 2) +
         26.4 * TMath::Power(alpha / TMath::Pi(), 3);
}

TGraph *create_graph(Double_t xd, Double_t xl, Float_t low, Float_t high,
                     Float_t stepsize) {
  std::vector<Double_t> res_x;
  std::vector<Double_t> res_y;
  std::vector<Double_t> res_n;
  for (Double_t x = low; x <= high; x = x + stepsize) {
    Double_t _temp = I_xu_xd_xl(x, xd, xl);
    if (TMath::IsNaN(_temp))
      continue;
    Double_t _temp0 = I_xu_xd_xl(x, 0, 0);
    if (TMath::IsNaN(_temp0))
      continue;
    res_x.emplace_back(x);
    res_y.emplace_back(_temp);
    res_n.emplace_back(_temp0);
  }

  for (UInt_t i = 0; i < res_y.size(); i++) {
    Double_t ratio = res_y.at(i) / res_n.at(i);
    if (i > 0 && (ratio - res_y.at(i - 1)) / ratio > 0.5)
      continue;
    res_y.at(i) = res_y.at(i) / res_n.at(i);
  }

  TGraph *g = new TGraph(res_x.size(), &(res_x[0]), &(res_y[0]));

  return g;
}

Double_t ctauToU2(std::shared_ptr<Config> cfg, Double_t target,
                  const std::vector<Lepton> &leptons,
                  const std::vector<Meson> &mesons, HNL &N, Double_t start,
                  Double_t tol, Double_t stepsize) {
  N.setAngle(start);
  Double_t ctau = gamma2ctau(cfg, N.getTotalWidth(cfg, leptons, mesons));
  Double_t found_angle = start;
  UInt_t iterations = 0;
  Double_t prev_ctau = ctau;

  for (Double_t angle = start; abs(target - ctau) > tol; angle *= stepsize) {
    iterations++;
    N.setAngle(angle);
    prev_ctau = ctau;
    ctau = gamma2ctau(cfg, N.getTotalWidth(cfg, leptons, mesons));
    if (abs(target - ctau) > abs(target - prev_ctau)) {
      LOG_WARNING("Stopped at angle="
                  << angle << " with ctau=" << ctau
                  << " as distance increased to previous ctau=" << prev_ctau);
      break;
    }
    found_angle = angle;
    if (iterations % 1000) {
      std::cout << "\rangle=" << angle << ", ctau=" << ctau << std::flush;
    }
  }
  return found_angle;
}

TString pdgIdToLaTeX(Int_t p) {
  std::map<Int_t, TString> label = {
      {11, "e^-"},           {12, "\\nu_e"},   {13, "\\mu^-"},
      {14, "\\nu_\\mu"},     {15, "\\tau^-"},  {16, "\\nu_\\tau"},
      {211, "\\pi^+"},       {321, "K^+"},     {411, "D^+"},
      {431, "D_s^+"},        {521, "B^+"},     {541, "B_c^+"},
      {111, "\\pi^0"},       {221, "\\eta"},   {331, "\\eta'"},
      {441, "\\eta_c"},      {213, "\\rho^+"}, {413, "D^{\\ast+}"},
      {431, "D^{\\ast+}_s"}, {113, "\\rho^0"}, {223, "\\omega"},
      {333, "\\phi"},        {443, "J/\\Psi"}};

  if (p > 0)
    return label.at(p);
  else if (p < 0) {
    TString ret = label.at(-p);
    if (ret.Contains("+"))
      ret.ReplaceAll("+", "-");
    else if (ret.Contains("-"))
      ret.ReplaceAll("-", "+");
    else
      ret = TString("\\bar{" + ret + "}");

    return ret;
  }
  return TString("");
}

std::vector<std::vector<Double_t> > parseFile(std::string name) {
  std::vector<std::vector<Double_t> >     data;

  std::ifstream file(name);

  std::string line;

  while(std::getline(file, line))
    {
      std::vector<Double_t>   lineData;
      std::stringstream  lineStream(line);

      Double_t value;
      while(lineStream >> value)
        {
          lineData.push_back(value);
        }
      data.push_back(lineData);
    }

  return data;
}

void plot_I() {
  TCanvas *c1 = new TCanvas("c1", "c1", 500, 400);
  TGraph *g = create_graph(0.5, 0, 0.1);
  TGraph *h = create_graph(0.25, 0.25, 0.1);
  TGraph *i = create_graph(0.1, 0.4, 0.1);
  g->SetLineColor(kBlue);
  g->SetLineWidth(2);
  h->SetLineColor(kRed);
  h->SetLineWidth(2);
  h->SetLineStyle(9);
  i->SetLineColor(kGreen);
  i->SetLineWidth(2);
  i->SetLineStyle(7);

  TLegend *leg = new TLegend(0.15, 0.15, 0.6, 0.3);
  leg->AddEntry(g, "u_{d} = 0.50, u_{l} = 0.00");
  leg->AddEntry(h, "u_{d} = 0.25, u_{l} = 0.25");
  leg->AddEntry(i, "u_{d} = 0.10, u_{l} = 0.40");
  c1->SetLogy();
  c1->SetLogx();

  g->Draw("AL");
  g->GetYaxis()->SetRangeUser(0.001, 1);
  g->GetXaxis()->SetRangeUser(0.1, 0.5);
  h->Draw("SAME");
  i->Draw("SAME");
  leg->Draw();
  c1->SaveAs("I_xu_xd_xl.pdf");
  delete c1;
  delete g;
  delete h;
  delete i;
  delete leg;
}
