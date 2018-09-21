// Copyright (C) 2018 - Fabian A.J. Thiele, <fabian.thiele@cern.ch>

#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObjString.h"

#include "HNL.h"
#include "Lepton.h"
#include "Meson.h"

#include "Config.h"
#include "Logger.h"
#include "ParticleCatalogue.h"
#include "auxfunctions.h"
#include "partialWidths.h"
#include "plots.h"
#include <fstream>
#include <getopt.h>
#include <gmp.h>
#include <iomanip>
#include <iostream>
#include <mpfr.h>
#include <stdlib.h>

Level gLOGLEVEL;
int main(int argc, char **argv) {
  auto cfg = std::make_shared<Config>(); // set BITS and initialize constants

  // initialize the three variables that define our HNL
  Int_t HNLmass = 0;
  Double_t angle = 0;
  std::vector<Lepton> mixes_with;

  // sensible default for loglevel
  gLOGLEVEL = Level::INFO;

  ParticleCatalogue pc;

  std::vector<Lepton> all_leptons = pc.getAllLeptons();
  std::vector<Meson> mesons = pc.getAllMesons();

  while (1) {
    int c;
    static struct option long_options[] = {
        /* These options set a flag. */
        {"loglevel", required_argument, 0, 'l'},
        {"generations", required_argument, 0, 'g'},
        {"angle", required_argument, 0, 'u'},
        {"mass", required_argument, 0, 'm'},
        {0, 0, 0, 0}};

    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long(argc, argv, "v:l:u:m:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c) {
    case 0:
      /* If this option set a flag, do nothing else now. */
      if (long_options[option_index].flag != 0)
        break;
      printf("option %s", long_options[option_index].name);
      if (optarg)
        printf(" with arg %s", optarg);
      printf("\n");
      break;

    case 'l': {
      TString x = TString(optarg);
      if (x.CompareTo("debug") == 0) {
        gLOGLEVEL = Level::DEBUG;

      } else if (x.CompareTo("warning") == 0) {
        gLOGLEVEL = Level::WARNING;

      } else if (x.CompareTo("error") == 0) {
        gLOGLEVEL = Level::ERROR;

      } else {
        gLOGLEVEL = Level::INFO;
      }
    }

    break;

    case 'g': {
      TString x = TString(optarg);
      TObjArray *tx = x.Tokenize(",");

      for (Int_t i = 0; i < tx->GetEntries(); i++) {
        TString iter_str = ((TObjString *)tx->At(i))->CopyString();

        for (auto l : all_leptons) {
          if (iter_str.Atoi() == l.getPdgId())
            mixes_with.emplace_back(l);
        }
      }
    }

    break;

    case 'u':
      angle = std::atof(optarg);
      break;

    case 'm':
      HNLmass = std::atoi(optarg);
      break;

    case '?':
      /* getopt_long already printed an error message. */
      break;

    default:
      abort();
    }
  }

  /* Print any remaining command line arguments (not options). */
  if (optind < argc) {
    printf("non-option ARGV-elements: ");
    while (optind < argc)
      printf("%s ", argv[optind++]);
    putchar('\n');
  }

  for (auto l : mixes_with) {
    LOG_DEBUG("HNL mixes with " << l.getName());
  }
  // this is the HNL we configure via command line options. If it mixes with
  // multiple generations it can only assume mixing angle ratios 1:1:..:1
  HNL N = HNL("HNL", HNLmass, angle, mixes_with);

  Double_t Gamma = N.getTotalWidth(cfg, all_leptons, mesons);
  LOG_INFO("mass=" << N.getMass() / 1000. << " GeV, "
                   << "ctau=" << gamma2ctau(cfg, Gamma) << " mm");

  auto ch = N.getDecayChannels();
  Int_t longest_channel = 0;
  for (auto it = ch.begin(); it != ch.end(); ++it) {
    if (longest_channel < (it->first).size()) {
      longest_channel = (it->first).size();
    }
  }

  std::ofstream lxf("decaychannels.txt");
  for (auto it = ch.begin(); it != ch.end(); ++it) {
    if (it->second > 0) {
      lxf << "\\Gamma(";
      for (Int_t index = 0; index < longest_channel; ++index) {
        if (index < (it->first).size()) {
          Int_t p = (it->first).at(index);
          std::cout << std::setw(9 * (1 + index)) << p;
          lxf << std::setw(9 * (1 + index)) << pdgIdToLaTeX(p);
        } else {
          std::cout << std::setw(9 * (1 + index)) << " ";
          lxf << std::setw(9 * (1 + index)) << " ";
        }
      }
      lxf << ") = \\SI{";
      std::cout << std::setw(9 * longest_channel) << it->second << std::endl;
      lxf << std::setw(9 * longest_channel) << it->second;
      lxf << "}{\\MeV} \\\\" << std::endl;
    }
  }
  lxf.close();

  plot_br_low(cfg, all_leptons, mesons, N, "BR_low.pdf", 50, 1000, 1);

  std::vector<Meson> vector_charged;
  for (auto m : mesons) {
    if (m.getMesonType() != MesonType::vector)
      continue;
    if (m.getCharge() != Charge::charged)
      continue;
    vector_charged.emplace_back(m);
  }

  TF1 *f = new TF1("#Delta_{QCD}", qcd_coupling, 1, 100, 0);
  Double_t qcd_corr = qcd_correction(f->Eval(1776 / 1000.));
  LOG_INFO("QCD correction at tau mass: " << qcd_corr);
  LOG_INFO(
      "QCD correction at 20 GeV: " << qcd_correction(f->Eval(20000 / 1000.)));
  LOG_INFO("QCD coupling at tau: " << qcd_coupling(1776. / 1000.));
  LOG_INFO("QCD coupling at Z mass: " << qcd_coupling(91.2));
  delete f;

  plot_I();
  plot_br(cfg, all_leptons, mesons, N, "BR.pdf", 1000, 5000);

  plot_qcd_correction("qcd_corr.pdf");
  plot_qcd_coupling("qcd_coupl.pdf");

  bool isMajorana = N.isMajorana();
  N.setMajorana(false); // disable majorana here because we want to compare the
                        // pw for one charge configuration only
  plot_meson_pw(cfg, pc.getLepton(13), vector_charged, N,
                "mesons_vector_charged.pdf", 500, 5000);
  N.setMajorana(isMajorana); // restore majorana value

  // EXAMPLE for further checks and configurations
  // ----------
  // N.setMass(3000);
  // LOG_INFO("pw: " << N.getPartialWidth(cfg, mu, rho));
  // LOG_INFO("Searching now for U2 for 0.1 mm ctau... This might take a
  // while!"); Double_t foundU2 = ctauToU2(cfg, 0.1, all_leptons, mesons, N);
  // LOG_INFO(std::endl << "FOUND U2: " << foundU2);

  return EXIT_SUCCESS;
}
