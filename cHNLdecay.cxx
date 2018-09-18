// Copyright (C) 2018 - Fabian A.J. Thiele, <fabian.thiele@cern.ch>

#include "TGraph.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TMath.h"
#include "TF1.h"
#include "TObjString.h"

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
#include <stdlib.h>
#include <getopt.h>

Level gLOGLEVEL;
int main (int argc, char **argv)
{
  auto cfg = std::make_shared<Config>(); // set BITS and initialize constants

  // initialize the three variables that define our HNL
  Int_t HNLmass = 0;
  Double_t angle = 0;
  std::vector<Lepton> mixes_with;

  // sensible default for loglevel
  gLOGLEVEL = Level::INFO;

  Double_t muonMass = 105.6583715; // MeV
  Double_t electronMass = 0.5109989461; // MeV
  Double_t tauMass = 1776.82; // MeV

  Lepton mu = Lepton(13, "\\mu", muonMass);
  Lepton el = Lepton(11, "e", electronMass);
  Lepton tau = Lepton(15, "\\tau", tauMass);
  std::vector<Lepton> all_leptons = {el, mu, tau};

  Meson pi = Meson(211,"\\pi^+", 139.57018, 130.2, MesonType::pseudoscalar, Charge::charged);
  Meson K = Meson(321, "K^+", 493.677, 155.6, MesonType::pseudoscalar, Charge::charged);
  Meson D = Meson(411, "D", 1869.62, 212, MesonType::pseudoscalar, Charge::charged);
  Meson Ds = Meson(431, "D_s", 1968.47, 249, MesonType::pseudoscalar, Charge::charged);
  Meson B = Meson(521, "B", 5279.29, 187, MesonType::pseudoscalar, Charge::charged);
  Meson Bc = Meson(541, "B_c", 6275.1, 434, MesonType::pseudoscalar, Charge::charged);

  Meson pi0 = Meson(111, "\\pi^0", 134.9766, 130.2, MesonType::pseudoscalar, Charge::neutral);
  Meson eta = Meson(221, "\\eta", 547.862, 81.7, MesonType::pseudoscalar, Charge::neutral);
  Meson etaprime = Meson(958, "\\eta'", 957.78, -94.7, MesonType::pseudoscalar, Charge::neutral);
  Meson etac = Meson(441, "\\eta_c", 2983.6, 237, MesonType::pseudoscalar, Charge::neutral);

  Meson rho = Meson(213, "\\rho", 775.11, 162000, MesonType::vector, Charge::charged);
  Meson Dstar = Meson(413, "D^{\\ast+}", 2010.26, 535000, MesonType::vector, Charge::charged);
  Meson Dstars = Meson(431, "D^{\\ast+}_s", 2112.1, 650000, MesonType::vector, Charge::charged);


  Double_t weinberg = 0.2223;
  Meson rho0 = Meson(113, "\\rho^0", 775.26, 162000, MesonType::vector, Charge::neutral);
  rho0.insertValue("kh", 1.-2.*weinberg);
  Meson omega = Meson(223, "\\omega", 782.65, 153000, MesonType::vector, Charge::neutral);
  omega.insertValue("kh", 4./3.*weinberg);
  Meson phi = Meson(333, "\\phi", 1019.461, 234000, MesonType::vector, Charge::neutral);
  phi.insertValue("kh", (4./3.*weinberg)-1.);
  Meson jpsi = Meson(443, "J/\\Psi", 3096.916, 1290000, MesonType::vector, Charge::neutral);
  jpsi.insertValue("kh", 1.-8./3.*weinberg);

  std::vector<Meson> mesons = {pi, K, D, Ds, B, Bc, pi0, eta, etaprime, etac, rho, Dstar, Dstars, rho0, omega, phi, jpsi};

  while (1)
    {
      int c;
      static struct option long_options[] =
        {
          /* These options set a flag. */
          {"loglevel",     required_argument, 0, 'l'},
          {"generations",  required_argument, 0, 'g'},
          {"angle",        required_argument, 0, 'u'},
          {"mass",         required_argument, 0, 'm'},
          {0, 0, 0, 0}
        };

      /* getopt_long stores the option index here. */
      int option_index = 0;

      c = getopt_long (argc, argv, "v:l:u:m:",
                       long_options, &option_index);

      /* Detect the end of the options. */
      if (c == -1)
        break;

      switch (c)
        {
        case 0:
          /* If this option set a flag, do nothing else now. */
          if (long_options[option_index].flag != 0)
            break;
          printf ("option %s", long_options[option_index].name);
          if (optarg)
            printf (" with arg %s", optarg);
          printf ("\n");
          break;

        case 'l':
          {
            TString x = TString(optarg);
            if(x.CompareTo("debug")==0) {
              gLOGLEVEL = Level::DEBUG;

            } else if (x.CompareTo("warning")==0) {
              gLOGLEVEL = Level::WARNING;

            } else if (x.CompareTo("error")==0) {
              gLOGLEVEL = Level::ERROR;

            } else {
              gLOGLEVEL = Level::INFO;
            }
          }

          break;

        case 'g':
          {
            TString x = TString(optarg);
            TObjArray *tx = x.Tokenize(",");

            for (Int_t i = 0; i < tx->GetEntries(); i++) {
              TString iter_str = ((TObjString*)tx->At(i))->CopyString();
              for(auto l : all_leptons) {
                if(iter_str.CompareTo(l.getName())==0) mixes_with.emplace_back(l);
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
          abort ();
        }
    }

  /* Print any remaining command line arguments (not options). */
  if (optind < argc)
    {
      printf ("non-option ARGV-elements: ");
      while (optind < argc)
        printf ("%s ", argv[optind++]);
      putchar ('\n');
    }


  for(auto l : mixes_with) {
    LOG_DEBUG("HNL mixes with " << l.getName());
  }
  // this is the HNL we configure via command line options. If it mixes with multiple generations it can only
  // assume mixing angle ratios 1:1:..:1
  HNL N = HNL("HNL", HNLmass, angle, mixes_with);

  Double_t Gamma = N.getTotalWidth(cfg, all_leptons, mesons);
  LOG_INFO("mass=" << N.getMass()/1000. << " GeV, " <<
           "ctau=" << gamma2ctau(cfg, Gamma) << " mm");

  auto ch = N.getDecayChannels();
  for (auto it=ch.begin(); it!=ch.end(); ++it) {
    if(it->second > 0)
      std::cout << "\\Gamma(" << it->first << ") = \\SI{" << it->second << "}{}\\\\\n";
  }


  plot_br_low(cfg, all_leptons, mesons, N, "BR_low.pdf", 50, 1000, 1);

  std::vector<Meson> vector_charged;
  for(auto m : mesons) {
    if(m.getMesonType() != MesonType::vector) continue;
    if(m.getCharge() != Charge::charged) continue;
    vector_charged.emplace_back(m);
  }

  TF1* f = new TF1("#Delta_{QCD}", qcd_coupling, 1, 100, 0);
  Double_t qcd_corr = qcd_correction(f->Eval(1776/1000.));
  LOG_INFO("QCD correction at tau mass: " << qcd_corr);
  LOG_INFO("QCD correction at 20 GeV: " << qcd_correction(f->Eval(20000/1000.)));
  LOG_INFO("QCD coupling at tau: " << qcd_coupling(1776./1000.));
  LOG_INFO("QCD coupling at Z mass: " << qcd_coupling(91.2));
  delete f;

  plot_I();
  plot_br(cfg, all_leptons, mesons, N, "BR.pdf", 1000, 5000);

  plot_qcd_correction("qcd_corr.pdf");
  plot_qcd_coupling("qcd_coupl.pdf");

  bool isMajorana = N.isMajorana();
  N.setMajorana(false); // disable majorana here because we want to compare the pw for one charge configuration only
  plot_meson_pw(cfg, mu, vector_charged, N, "mesons_vector_charged.pdf", 500, 5000);
  N.setMajorana(isMajorana); // restore majorana value

  LOG_INFO("Searching now for U2 for 0.1 mm ctau... This might take a while!");
  Double_t foundU2 = ctauToU2(cfg, 0.1, all_leptons, mesons, N);
  LOG_INFO(std::endl << "FOUND U2: " << foundU2);

  // EXAMPLE for further checks and configurations
  // ----------
  // N.setMass(3000);
  // LOG_INFO("pw: " << N.getPartialWidth(cfg, mu, rho));


  return EXIT_SUCCESS;
}
