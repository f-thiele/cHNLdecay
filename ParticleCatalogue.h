// Copyright (C) 2018 - Fabian A.J. Thiele, <fabian.thiele@cern.ch>

#ifndef PARTICLECATALOGUE_H
#define PARTICLECATALOGUE_H
#include "Lepton.h"
#include "Meson.h"
#include <map>

class ParticleCatalogue {
public:
  ParticleCatalogue() {
    Double_t muonMass = 105.6583715;      // MeV
    Double_t electronMass = 0.5109989461; // MeV
    Double_t tauMass = 1776.82;           // MeV

    Lepton el = Lepton(11, electronMass);
    addTo(el);
    Lepton mu = Lepton(13, muonMass);
    addTo(mu);
    Lepton tau = Lepton(15, tauMass);
    addTo(tau);

    Meson pi =
        Meson(211, 139.57018, 130.2, MesonType::pseudoscalar, Charge::charged);
    addTo(pi);
    Meson K =
        Meson(321, 493.677, 155.6, MesonType::pseudoscalar, Charge::charged);
    addTo(K);
    Meson D =
        Meson(411, 1869.62, 212, MesonType::pseudoscalar, Charge::charged);
    addTo(D);
    Meson Ds =
        Meson(431, 1968.47, 249, MesonType::pseudoscalar, Charge::charged);
    addTo(Ds);
    Meson B =
        Meson(521, 5279.29, 187, MesonType::pseudoscalar, Charge::charged);
    addTo(B);
    Meson Bc =
        Meson(541, 6275.1, 434, MesonType::pseudoscalar, Charge::charged);
    addTo(Bc);

    Meson pi0 =
        Meson(111, 134.9766, 130.2, MesonType::pseudoscalar, Charge::neutral);
    addTo(pi0);
    Meson eta =
        Meson(221, 547.862, 81.7, MesonType::pseudoscalar, Charge::neutral);
    addTo(eta);
    Meson etaprime =
        Meson(331, 957.78, -94.7, MesonType::pseudoscalar, Charge::neutral);
    addTo(etaprime);
    Meson etac =
        Meson(441, 2983.6, 237, MesonType::pseudoscalar, Charge::neutral);
    addTo(etac);

    Meson rho = Meson(213, 775.11, 162000, MesonType::vector, Charge::charged);
    addTo(rho);
    Meson Dstar =
        Meson(413, 2010.26, 535000, MesonType::vector, Charge::charged);
    addTo(Dstar);
    Meson Dstars =
        Meson(431, 2112.1, 650000, MesonType::vector, Charge::charged);
    addTo(Dstars);

    Double_t weinberg = 0.2223;
    Meson rho0 = Meson(113, 775.26, 162000, MesonType::vector, Charge::neutral);
    rho0.insertValue("kh", 1. - 2. * weinberg);
    addTo(rho0);
    Meson omega =
        Meson(223, 782.65, 153000, MesonType::vector, Charge::neutral);
    omega.insertValue("kh", 4. / 3. * weinberg);
    addTo(omega);
    Meson phi =
        Meson(333, 1019.461, 234000, MesonType::vector, Charge::neutral);
    phi.insertValue("kh", (4. / 3. * weinberg) - 1.);
    addTo(phi);
    Meson jpsi =
        Meson(443, 3096.916, 1290000, MesonType::vector, Charge::neutral);
    jpsi.insertValue("kh", 1. - 8. / 3. * weinberg);
    addTo(jpsi);
  }

  Meson getMeson(Int_t p) const { return mesons.at(p); }

  std::vector<Meson> getAllMesons() const {
    std::vector<Meson> ret;
    for (auto m : mesons) {
      ret.emplace_back(m.second);
    }
    return ret;
  }

  Lepton getLepton(Int_t p) const { return leptons.at(p); }

  std::vector<Lepton> getAllLeptons() const {
    std::vector<Lepton> ret;
    for (auto l : leptons) {
      ret.emplace_back(l.second);
    }
    return ret;
  }

  void addTo(const Meson &p) {
    mesons.insert(std::make_pair<Int_t, Meson>(p.getPdgId(), Meson(p)));
  }

  void addTo(const Lepton &p) {
    leptons.insert(std::make_pair<Int_t, Lepton>(p.getPdgId(), Lepton(p)));
  }

private:
  std::map<Int_t, Meson> mesons;
  std::map<Int_t, Lepton> leptons;
};
#endif
