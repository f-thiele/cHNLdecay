// Copyright (C) 2018 - Fabian A.J. Thiele, <fabian.thiele@cern.ch>

#ifndef   PARTICLE_H
#define   PARTICLE_H

#include "TString.h"

class Particle {
 public:
  Particle() {
    pdgid = 0;
    name = "";
    mass = 0;
  }
  Particle(Int_t p, TString n, Double_t m) {
    pdgid = p;
    name = n;
    mass = m;
  }
  Particle(const Particle &obj) {
    pdgid = obj.getPdgId();
    name  = obj.getName();
    mass  = obj.getMass();
  }

  Double_t getMass() const {
    return mass;
  }
  TString getName() const {
    return name;
  }

  Int_t getPdgId() const {
    return pdgid;
  }

  bool operator==(const Particle& a) const {
    return pdgid == a.getPdgId(); // true if they have the same pdg IDs
  }

 protected:
  Int_t pdgid;
  TString name;
  Double_t mass;
};
#endif
