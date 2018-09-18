// Copyright (C) 2018 - Fabian A.J. Thiele, <fabian.thiele@cern.ch>

#ifndef   LEPTON_H
#define   LEPTON_H

#include "TString.h"
#include "Particle.h"

class Lepton: public Particle {
  public:
    Lepton() : Particle() {}
    Lepton(Int_t p, TString n, Double_t m) : Particle(p, n, m) {}
    Lepton(const Lepton &obj) : Particle(obj) {}
};
#endif
