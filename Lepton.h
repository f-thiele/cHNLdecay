// Copyright (C) 2018 - Fabian A.J. Thiele, <fabian.thiele@cern.ch>

#ifndef LEPTON_H
#define LEPTON_H

#include "Particle.h"
#include "TString.h"

class Lepton : public Particle {
public:
  Lepton() : Particle() {}
  Lepton(Int_t p, Double_t m) : Particle(p, m) {}
  Lepton(const Lepton &obj) : Particle(obj) {}
};
#endif
