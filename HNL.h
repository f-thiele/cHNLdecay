#ifndef   HNL_H
#define   HNL_H

#include <iostream>
#include "TString.h"
#include "Config.h"
#include "Lepton.h"
#include "Meson.h"

class HNL {
 public:
  HNL(TString n, Double_t m, Double_t U2, const Lepton& a) {
    name = n;
    mass = m;
    angle = U2;
    generation = Lepton(a.getName(), a.getMass());
  }

  bool mixesWith(const Lepton& a) const {
    return generation==a;
  }

  Double_t getMass() const {
    return mass;
  }

  Double_t getAngle() const {
    return angle;
  }

  void setMass(Double_t m) {
    mass = m;
  }

  void setAngle(Double_t a) {
    angle = a;
  }

  TString getName() const {
    return name;
  }

  Lepton getGeneration() const {
    return Lepton(generation);
  }

  Double_t getPartialWidth(std::shared_ptr<Config> cfg, Lepton &alpha, Lepton &beta) const;
  Double_t getPartialWidth(std::shared_ptr<Config> cfg, Lepton &alpha, Meson &m) const;

 private:
  TString name;
  Double_t mass;
  Double_t angle;
  Lepton generation;
};
#endif
