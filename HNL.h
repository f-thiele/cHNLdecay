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
    generation.emplace_back(Lepton(a.getName(), a.getMass()));
  }
  HNL(TString n, Double_t m, Double_t U2, const std::vector<Lepton> &a) {
    name = n;
    mass = m;
    angle = U2;
    for(auto l : a)
      generation.emplace_back(Lepton(l.getName(), l.getMass()));
  }

  bool mixesWith(const Lepton& a) const {
    bool mixes = false;
    for(auto g : generation) {
      if(g==a) mixes = true;
    }
    return mixes;
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

  std::vector<Lepton> getGeneration() const {
    return generation;
  }

  Double_t getPartialWidth(std::shared_ptr<Config> cfg, const Lepton &alpha, const Lepton &beta, bool invisible=true) const;
  Double_t getPartialWidthInv(std::shared_ptr<Config> cfg, const Lepton &alpha, const Lepton &beta) const;
  Double_t getPartialWidth(std::shared_ptr<Config> cfg, const Lepton &alpha, const Meson &m) const;

 private:
  TString name;
  Double_t mass;
  Double_t angle;
  std::vector<Lepton> generation;
};
#endif
