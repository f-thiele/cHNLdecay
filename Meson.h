// Copyright (C) 2018 - Fabian A.J. Thiele, <fabian.thiele@cern.ch>

#ifndef   MESON_H
#define   MESON_H
#include <map>
#include "TString.h"
#include "Particle.h"

enum class MesonType { Unknown, pseudoscalar, vector };
enum class Charge { Unknown, charged, neutral };

class Meson: public Particle {
public:
  Meson() : Particle() {
    decayConstant = 0;
    vals = {};
    charge = Charge::Unknown;
    type = MesonType::Unknown;
  }
  Meson(Int_t p, TString n, Double_t m, Double_t c, MesonType t, Charge q) : Particle(p, n, m) {
    decayConstant = c;
    vals = {};
    charge = q;
    type = t;
  }
  Meson(const Meson &obj) : Particle(obj) {
    decayConstant = obj.getDecayConstant();
    vals          = obj.getValueMap();
    charge        = obj.getCharge();
    type          = obj.getMesonType();
  }

  MesonType getMesonType() const {
    return type;
  }

  Charge getCharge() const {
    return charge;
  }

  Double_t getDecayConstant() const {
    return decayConstant;
  }

  bool hasValue(TString name) const {
    return vals.count(name)>0;
  }

  Double_t getValue(TString name) const {
    return vals.at(name);
  }

  void insertValue(TString name, Double_t value) {
    vals.insert( std::pair<TString, Double_t>(name, value) );
  }

  std::map<TString, Double_t> getValueMap() const {
          return vals;
  }

private:
  Double_t decayConstant;
  std::map<TString, Double_t> vals;
  MesonType type;
  Charge charge;
};
#endif
