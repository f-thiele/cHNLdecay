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

#ifndef MESON_H
#define MESON_H
#include "Particle.h"
#include "TString.h"
#include <map>

enum class MesonType { Unknown, pseudoscalar, vector };
enum class Charge { Unknown, charged, neutral };

class Meson : public Particle {
public:
  Meson() : Particle() {
    decayConstant = 0;
    charge = Charge::Unknown;
    type = MesonType::Unknown;
  }
  Meson(Int_t p, Double_t m, Double_t c, MesonType t, Charge q)
      : Particle(p, m) {
    decayConstant = c;
    charge = q;
    type = t;
  }
  Meson(const Meson &obj) : Particle(obj) {
    decayConstant = obj.getDecayConstant();
    vals = obj.getValueMap();
    charge = obj.getCharge();
    type = obj.getMesonType();
  }

  Meson& operator=(const Meson &obj) {
    Particle::operator=(obj);

    decayConstant = obj.getDecayConstant();
    vals = obj.getValueMap();
    charge = obj.getCharge();
    type = obj.getMesonType();

    return *this;
  }

  MesonType getMesonType() const { return type; }

  Charge getCharge() const { return charge; }

  Double_t getDecayConstant() const { return decayConstant; }

  bool hasValue(TString name) const { return vals.count(name) > 0; }

  Double_t getValue(TString name) const { return vals.at(name); }

  void insertValue(TString name, Double_t value) {
    vals.insert(std::pair<TString, Double_t>(name, value));
  }

  std::map<TString, Double_t> getValueMap() const { return vals; }

private:
  Double_t decayConstant;
  std::map<TString, Double_t> vals;
  MesonType type;
  Charge charge;
};
#endif
