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

#ifndef QUARK_H
#define QUARK_H
#include "Particle.h"
#include "TString.h"
#include "Meson.h"
#include <map>

class Quark : public Particle {
public:
  Quark() : Particle() {
    q = Quark_Type::Unknown;
  }
  Quark(Int_t p, Double_t m, Quark_Type type = Quark_Type::Unknown) : Particle(p, m) {
    q = type;
  }
  Quark(const Quark &obj) : Particle(obj) {
    q = obj.getQuarkType();

  }

  Quark &operator=(const Quark &obj) {
    Particle::operator=(obj);

    q = obj.getQuarkType();

    return *this;
  }

  Quark_Type getQuarkType() const {
    return q;
  }

private:
  Quark_Type q;
};
#endif
