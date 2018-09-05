#ifndef   LEPTON_H
#define   LEPTON_H

#include "TString.h"

class Lepton {
 public:
  Lepton() {
    name = "";
    mass = 0;
  }
  Lepton(TString n, Double_t m) {
    name = n;
    mass = m;
  }
  Lepton(const Lepton &obj) {
    name = obj.getName();
    mass = obj.getMass();
  }

  Double_t getMass() const {
    return mass;
  }
  TString getName() const {
    return name;
  }

  bool operator==(const Lepton& a) const {
    return getMass() == a.getMass(); // true if they have same masses false otherwise
  }

 private:
  TString name;
  Double_t mass;

};
#endif
