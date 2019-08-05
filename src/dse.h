#ifndef dse_H
#define dse_H

#include "global.h"
#include "utility/utility.h"
#include <array>
#include <string>
#include <vector>

extern parameter Para;

namespace dse {
using namespace std;

struct green {
  momentum *Mom;
  int Tau[2];
  double Weight;
};

struct ver4 {
  int Level;
  ver4 *SubVer[2];
  momentum *Legs[4];
  int ExtTau[4];
  green *Internal[MaxOrder * 2];
  double Weight[MaxOrder * MaxOrder * 4];
};

struct pool {
  momentum Mom[MaxDiagNum * MaxOrder * 2];
  int MomIndex;
  green G[MaxDiagNum * MaxOrder * 2];
  int GIndex;
  ver4 Ver4[MaxDiagNum];
  int Ver4Index;
};

ver4 *Bubble()

} // namespace dse

#endif