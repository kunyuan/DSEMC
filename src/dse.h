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
  int DiagNum;
  bool IsProjection;
  ver4 *SubVer[2];
  momentum *LegK[4];
  int LegT[4];
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

typedef array<bool, 3> channel;

class verDiagram {
public:
  void Build(int LoopNum);
  pool Pool;

private:
  ver4 *Ver4AtOrder(
      array<int, 4> LegK, array<int, 2> InLegT,
      int Type // -1: DSE diagrams, 0: normal digrams, 1: renormalzied diagrams
  );

  ver4 *Ver0(array<int, 4> LegK, array<int, 2> InLegT, bool IsBare = false);

  ver4 *Bubble(array<int, 4> LegK, array<int, 2> InLegT,
               array<int, 2> SubVerType, array<int, 2> SubVerLoopNum,
               array<channel, 2> SubVerChannel, int TauIndex, int LoopIndex,
               bool IsProjection);

  channel ALL = {true, true, true};
  channel US[3] = {false, true, true};
  channel UT[3] = {true, true, false};
  channel ST[3] = {true, false, true};
  channel T[3] = {true, false, false};
  channel U[3] = {false, true, false};
  channel S[3] = {false, false, true};
};

} // namespace dse
#endif