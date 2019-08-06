#ifndef dse_H
#define dse_H

#include "global.h"
#include "utility/utility.h"
#include "vertex.h"
#include <array>
#include <string>
#include <vector>

extern parameter Para;

namespace dse {
using namespace std;

struct green {
  momentum *Mom;
  int Type;
  int Tau[2];
  double Weight;
};

const int MAXDIAG = MaxOrder * MaxOrder * 4;
enum vertype { BARE, DYNAMIC, NORMAL, PROJECTED, RENORMALIZED };

struct ver4 {
  int Level;
  int DiagNum;
  vertype Type[MAXDIAG]; //-2: bare coupling, -1: dynamic coupling, 0: normal
                         // vertex, 1:
                         // projected vertex, 2: renormalized vertex
  int LoopNum;
  bool IsProjection;
  ver4 *SubVer[2];
  int LegK[4];
  int InT[2];
  int OutT[MAXDIAG][2];
  green *Internal[MaxOrder * 2];
  // double Weight[MaxOrder * MaxOrder * 4];
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

class verDiag {
public:
  // verDiag(array<momentum, MaxLoopNum> &loopmom, array<double, MaxTauNum>
  // &tau)
  //     : LoopMom(loopmom), Tau(tau){};
  void Build(int LoopNum);
  pool Pool;

private:
  // array<momentum, MaxLoopNum> &LoopMom; // all momentum loop variables
  // array<double, MaxTauNum> &Tau;        // all tau variables

  // verQTheta VerWeight;

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