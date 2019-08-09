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
  int Channel; // 0: I, 1: T, 2: U, 3: S, 23: IUS, 123: ITUS
  int Side;    // 0: left, 1: right
  int Level;
  int DiagNum;
  int LoopNum;
  int TauIndex;
  int LoopIndex;
  ver4 *LVer;
  ver4 *RVer;
  int LegK[4];
  int InT[2];
  momentum Internal[2];
  vector<double> GWeight;
  vector<vertype> Type; //-1: bare/dynamic coupling,  0: normal vertex,
                        // 1: projected vertex, 2: renormalized vertex
  vector<array<int, 2>> OutT[2];
  vector<double> Weight;
  // double Weight[MaxOrder * MaxOrder * 4];
};

struct pool {
  vector<momentum> Mom;
  vector<ver4> Ver4;
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