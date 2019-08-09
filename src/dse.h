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
enum channel { IDIR, IEX, T, U, S };

struct ver4 {
  // int Channel; // 0: I, 1: T, 2: U, 3: S, 23: IUS, 123: ITUS
  int Side; // 0: left, 1: right
  int Level;
  int DiagNum;
  int LoopNum;
  int TauIndex;
  int LoopIndex;
  vector<ver4 *> LVer; // left vertex has three channels
  vector<ver4 *> RVer; // right vertex has four channels
  int LegK[4];
  int InT[2];
  momentum Internal[2];
  vector<double> GWeight;

  vector<channel> Channel;
  vector<vertype> Type;
  vector<array<int, 2>> OutT;
  array<double, 4 * MaxOrder * MaxOrder> Weight;
  // double Weight[MaxOrder * MaxOrder * 4];
};

// struct pool {
//   vector<momentum> Mom;
//   vector<ver4> Ver4;
// };

class verDiag {
public:
  // verDiag(array<momentum, MaxLoopNum> &loopmom, array<double, MaxTauNum>
  // &tau)
  //     : LoopMom(loopmom), Tau(tau){};
  void Build(int LoopNum);
  vector<ver4> VerPool;
  // pool Pool;

private:
  // array<momentum, MaxLoopNum> &LoopMom; // all momentum loop variables
  // array<double, MaxTauNum> &Tau;        // all tau variables

  // verQTheta VerWeight;

  ver4 *Ver4AtOrder(
      array<int, 4> LegK, array<int, 2> InLegT,
      int Type // -1: DSE diagrams, 0: normal digrams, 1: renormalzied diagrams
  );

  ver4 *Ver0(array<int, 2> InT, int LoopNum, vertype Type);
  ver4 *ChanI(array<int, 2> InT, int LoopNum, vertype Type);

  ver4 *Bubble(array<int, 2> InT, int LoopNum, vector<channel> Channel,
               vertype Type);
};

} // namespace dse
#endif