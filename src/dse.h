#ifndef dse_H
#define dse_H

#include "global.h"
#include "utility/utility.h"
#include "vertex.h"
#include <array>
#include <string>
#include <tuple>
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
enum channel { I, T, U, S };

struct pair {
  int LVer;
  int RVer;
  channel Chan;
  vertype Type;
  vector<vector<int>> Map;
};

struct ver4 {
  // int Channel; // 0: I, 1: T, 2: U, 3: S, 23: IUS, 123: ITUS
  int ID;
  int Level;
  int LoopNum;
  int TauIndex;
  int LoopIndex;
  vertype Type;
  vector<channel> Channel;
  momentum Internal2, ExQ;
  momentum VerLInL, VerLInR, VerLDiTran, VerRInL, VerRInR, VerRDiTran;
  vector<double> GL2R, GR2L;

  // vector<channel> Channel;
  vector<ver4> SubVer; // subver list
  vector<pair> Pairs;
  vector<array<int, 4>> T;
  vector<double> Weight;
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
  ver4 Build(int LoopNum, vector<channel> Channel, vertype Type);
  vector<ver4> VerPool;
  string ToString(const ver4 &Vertex);
  // pool Pool;

private:
  // array<momentum, MaxLoopNum> &LoopMom; // all momentum loop variables
  // array<double, MaxTauNum> &Tau;        // all tau variables

  // verQTheta VerWeight;
  int DiagNum = 0;

  ver4 Ver0(int InTL, vertype Type);
  ver4 ChanI(int InTL, int LoopNum, vertype Type);

  ver4 Bubble(int InTL, int LoopNum, vector<channel> Channel, vertype Type);
};

bool verTest();

} // namespace dse
#endif