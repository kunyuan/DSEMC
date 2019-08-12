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
enum channel { I = 0, T, U, S };

class map {
public:
  map(int l, int r) {
    LNum = l;
    RNum = r;
    _Index.resize(l * r);
  }
  const int Get(int l, int r) { return _Index[l * RNum + r]; }
  int Set(int l, int r, int index) { _Index[l * RNum + r] = index; }

private:
  int LNum;
  int RNum;
  vector<int> _Index;
};

struct pair;

struct ver4 {
  // int Channel; // 0: I, 1: T, 2: U, 3: S, 23: IUS, 123: ITUS
  int ID;
  int Level;
  int LoopNum;
  int TauNum;
  vertype Type;
  vector<channel> Channel;

  vector<pair> Pairs;
  vector<array<int, 4>> T;

  int K1 = -1; // internal K index for t, u, s
  array<int, 4> K2;
  array<int, 4> LegK; // legK index

  array<int, MaxTauNum * MaxTauNum> G1; // size: TauNum*TauNum
  array<array<int, MaxTauNum * MaxTauNum>, 4> G2;
  vector<double> Weight; // size: equal to T.size()
};

struct pair {
  ver4 LVer;
  ver4 RVer;
  vector<array<int, 4>> InterT; // G1Start, G1End, G2Start, G2End
  map Map; // map LVer T index and RVer T index to merged T index
  channel Chan;
  int SymFactor;
};

// struct pool {
//   vector<momentum> Mom;
//   vector<ver4> Ver4;
// };

class verDiag {
public:
  ver4 Build(int LoopNum, vector<channel> Channel, vertype Type);
  string ToString(const ver4 &Vertex);

private:
  // array<momentum, MaxLoopNum> &LoopMom; // all momentum loop variables
  // array<double, MaxTauNum> &Tau;        // all tau variables

  // verQTheta VerWeight;
  int DiagNum = 0;
  int MomNum = MaxLoopNum;

  ver4 Vertex(array<int, 4> LegK, int InTL, int LoopNum, int LoopIndex,
              vector<channel> Channel, vertype Type, int Side);

  ver4 Ver0(ver4 Ver4, int InTL, vertype Type, int Side);
  ver4 ChanI(ver4 Ver4, int InTL, int LoopNum, int LoopIndex, vertype Type,
             int Side);
  ver4 ChanUST(ver4 Ver4, int InTL, int LoopNum, int LoopIndex, channel Channel,
               vertype Type, int Side);
  int NextMom();
};

bool verTest();

} // namespace dse
#endif