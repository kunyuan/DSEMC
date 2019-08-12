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

enum vertype { BARE, DYNAMIC, NORMAL, PROJECTED, RENORMALIZED };
enum channel { I = 0, T, U, S };

struct pair;
class gMatrix {
public:
  gMatrix() {
    _TauNum = 0;
    _InTL = 0;
  }
  gMatrix(int TauNum, int InTL) {
    _TauNum = TauNum;
    _InTL = InTL;
    _G.resize(TauNum * TauNum);
    for (auto &g : _G)
      g = 0.0;
  }
  double &operator()(int l, int r) {
    return _G[(l - _InTL) * _TauNum + r - _InTL];
  }

private:
  int _TauNum;
  int _InTL;
  vector<double> _G;
};

template <typename T> class map {
  // IMPORTANT: please make sure all elements are initialized!
public:
  map(int l, int r) {
    LNum = l;
    RNum = r;
    _Index.resize(l * r);
  }
  T &operator()(int l, int r) { return _Index[l * RNum + r]; }
  // int Set(int l, int r, int index) { _Index[l * RNum + r] = index; }

private:
  int LNum;
  int RNum;
  vector<T> _Index;
};

struct ver4 {
  // int Channel; // 0: I, 1: T, 2: U, 3: S, 23: IUS, 123: ITUS
  int ID;
  int LoopNum;
  int InTL;
  int TauNum;
  vertype Type;
  vector<channel> Channel;
  array<int, 4> LegK; // legK index
  vector<array<int, 4>> T;

  // bubble diagram
  vector<pair> Pairs;
  int K1; // internal K index for t, u, s
  gMatrix G1;
  array<int, 4> K2; // internal K2 for four different channel
  array<gMatrix, 4> G2;

  // TODO: envelope diagram

  vector<double> Weight; // size: equal to T.size()
};

struct pair {
  ver4 LVer;
  ver4 RVer;
  // map LVer T index and RVer T index to Internal T for G1 and G2
  map<array<int, 2>> IntT1;
  map<array<int, 2>> IntT2;
  // map LVer T index and RVer T index to merged T index
  map<int> Map;
  channel Chan;
  double SymFactor;
};

class verDiag {
public:
  ver4 Build(int LoopNum, vector<channel> Channel, vertype Type);
  string ToString(const ver4 &Vertex);

private:
  int DiagNum = 0;
  int MomNum = MaxLoopNum;

  ver4 Vertex(array<int, 4> LegK, int InTL, int LoopNum, int LoopIndex,
              vector<channel> Channel, vertype Type, int Side);

  ver4 Ver0(ver4 Ver4, int InTL, vertype Type, int Side);
  ver4 ChanI(ver4 Ver4, int InTL, int LoopNum, int LoopIndex, vertype Type,
             int Side);
  ver4 ChanUST(ver4 Ver4, int InTL, int LoopNum, int LoopIndex, channel Channel,
               vertype Type, int Side);
  int NewMom();
};

bool verTest();

} // namespace dse
#endif