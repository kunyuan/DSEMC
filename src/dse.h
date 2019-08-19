#ifndef dse_H
#define dse_H

#include "global.h"
#include "utility/utility.h"
#include "vertex.h"
#include <array>
#include <map>
#include <string>
#include <tuple>
#include <vector>

extern parameter Para;

namespace dse {
using namespace std;

enum caltype { BARE, RG, PARQUET, RENORMALIZED };
enum channel { I = 0, T, U, S };

struct pair;
struct bubble;
struct envelope;

class gMatrix {
public:
  gMatrix() {
    _TauNum = 0;
    _InTL = 0;
  }
  gMatrix(int TauNum, int InTL, momentum *k) {
    _TauNum = TauNum;
    _InTL = InTL;
    K = k;
    _G.resize(TauNum * TauNum);
    for (auto &g : _G)
      g = 0.0;
  }
  double &operator()(int l, int r) {
    return _G[(l - _InTL) * _TauNum + r - _InTL];
  }

  double &operator()(const array<int, 2> &t) {
    return _G[(t[IN] - _InTL) * _TauNum + t[OUT] - _InTL];
  }

  momentum *K;

private:
  int _TauNum;
  int _InTL;
  vector<double> _G;
};

struct ver4 {
  int ID;
  int InTL;
  int LoopNum;
  int TauNum;
  caltype Type;
  bool ReExpandBare;
  bool ReExpandVer4;
  vector<channel> Channel;
  array<momentum *, 4> LegK; // legK index
  vector<array<int, 4>> T;   // external T list

  // bubble diagram
  vector<pair> Pairs; // all two-particle reducible diagrams
  array<gMatrix, 4> G;

  vector<bubble> Bubble;     // bubble diagrams
  vector<envelope> Envelope; // envelop diagrams at order 4

  vector<double> Weight; // size: equal to T.size()
};

struct mapT {
  int LVerT;
  int RVerT;
  // LVer T index and RVer T index to Internal T for G1 and G2
  array<int, 2> G1T;
  array<int, 2> G2T;
  // map LVer T index and RVer T index to merged T index
  int T;
};

struct mapT2 {
  int LVerT;
  int RVerT;
  // LVer T index and RVer T index to Internal T for G1 and G2
  array<array<int, 2>, 4> GT; // four indepdent G
  // map LVer T index and RVer T index to merged T index
  array<int, 3> T; // three channels
};

struct bubble {
  array<ver4, 3> LVer;
  array<ver4, 3> RVer;
  array<double, 3> SymFactor;
  array<gMatrix, 4> G;
  vector<mapT2> Map;
};

struct pair {
  ver4 LVer;
  ver4 RVer;
  channel Chan;
  double SymFactor;
  vector<mapT> Map;
};

class g2Matrix {
public:
  g2Matrix() {
    _InShift = 0;
    _OutShift = 0;
  }
  g2Matrix(int InShift, int OutShift, momentum *k) {
    _InShift = InShift;
    _OutShift = OutShift;
    K = k;
    _G.resize(2 * 2);
    for (auto &g : _G)
      g = 0.0;
    InT = {_InShift, _InShift + 1};
    OutT = {_OutShift, _OutShift + 1};
  }

  double &operator()(int in, int out) {
    return _G[(in - _InShift) * 2 + out - _OutShift];
  }

  double &operator()(const array<int, 2> &t) {
    return _G[(t[IN] - _InShift) * 2 + t[OUT] - _OutShift];
  }

  array<int, 2> InT, OutT; // possible InT and OutT
  momentum *K;             // momentum on G

private:
  int _InShift;
  int _OutShift;
  vector<double> _G;
};

struct mapT4 {
  int LDVerT;
  int RDVerT;
  int LUVerT;
  int RUVerT;
  // LVer T index and RVer T index to Internal T for G1 and G2
  array<array<int, 2>, 9> GT;
  // map LVer T index and RVer T index to merged T index
  array<int, 4> T; // external T for four envelop diagrams
};

struct envelope {
  array<ver4, 10> Ver;
  array<g2Matrix, 9> G;
  vector<mapT4> Map;
  array<double, 4> SymFactor;
};

class verDiag {
public:
  ver4 Build(array<momentum, MaxMomNum> &loopmom, int LoopNum,
             vector<channel> Channel, caltype Type);
  string ToString(const ver4 &Vertex);

private:
  int DiagNum = 0;
  int MomNum = MaxLoopNum;
  array<momentum, MaxMomNum> *LoopMom; // all momentum loop variables

  ver4 Vertex(array<momentum *, 4> LegK, int InTL, int LoopNum, int LoopIndex,
              vector<channel> Channel, caltype Type, int Side);

  ver4 Ver0(ver4 Ver4, int InTL, int Side);
  ver4 ChanI(ver4 Ver4, int InTL, int LoopNum, int LoopIndex, int Side);
  ver4 ChanUST(ver4 Ver4, int InTL, int LoopNum, int LoopIndex, channel Channel,
               int Side);
  momentum *NextMom();
};

bool verTest();

} // namespace dse
#endif