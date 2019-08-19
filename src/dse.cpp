#include "dse.h"
#include "global.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/vector.h"
#include <array>
#include <iostream>
#include <string>

using namespace dse;
using namespace std;

double Sym(channel chan) {
  if (chan == T)
    return -1.0;
  else if (chan == U)
    return 1.0;
  else if (chan == S)
    return 0.5;
  else
    return 1.0;
}

int AddToTList(vector<array<int, 4>> &TList, const array<int, 4> T) {
  // find the T array in the list, if failed, create a new array
  for (int i = 0; i < TList.size(); i++) {
    auto t = TList[i];
    ASSERT_ALLWAYS(t[INL] == T[INL],
                   "left Tin must be the same for all subvertex!");
    if (t[OUTL] == T[OUTL] && t[INR] == T[INR] && t[OUTR] == T[OUTR])
      return i;
  }
  TList.push_back(T);
  return TList.size() - 1;
}

momentum *verDiag::NextMom() {
  MomNum += 1;
  ASSERT_ALLWAYS(MomNum < MaxMomNum, "Too many momentum variables! " << MomNum);
  return &(*LoopMom)[MomNum - 1];
}

ver4 verDiag::Build(array<momentum, MaxMomNum> &loopMom, int LoopNum,
                    vector<channel> Channel, caltype Type) {
  ASSERT_ALLWAYS(LoopNum > 0, "LoopNum must be larger than zero!");
  DiagNum = 0;
  MomNum = MaxLoopNum;
  LoopMom = &loopMom;
  array<momentum *, 4> LegK = {&(*LoopMom)[1], NextMom(), &(*LoopMom)[2],
                               NextMom()};
  return Vertex(LegK, 0, LoopNum, 3, Channel, Type, LEFT);
}

ver4 verDiag::Vertex(array<momentum *, 4> LegK, int InTL, int LoopNum,
                     int LoopIndex, vector<channel> Channel, caltype Type,
                     int Side) {
  ver4 Ver4;
  Ver4.ID = DiagNum;
  DiagNum++;
  Ver4.LoopNum = LoopNum;
  Ver4.TauNum = 2 * (LoopNum + 1);
  Ver4.Type = Type;
  Ver4.LegK = LegK;

  if (Type == caltype::BARE) {
    Ver4.ReExpandBare = false;
    Ver4.ReExpandVer4 = false;
  } else if (Type == caltype::RG && Type == caltype::RENORMALIZED) {
    // In RG and renormalization calculation, the projected vertex will be
    // measured
    Ver4.ReExpandBare = true;
    Ver4.ReExpandVer4 = true;
  } else if (Type == caltype::PARQUET) {
    Ver4.ReExpandBare = false;
    Ver4.ReExpandVer4 = true;
  }

  if (LoopNum == 0) {
    ASSERT_ALLWAYS(Channel[0] == I, "Only I channel has zero loop vertex!");
    Ver4 = Ver0(Ver4, InTL, Side);
  } else {
    vector<channel> UST;
    for (auto &chan : Channel) {
      if (chan == I)
        Ver4 = ChanI(Ver4, InTL, LoopNum, LoopIndex, Side);
      else
        UST.push_back(chan);
    }
    Ver4 = ChanUST(Ver4, UST, InTL, LoopNum, LoopIndex, Side);
  }

  Ver4.Weight.resize(Ver4.T.size());
  for (auto &d : Ver4.Weight)
    d = 0.0;

  return Ver4;
}

ver4 verDiag::Ver0(ver4 Ver4, int InTL, int Side) {
  ////////////// bare interaction ///////////
  if (Side == LEFT)
    // Side==left, then make sure INL Tau are the last TauIndex
    Ver4.T.push_back({InTL, InTL, InTL, InTL});
  else
    // Side==right, then make sure INR Tau are the last TauIndex
    Ver4.T.push_back({InTL + 1, InTL + 1, InTL + 1, InTL + 1});

  if (Ver4.Type != caltype::BARE) {
    //////////// dressed interaction ///////////
    Ver4.T.push_back({InTL, InTL, InTL + 1, InTL + 1});
    Ver4.T.push_back({InTL, InTL + 1, InTL + 1, InTL});
  }
  return Ver4;
}

vector<mapT2> CreateMapT2(ver4 &Ver4, ver4 LVer, ver4 RVer) {
  ///////////   External and Internal Tau  ////////////////
  vector<mapT2> Map;
  array<array<int, 2>, 4> GT;
  array<array<int, 4>, 3> LegT;
  array<int, 3> Tidx;

  for (int lt = 0; lt < LVer.T.size(); ++lt)
    for (int rt = 0; rt < RVer.T.size(); ++rt) {

      auto &LvT = LVer.T[lt];
      auto &RvT = RVer.T[rt];

      GT[0] = {LvT[OUTR], RvT[INL]};
      GT[1] = {RvT[OUTL], LvT[INR]};
      GT[2] = {RvT[OUTL], LvT[INR]};
      GT[3] = {LvT[OUTL], RvT[INL]};

      LegT[0] = {LvT[INL], LvT[OUTL], RvT[INR], RvT[OUTR]};
      LegT[1] = {LvT[INL], RvT[OUTR], RvT[INR], LvT[OUTL]};
      LegT[2] = {LvT[INL], RvT[OUTL], LvT[INR], RvT[OUTR]};

      // add T array into the T pool of the vertex
      for (int i = 0; i < 3; i++)
        Tidx[i] = AddToTList(Ver4.T, LegT[i]);
      Map.push_back(mapT2{lt, rt, Tidx, GT});
    }
  return Map;
}

ver4 verDiag::ChanUST(ver4 Ver4, vector<channel> Channel, int InTL, int LoopNum,
                      int LoopIndex, int Side, bool IsProjected) {
  bubble Bubble;
  Bubble.IsProjected = IsProjected;
  Bubble.Channel = Channel;
  auto &LegK = Ver4.LegK;
  caltype Type = Ver4.Type;
  if (IsProjected == false)
    Bubble.LegK = LegK;

  auto &G = Bubble.G;

  G[0] = gMatrix(Ver4.TauNum, InTL, &(*LoopMom)[LoopIndex]);
  G[1] = gMatrix(Ver4.TauNum, InTL, NextMom());
  G[2] = gMatrix(Ver4.TauNum, InTL, NextMom());
  G[3] = gMatrix(Ver4.TauNum, InTL, NextMom());

  for (int ol = 0; ol < LoopNum; ol++) {
    pair Pair;

    // left and right vertex external LegK
    array<momentum *, 4> LLegK[4], RLegK[4];

    ////////////////// T channel ////////////////////////////
    LLegK[T] = {LegK[INL], LegK[OUTL], G[T].K, G[0].K};
    RLegK[T] = {G[0].K, G[T].K, LegK[INR], LegK[OUTR]};

    ////////////////// U channel ////////////////////////////
    LLegK[U] = {LegK[INL], LegK[OUTR], G[U].K, G[0].K};
    RLegK[U] = {G[0].K, G[U].K, LegK[INR], LegK[OUTL]};

    ////////////////// S channel ////////////////////////////
    LLegK[S] = {LegK[INL], G[S].K, LegK[INR], G[0].K};
    RLegK[S] = {G[0].K, LegK[OUTL], G[S].K, LegK[OUTR]};

    for (auto &c : Bubble.Channel) {
      ////////////////////   Right SubVer  ///////////////////
      int oR = LoopNum - 1 - ol;
      int RInTL = InTL + 2 * (ol + 1);
      int Rlopidx = LoopIndex + 1 + ol;

      if (c == U || c == T) {
        Pair.LVer[c] =
            Vertex(LLegK[c], InTL, ol, LoopIndex + 1, {I, U, S}, Type, LEFT);
        Pair.RVer[c] =
            Vertex(RLegK[c], RInTL, oR, Rlopidx, {I, U, S, T}, Type, RIGHT);
      } else if (c == S) {
        Pair.LVer[c] =
            Vertex(LLegK[c], InTL, ol, LoopIndex + 1, {I, U, T}, Type, LEFT);
        Pair.RVer[c] =
            Vertex(RLegK[c], RInTL, oR, Rlopidx, {I, U, S, T}, Type, RIGHT);
      }
    }

    Pair.Map = CreateMapT2(Ver4, Pair.LVer[0], Pair.RVer[0]);
    Pair.SymFactor = {-1.0, 1.0, -0.5};
    Bubble.Pair.push_back(Pair);
  }

  Ver4.Bubble.push_back(Bubble);
  return Ver4;
}

vector<mapT4> CreateMapT4(ver4 &Ver4, ver4 LDVer, ver4 LUVer, ver4 RDVer,
                          ver4 RUVer) {
  vector<mapT4> Map;
  array<array<int, 2>, 9> GT; // G Tau pair
  array<int, 4> LegT[4], Tidx;

  for (int ldt = 0; ldt < LDVer.T.size(); ldt++)
    for (int lut = 0; lut < LUVer.T.size(); lut++)
      for (int rdt = 0; rdt < RDVer.T.size(); rdt++)
        for (int rut = 0; rut < RUVer.T.size(); rut++) {
          auto &ldT = LDVer.T[ldt];
          auto &luT = LUVer.T[lut];
          auto &rdT = RDVer.T[rdt];
          auto &ruT = RUVer.T[rut];

          // Tau Index for all possible internal G
          GT[0] = {ldT[OUTR], rdT[INL]};
          GT[1] = {ldT[OUTL], luT[INL]};
          GT[2] = {rdT[OUTL], luT[INR]};
          GT[3] = {ruT[OUTL], ldT[INR]};
          GT[4] = {luT[OUTR], ruT[INL]};
          GT[5] = {rdT[OUTR], ruT[INR]};
          GT[6] = {luT[OUTR], ruT[INL]};
          GT[7] = {ruT[OUTR], rdT[INR]};
          GT[8] = {ruT[OUTR], rdT[INR]};

          // external T for four envelope diagram
          // INL, OUTL, INR, OUTR
          LegT[0] = {ldT[INL], luT[OUTL], rdT[INR], ruT[OUTR]};
          LegT[1] = {ldT[INL], ruT[OUTR], rdT[INR], luT[OUTL]};
          LegT[2] = {ldT[INL], luT[OUTL], ruT[INR], rdT[OUTR]};
          LegT[3] = {ldT[INL], rdT[OUTR], ruT[INR], luT[OUTL]};

          for (int i = 0; i < 4; i++)
            Tidx[i] = AddToTList(Ver4.T, LegT[i]);

          Map.push_back(mapT4{ldt, lut, rdt, rut, Tidx, GT});
        }
  return Map;
}

ver4 verDiag::ChanI(ver4 Ver4, int InTL, int LoopNum, int LoopIndex, int Side,
                    bool IsProjected) {

  if (LoopNum != 3)
    return Ver4;

  envelope Env;
  Env.IsProjected = IsProjected;
  if (IsProjected == false)
    Env.LegK = Ver4.LegK;
  auto &G = Env.G;

  int LDInTL = InTL;
  int LUInTL = InTL + 2;
  int RDInTL = InTL + 4;
  int RUInTL = InTL + 6;

  /////// Initialize G Tau and K Table  /////////
  G[0] = g2Matrix(LDInTL, RDInTL, &(*LoopMom)[LoopIndex]);
  G[1] = g2Matrix(LDInTL, LUInTL, &(*LoopMom)[LoopIndex + 1]);
  G[2] = g2Matrix(RDInTL, LUInTL, &(*LoopMom)[LoopIndex + 2]);
  G[3] = g2Matrix(RUInTL, LDInTL, NextMom());
  G[4] = g2Matrix(LUInTL, RUInTL, NextMom());
  G[5] = g2Matrix(RDInTL, RUInTL, NextMom());
  G[6] = g2Matrix(LUInTL, RUInTL, NextMom());
  G[7] = g2Matrix(RUInTL, RDInTL, NextMom());
  G[8] = g2Matrix(RUInTL, RDInTL, NextMom());

  momentum *InL = Ver4.LegK[INL];
  momentum *OutL = Ver4.LegK[OUTL];
  momentum *InR = Ver4.LegK[INR];
  momentum *OutR = Ver4.LegK[OUTR];

  //////// Initialize all sub-vertex ///////////

  array<momentum *, 4> LegK[10];
  vector<channel> ALL = {I, U, S, T};

  // LD Vertex
  LegK[0] = {InL, G[1].K, G[3].K, G[0].K};
  Env.Ver[0] = Vertex(LegK[0], LDInTL, 0, LoopIndex, ALL, Ver4.Type, LEFT);

  // LU Vertex
  LegK[1] = {G[1].K, OutL, G[2].K, G[4].K};
  LegK[2] = {G[1].K, OutR, G[2].K, G[6].K};
  Env.Ver[1] = Vertex(LegK[1], LUInTL, 0, LoopIndex, ALL, Ver4.Type, LEFT);
  Env.Ver[2] = Vertex(LegK[2], LUInTL, 0, LoopIndex, ALL, Ver4.Type, LEFT);

  // RD Vertex
  LegK[3] = {G[0].K, G[2].K, InR, G[5].K};
  LegK[4] = {G[0].K, G[2].K, G[7].K, OutR};
  LegK[5] = {G[0].K, G[2].K, G[8].K, OutL};
  for (int i = 3; i <= 5; i++)
    Env.Ver[i] = Vertex(LegK[i], RDInTL, 0, LoopIndex, ALL, Ver4.Type, RIGHT);

  // RU Vertex
  LegK[6] = {G[4].K, G[3].K, G[5].K, OutR};
  LegK[7] = {G[6].K, G[3].K, G[5].K, OutL};
  LegK[8] = {G[4].K, G[3].K, InR, G[7].K};
  LegK[9] = {G[6].K, G[3].K, InR, G[8].K};
  for (int i = 6; i <= 9; i++)
    Env.Ver[i] = Vertex(LegK[i], RUInTL, 0, LoopIndex, ALL, Ver4.Type, RIGHT);

  //////// T map (for all four envelope diagram) //////
  // four diagrams have the same sub-vertex Tau configuration
  // so here we just use the first diagram
  Env.Map = CreateMapT4(Ver4, Env.Ver[0], Env.Ver[1], Env.Ver[3], Env.Ver[6]);

  Env.SymFactor = {-1.0, 1.0, -1.0, 1.0};
  Ver4.Envelope.push_back(Env);
  return Ver4;
}

string verDiag::ToString(const ver4 &Ver4) {
  string Info = fmt::format("Root: \n  ID: {0}; T: ", Ver4.ID);
  for (auto &t : Ver4.T)
    Info +=
        fmt::format("({0}, {1}, {2}, {3}), ", t[INL], t[OUTL], t[INR], t[OUTR]);
  Info += "\n";
  for (auto &bubble : Ver4.Bubble) {
    Info += "SubVer: \n";
    for (int p = 0; p < bubble.Pair.size(); p++) {
      pair pp = bubble.Pair[p];
      Info += fmt::format("  LVer ID: {0}, T: ", pp.LVer[0].ID);
      for (auto &t : pp.LVer[0].T)
        Info += fmt::format("({0}, {1}, {2}, {3}), ", t[INL], t[OUTL], t[INR],
                            t[OUTR]);
      Info += "\n";

      Info += fmt::format("  RVer ID: {0}, T: ", pp.RVer[0].ID);
      for (auto &t : pp.RVer[0].T)
        Info += fmt::format("({0}, {1}, {2}, {3}), ", t[INL], t[OUTL], t[INR],
                            t[OUTR]);
      Info += "\n";

      // Info += fmt::format("  G1 Internal T Map: ");
      // for (auto &m : pp.Map)
      //   Info += fmt::format("({0}, {1}): {2}-{3}, ", m.LVerTidx, m.RVerTidx,
      //   m.GT[0],
      //                       m.G1T[1]);
      // Info += "\n";

      // Info += fmt::format("  G2 Internal T Map: ");
      // for (auto &m : pp.Map)
      //   Info += fmt::format("({0}, {1}): {2}-{3}, ", m.LVerT, m.RVerT,
      //   m.G2T[0],
      //                       m.G2T[1]);
      // Info += "\n";

      Info += fmt::format("         Map:        ");
      for (auto &m : pp.Map)
        Info +=
            fmt::format("({0}, {1}) => {2}, ", m.LVerTidx, m.RVerTidx, m.Tidx);
      Info += "\n";
    }
  }
  return Info;
}

bool verTest() {
  verDiag VerDiag;
  vector<channel> Chan = {T, U, S};
  //   ver4 *Root = VerDiag.Build(1, Chan, NORMAL);
  //   cout << VerDiag.ToString(*Root) << endl;
  return true;
}
