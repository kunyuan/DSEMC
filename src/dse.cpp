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

int verDiag::NewMom() {
  MomNum += 1;
  ASSERT_ALLWAYS(MomNum < MaxMomNum, "Too many momentum variables! " << MomNum);
  return MomNum - 1;
}

ver4 verDiag::Build(int LoopNum, vector<channel> Channel, vertype Type) {
  ASSERT_ALLWAYS(LoopNum > 0, "LoopNum must be larger than zero!");
  DiagNum = 0;
  MomNum = MaxLoopNum;
  array<int, 4> LegK = {1, NewMom(), 2, NewMom()};
  return Vertex(LegK, 0, LoopNum, 3, Channel, Type, LEFT);
}

ver4 verDiag::Vertex(array<int, 4> LegK, int InTL, int LoopNum, int LoopIndex,
                     vector<channel> Channel, vertype Type, int Side) {
  ver4 Ver4;
  Ver4.ID = DiagNum;
  DiagNum++;
  Ver4.LoopNum = LoopNum;
  Ver4.TauNum = 2 * (LoopNum + 1);
  Ver4.Type = Type;
  Ver4.LegK = LegK;

  if (LoopNum == 0) {
    ASSERT_ALLWAYS(Channel[0] == I, "Only I channel has zero loop vertex!");
    Ver4 = Ver0(Ver4, InTL, Type, Side);
  } else {

    Ver4.Channel = Channel;
    Ver4.K[0] = LoopIndex;
    Ver4.G[0] = gMatrix(Ver4.TauNum, InTL);

    for (auto &chan : Channel) {
      if (chan == I)
        Ver4 = ChanI(Ver4, InTL, LoopNum, LoopIndex, Type, Side);
      else
        Ver4 = ChanUST(Ver4, InTL, LoopNum, LoopIndex, chan, Type, Side);
    }
  }

  Ver4.Weight.resize(Ver4.T.size());
  for (auto &d : Ver4.Weight)
    d = 0.0;

  return Ver4;
}

ver4 verDiag::Ver0(ver4 Ver4, int InTL, vertype Type, int Side) {
  ////////////// bare interaction ///////////
  if (Side == LEFT)
    // Side==left, then make sure INL Tau are the last TauIndex
    Ver4.T.push_back({InTL, InTL, InTL, InTL});
  else
    // Side==right, then make sure INR Tau are the last TauIndex
    Ver4.T.push_back({InTL + 1, InTL + 1, InTL + 1, InTL + 1});

  if (Type != BARE) {
    //////////// dressed interaction ///////////
    Ver4.T.push_back({InTL, InTL, InTL + 1, InTL + 1});
    Ver4.T.push_back({InTL, InTL + 1, InTL + 1, InTL});
  }
  return Ver4;
}

ver4 verDiag::ChanUST(ver4 Ver4, int InTL, int LoopNum, int LoopIndex,
                      channel chan, vertype Type, int Side) {

  ASSERT_ALLWAYS(chan != I, "ChanUST can not process I channel!");

  Ver4.K[chan] = NewMom();
  Ver4.G[chan] = gMatrix(Ver4.TauNum, InTL);

  array<int, 4> LLegK, RLegK; // left and right vertex external LegK
  if (chan == T) {

    LLegK = {Ver4.LegK[INL], Ver4.LegK[OUTL], Ver4.K[chan], Ver4.K[0]};
    RLegK = {Ver4.K[0], Ver4.K[chan], Ver4.LegK[INR], Ver4.LegK[OUTR]};

  } else if (chan == U) {

    LLegK = {Ver4.LegK[INL], Ver4.LegK[OUTR], Ver4.K[chan], Ver4.K[0]};
    RLegK = {Ver4.K[0], Ver4.K[chan], Ver4.LegK[INR], Ver4.LegK[OUTL]};

  } else if (chan == S) {

    LLegK = {Ver4.LegK[INL], Ver4.K[0], Ver4.LegK[INR], Ver4.K[chan]};
    RLegK = {Ver4.K[0], Ver4.LegK[OUTL], Ver4.K[chan], Ver4.LegK[OUTR]};
  }

  ver4 LVer, RVer;
  for (int ol = 0; ol < LoopNum; ol++) {

    ////////////////////   Left SubVer  ///////////////////
    if (chan == U || chan == T)
      LVer = Vertex(LLegK, InTL, ol, LoopIndex + 1, {I, U, S}, Type, LEFT);
    else
      LVer = Vertex(LLegK, InTL, ol, LoopIndex + 1, {I, U, T}, Type, LEFT);

    ////////////////////   Right SubVer  ///////////////////
    int oR = LoopNum - 1 - ol;
    int RInTL = InTL + 2 * (ol + 1);
    int RLoopNum = LoopIndex + 1 + ol;
    RVer = Vertex(RLegK, RInTL, oR, RLoopNum, {I, U, S, T}, Type, RIGHT);

    ///////////   External and Internal Tau  ////////////////
    array<int, 2> G1T, G2T;
    vector<mapT> Map;

    for (int lt = 0; lt < LVer.T.size(); ++lt)
      for (int rt = 0; rt < RVer.T.size(); ++rt) {

        auto &LvT = LVer.T[lt];
        auto &RvT = RVer.T[rt];
        array<int, 4> LegT;

        if (chan == T) {

          LegT = {LvT[INL], LvT[OUTL], RvT[INR], RvT[OUTR]};
          G1T = {LvT[OUTR], RvT[INL]};
          G2T = {RvT[OUTL], LvT[INR]};

        } else if (chan == U) {

          LegT = {LvT[INL], RvT[OUTR], RvT[INR], LvT[OUTL]};
          G1T = {LvT[OUTR], RvT[INL]};
          G2T = {RvT[OUTL], LvT[INR]};

        } else if (chan == S) {

          LegT = {LvT[INL], RvT[OUTL], LvT[INR], RvT[OUTR]};
          G1T = {LvT[OUTL], RvT[INL]};
          G2T = {LvT[OUTR], RvT[INR]};
        }

        // add T array into the T pool of the vertex
        int Index = AddToTList(Ver4.T, LegT);
        Map.push_back(mapT{lt, rt, G1T, G2T, Index});
      }

    Ver4.Pairs.push_back(pair{LVer, RVer, chan, Sym(chan), Map});
  }
  return Ver4;
}

void CreateMapT(ver4 &Ver4, ver4 LDVer, ver4 LUVer, ver4 RDVer, ver4 RUVer,
                int Index) {
  vector<mapT4> Map;

  for (int ldt = 0; ldt < LDVer.T.size(); ldt++)
    for (int lut = 0; lut < LUVer.T.size(); lut++)
      for (int rdt = 0; rdt < RDVer.T.size(); rdt++)
        for (int rut = 0; rut < RUVer.T.size(); rut++) {
          auto &LDvT = LDVer.T[ldt];
          auto &LUvT = LUVer.T[lut];
          auto &RDvT = RDVer.T[rdt];
          auto &RUvT = RUVer.T[rut];
          array<int, 4> LegT = {LDvT[INL], LUvT[OUTL], RDvT[INR], RUvT[OUTR]};
          array<array<int, 2>, 6> GT;
          if (Index == 0) {
            // LD: 1, LU: 0, RD: 0, RU: 1
            GT[0] = {LDvT[OUTR], RDvT[INL]};
            GT[1] = {LDvT[OUTL], LUvT[INL]};
            GT[2] = {RUvT[OUTL], LDvT[INR]};
            GT[3] = {RDvT[OUTL], LUvT[INR]};
            GT[4] = {LUvT[OUTR], RUvT[INL]};
            GT[5] = {RDvT[OUTR], RUvT[INR]};
          } else if (Index == 1) {
            // LD: 0, LU: 0, RD: 1, RU: 2
            GT[0] = {RDvT[OUTR], LDvT[INR]};
            GT[1] = {LDvT[OUTL], LUvT[INL]};
            GT[2] = {LDvT[OUTR], RUvT[INR]};
            GT[3] = {RDvT[OUTL], LUvT[INR]};
            GT[4] = {LUvT[OUTR], RUvT[INL]};
            GT[5] = {RUvT[OUTL], RDvT[INL]};
          } else if (Index == 2) {
            // LD: 0, LU: 1, RD: 2, RU: 0
            GT[0] = {RDvT[OUTR], LDvT[INR]};
            GT[1] = {LDvT[OUTL], LUvT[INL]};
            GT[2] = {LDvT[OUTR], RUvT[INL]};
            GT[3] = {LUvT[OUTR], RDvT[INL]};
            GT[4] = {LUvT[OUTR], RUvT[INL]};
            GT[5] = {RDvT[OUTL], RUvT[INR]};
          } else if (Index == 3) {
            // LD: 2, LU: 2, RD: 0, RU: 0
            GT[0] = {LDvT[OUTR], RDvT[INL]};
            GT[1] = {LUvT[OUTR], LDvT[INR]};
            GT[2] = {LDvT[OUTL], RUvT[INL]};
            GT[3] = {RDvT[OUTL], LUvT[INL]};
            GT[4] = {RUvT[OUTL], LUvT[INR]};
            GT[5] = {RDvT[OUTR], RUvT[INR]};
          }

          int TIndex = AddToTList(Ver4.T, LegT);
          Map.push_back(mapT4{ldt, lut, rdt, rut, GT, TIndex});
        }
  Ver4.Pair4s.push_back(pair4{LDVer, LUVer, RDVer, RUVer, Map});
}

#define GETREF(ld, lu, rd, ru)

ver4 verDiag::ChanI(ver4 Ver4, int InTL, int LoopNum, int LoopIndex,
                    vertype Type, int Side) {

  if (LoopNum != 3)
    return;
  for (int i = 0; i < 6; i++) {
    // there are 14 independent G
    Ver4.Kip.push_back(NewMom());
    Ver4.Kin.push_back(NewMom());
    Ver4.Gi.push_back(gMatrix(Ver4.TauNum, InTL));
  }

  int InL = Ver4.LegK[INL];
  int OutL = Ver4.LegK[OUTL];
  int InR = Ver4.LegK[INR];
  int OutR = Ver4.LegK[OUTR];

  array<int, 4> LDLegK[4], LULegK[4], RDLegK[4], RULegK[4];

  LDLegK[0] = {InL, Ver4.Kip[1], Ver4.Kip[2], Ver4.Kip[0]};
  LDLegK[1] = {InL, Ver4.Kip[1], Ver4.Kin[0], Ver4.Kin[2]};
  LDLegK[2] = {InL, Ver4.Kip[1], Ver4.Kin[0], Ver4.Kin[2]};
  LDLegK[3] = {InL, Ver4.Kin[2], Ver4.Kin[1], Ver4.Kip[0]};

  LULegK[0] = {Ver4.Kip[1], OutL, Ver4.Kin[3], Ver4.Kip[4]};
  LULegK[1] = {Ver4.Kip[1], OutL, Ver4.Kin[3], Ver4.Kip[4]};
  LULegK[2] = {Ver4.Kip[1], OutL, Ver4.Kin[4], Ver4.Kip[3]};
  LULegK[3] = {Ver4.Kin[3], OutL, Ver4.Kin[4], Ver4.Kin[1]};

  RDLegK[0] = {Ver4.Kip[0], Ver4.Kin[3], InR, Ver4.Kip[5]};
  RDLegK[1] = {Ver4.Kin[5], Ver4.Kin[3], InR, Ver4.Kin[0]};
  RDLegK[2] = {Ver4.Kip[3], Ver4.Kip[5], InR, Ver4.Kin[0]};
  RDLegK[3] = {Ver4.Kip[0], Ver4.Kin[3], InR, Ver4.Kip[5]};

  RULegK[0] = {Ver4.Kip[4], Ver4.Kip[2], Ver4.Kip[5], OutL};
  RULegK[1] = {Ver4.Kip[4], Ver4.Kin[5], Ver4.Kin[2], OutL};
  RULegK[2] = {Ver4.Kin[2], Ver4.Kin[4], Ver4.Kip[5], OutL};
  RULegK[3] = {Ver4.Kin[2], Ver4.Kin[4], Ver4.Kip[5], OutL};

  vector<channel> ALL = {I, U, S, T};

  ver4 LDVer, LUVer, RDVer, RUVer;
  for (int ld = 0; ld <= LoopNum - 3; ld++)
    for (int lu = 0; lu <= LoopNum - 3; lu++)
      for (int rd = 0; rd <= LoopNum - 3; rd++)
        for (int ru = 0; ru <= LoopNum - 3; ru++) {
          int ldloopidx = LoopIndex + 3;
          int luloopidx = LoopIndex + 3 + ld;
          int rdloopidx = LoopIndex + 3 + ld + lu;
          int ruloopidx = LoopIndex + 3 + ld + lu + rd;
          int ldTL = InTL;
          int luTL = InTL + 2 * (ld + 1);
          int rdTL = InTL + 2 * (ld + lu + 1);
          int ruTL = InTL + 2 * (ld + lu + ru + 1);
          for (int t = 0; t < 4; t++) {
            LDVer = Vertex(LDLegK[t], ldTL, ld, ldloopidx, ALL, Type, LEFT);
            LUVer = Vertex(LULegK[t], luTL, lu, luloopidx, ALL, Type, LEFT);
            RDVer = Vertex(RDLegK[t], rdTL, rd, rdloopidx, ALL, Type, RIGHT);
            RUVer = Vertex(RULegK[t], ruTL, ru, ruloopidx, ALL, Type, RIGHT);
            CreateMapT(Ver4, LDVer, LUVer, RDVer, RUVer, t);
          }
        }
  return Ver4;
}

string verDiag::ToString(const ver4 &Ver4) {
  string Info = fmt::format("Root: \n  ID: {0}; T: ", Ver4.ID);
  for (auto &t : Ver4.T)
    Info +=
        fmt::format("({0}, {1}, {2}, {3}), ", t[INL], t[OUTL], t[INR], t[OUTR]);
  Info += "\n";
  Info += "SubVer: \n";
  for (int p = 0; p < Ver4.Pairs.size(); p++) {
    pair pp = Ver4.Pairs[p];
    Info += fmt::format("  LVer ID: {0}, T: ", pp.LVer.ID);
    for (auto &t : pp.LVer.T)
      Info += fmt::format("({0}, {1}, {2}, {3}), ", t[INL], t[OUTL], t[INR],
                          t[OUTR]);
    Info += "\n";

    Info += fmt::format("  RVer ID: {0}, T: ", pp.RVer.ID);
    for (auto &t : pp.RVer.T)
      Info += fmt::format("({0}, {1}, {2}, {3}), ", t[INL], t[OUTL], t[INR],
                          t[OUTR]);
    Info += "\n";

    Info += fmt::format("  G1 Internal T Map: ");
    for (auto &m : pp.Map)
      Info += fmt::format("({0}, {1}): {2}-{3}, ", m.LVerT, m.RVerT, m.G1T[0],
                          m.G1T[1]);
    Info += "\n";

    Info += fmt::format("  G2 Internal T Map: ");
    for (auto &m : pp.Map)
      Info += fmt::format("({0}, {1}): {2}-{3}, ", m.LVerT, m.RVerT, m.G2T[0],
                          m.G2T[1]);
    Info += "\n";

    Info += fmt::format("         Map:        ");
    for (auto &m : pp.Map)
      Info += fmt::format("({0}, {1}) => {2}, ", m.LVerT, m.RVerT, m.T);
    Info += "\n";
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
