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

vector<mapT4> CreateMapT(ver4 &Ver4, ver4 LDVer, ver4 LUVer, ver4 RDVer,
                         ver4 RUVer) {
  vector<mapT4> Map;
  array<array<int, 2>, 9> GT;
  for (int ldt = 0; ldt < 2; ldt++)
    for (int lut = 0; lut < 2; lut++)
      for (int rdt = 0; rdt < 2; rdt++)
        for (int rut = 0; rut < 2; rut++) {
          array<int, 4> LegT[4], Tindex;
          auto &ldT = LDVer.T[ldt];
          auto &luT = LUVer.T[lut];
          auto &rdT = RDVer.T[rdt];
          auto &ruT = RUVer.T[rut];

          GT[0] = {ldT[OUTR], rdT[INL]};
          GT[1] = {ldT[OUTL], luT[INL]};
          GT[2] = {ruT[OUTL], ldT[INR]};
          GT[3] = {rdT[OUTL], luT[INR]};
          GT[4] = {luT[OUTR], ruT[INL]};
          GT[5] = {rdT[OUTR], ruT[INR]};
          GT[6] = {luT[OUTR], ruT[INL]};
          GT[7] = {ruT[OUTR], rdT[INR]};
          GT[8] = {ruT[OUTR], rdT[INR]};

          LegT[0] = {ldT[INL], luT[OUTL], rdT[INR], ruT[OUTR]};
          LegT[1] = {ldT[INL], ruT[OUTR], rdT[INR], luT[OUTL]};
          LegT[2] = {ldT[INL], luT[OUTL], ruT[INR], rdT[OUTR]};
          LegT[3] = {ldT[INL], rdT[OUTR], ruT[INR], luT[OUTL]};

          for (int i = 0; i < 4; i++)
            Tindex[i] = AddToTList(Ver4.T, LegT[i]);

          Map.push_back(mapT4{ldt, lut, rdt, rut, GT, Tindex});
        }
  return Map;
}

ver4 verDiag::ChanI(ver4 Ver4, int InTL, int LoopNum, int LoopIndex,
                    vertype Type, int Side) {

  if (LoopNum != 3)
    return;
  envelope Env;
  Env.K[0] = LoopIndex;
  Env.K[1] = LoopIndex + 1;
  Env.K[2] = LoopIndex + 2;
  for (int i = 3; i < 9; i++)
    Env.K[i] = NewMom();

  int InL = Ver4.LegK[INL];
  int OutL = Ver4.LegK[OUTL];
  int InR = Ver4.LegK[INR];
  int OutR = Ver4.LegK[OUTR];

  array<int, 4> LegK[10];
  vector<channel> ALL = {I, U, S, T};

  // LD Vertex
  LegK[0] = {InL, Env.K[1], Env.K[2], Env.K[0]};
  Env.Ver[0] = Vertex(LegK[0], InTL, 0, LoopIndex, ALL, Type, LEFT);

  // LU Vertex
  LegK[1] = {Env.K[1], OutL, Env.K[3], Env.K[4]};
  LegK[2] = {Env.K[1], OutR, Env.K[3], Env.K[6]};
  Env.Ver[1] = Vertex(LegK[1], InTL + 2, 0, LoopIndex, ALL, Type, LEFT);
  Env.Ver[2] = Vertex(LegK[2], InTL + 2, 0, LoopIndex, ALL, Type, LEFT);

  // RD Vertex
  LegK[3] = {Env.K[0], Env.K[3], InR, Env.K[5]};
  LegK[4] = {Env.K[0], Env.K[3], Env.K[7], OutR};
  LegK[5] = {Env.K[0], Env.K[3], Env.K[8], OutL};
  for (int i = 3; i <= 5; i++)
    Env.Ver[i] = Vertex(LegK[i], InTL + 4, 0, LoopIndex, ALL, Type, RIGHT);

  // RU Vertex
  LegK[6] = {Env.K[4], Env.K[2], Env.K[5], OutR};
  LegK[7] = {Env.K[6], Env.K[2], Env.K[5], OutL};
  LegK[8] = {Env.K[4], Env.K[2], InR, Env.K[7]};
  LegK[9] = {Env.K[6], Env.K[2], InR, Env.K[8]};
  for (int i = 6; i <= 7; i++)
    Env.Ver[i] = Vertex(LegK[i], InTL + 6, 0, LoopIndex, ALL, Type, RIGHT);

  //T map
  Env.Map = CreateMapT(Ver4, Env.Ver[0], Env.Ver[1], Env.Ver[3], Env.Ver[6]);

  Ver4.Envelope.push_back(Env);
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
