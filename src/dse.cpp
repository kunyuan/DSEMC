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

int Sym(channel chan) {
  if (chan == T)
    return -1.0;
  else if (chan == U)
    return 1.0;
  else if (chan == S)
    return 0.5;
  else
    return 1.0;
}

int verDiag::NextMom() {
  MomNum += 1;
  return MomNum - 1;
}

ver4 verDiag::Build(int LoopNum, vector<channel> Channel, vertype Type) {
  ASSERT_ALLWAYS(LoopNum > 0, "LoopNum must be larger than zero!");
  DiagNum = 0;
  MomNum = MaxLoopNum;
  array<int, 4> LegK = {1, NextMom(), 2, NextMom()};
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
    Ver4.K1 = LoopIndex;
    //     Ver4.G1.resize(pow(Ver4.TauNum - 2, 2));

    for (auto &chan : Channel) {
      if (chan == I)
        Ver4 = ChanI(Ver4, InTL, LoopNum, LoopIndex, Type, Side);
      else {
        Ver4 = ChanUST(Ver4, InTL, LoopNum, LoopIndex, chan, Type, Side);
        Ver4.K2[chan] = NextMom();
        //  Ver4.G2[chan].resize(pow(Ver4.TauNum - 2, 2));
      }
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

ver4 verDiag::ChanI(ver4 Ver4, int InTL, int LoopNum, int LoopIndex,
                    vertype Type, int Side) {
  return Ver4;
}

ver4 verDiag::ChanUST(ver4 Ver4, int InTL, int LoopNum, int LoopIndex,
                      channel chan, vertype Type, int Side) {
  ASSERT_ALLWAYS(chan != I, "ChanUST can not process I channel!");
  ver4 LVer, RVer;
  array<int, 4> LLegK, RLegK;
  if (chan == T) {
    LLegK = {Ver4.LegK[INL], Ver4.LegK[OUTL], Ver4.K2[chan], Ver4.K1};
    RLegK = {Ver4.K1, Ver4.K2[chan], Ver4.LegK[INR], Ver4.LegK[OUTR]};
  } else if (chan == U) {
    LLegK = {Ver4.LegK[INL], Ver4.LegK[OUTR], Ver4.K2[chan], Ver4.K1};
    RLegK = {Ver4.K1, Ver4.K2[chan], Ver4.LegK[INR], Ver4.LegK[OUTL]};
  } else if (chan == S) {
    LLegK = {Ver4.LegK[INL], Ver4.K1, Ver4.LegK[INR], Ver4.K2[chan]};
    RLegK = {Ver4.K1, Ver4.LegK[OUTL], Ver4.K2[chan], Ver4.LegK[OUTR]};
  }

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

    ////////////////////   External Tau  ///////////////////
    map Map(LVer.T.size(), RVer.T.size());
    vector<array<int, 4>> InterTList;
    InterTList.resize(LVer.T.size() * RVer.T.size());
    for (int lt = 0; lt < LVer.T.size(); ++lt)
      for (int rt = 0; rt < RVer.T.size(); ++rt) {
        auto &LVerT = LVer.T[lt];
        auto &RVerT = RVer.T[rt];
        array<int, 4> LegT, InterT;
        if (chan == T) {
          LegT = {LVerT[INL], LVerT[OUTL], RVerT[INR], RVerT[OUTR]};
          InterT = {LVerT[OUTR], RVerT[INL], RVerT[OUTL], LVerT[INR]};
        } else if (chan == U) {
          LegT = {LVerT[INL], RVerT[OUTR], RVerT[INR], LVerT[OUTL]};
          InterT = {LVerT[OUTR], RVerT[INL], RVerT[OUTL], LVerT[INR]};
        } else if (chan == S) {
          LegT = {LVerT[INL], RVerT[OUTL], LVerT[INR], RVerT[OUTR]};
          InterT = {LVerT[OUTL], RVerT[INL], LVerT[OUTR], RVerT[INR]};
        }
        InterTList[lt * RVer.T.size() + rt] = InterT;
        int Index = AddToTList(Ver4.T, LegT);
        Map.Set(lt, rt, Index);
      }

    Ver4.Pairs.push_back(pair{LVer, RVer, InterTList, Map, chan, Sym(chan)});
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

    Info += fmt::format("  Map: ");
    for (int i = 0; i < pp.LVer.T.size(); i++)
      for (int j = 0; j < pp.RVer.T.size(); j++)
        Info += fmt::format("({0}, {1})>{2}, ", i, j, pp.Map.Get(i, j));
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