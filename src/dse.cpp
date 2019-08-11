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

ver4 verDiag::Build(int LoopNum, vector<channel> Channel, vertype Type) {
  ASSERT_ALLWAYS(LoopNum > 0, "LoopNum must be larger than zero!");
  return Vertex(0, LoopNum, 3, Channel, Type, LEFT);
}

ver4 verDiag::Vertex(int InTL, int LoopNum, int LoopIndex,
                     vector<channel> Channel, vertype Type, int Side) {
  ver4 Ver4;
  Ver4.ID = DiagNum;
  DiagNum++;
  Ver4.LoopNum = LoopNum;
  Ver4.LoopIndex = LoopIndex;
  Ver4.TauNum = 2 * (LoopNum + 1);
  Ver4.Type = Type;

  if (LoopNum == 0)
    Ver4 = Ver0(Ver4, InTL, Type, Side);
  else {

    Ver4.Channel = Channel;
    for (auto &chan : Channel) {
      if (chan == I)
        Ver4 = ChanI(Ver4, InTL, LoopNum, LoopIndex, Type, Side);
      else
        Ver4 = Bubble(Ver4, InTL, LoopNum, LoopIndex, chan, Type, Side);
    }
  }
  Ver4.G1.resize(pow(Ver4.TauNum - 2, 2));
  Ver4.G2.resize(pow(Ver4.TauNum - 2, 2));
  Ver4.Weight.resize(Ver4.T.size());
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

  if (Type == RENORMALIZED) {
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

ver4 verDiag::Bubble(ver4 Ver4, int InTL, int LoopNum, int LoopIndex,
                     channel chan, vertype Type, int Side) {
  if (chan == I)
    return Ver4;

  for (int i = 0; i < pow(2 * (LoopNum + 1), 3); i++)
    Ver4.Weight.push_back(0.0);

  for (int ol = 0; ol < LoopNum; ol++) {
    int LTauNum = 2 * (ol + 1);

    ver4 LVer, RVer;

    ////////////////////   Left SubVer  ///////////////////
    if (chan == U || chan == T)
      LVer = Vertex(InTL, ol, LoopIndex + 1, {I, U, S}, RENORMALIZED, LEFT);
    else
      LVer = Vertex(InTL, ol, LoopIndex + 1, {I, U, T}, RENORMALIZED, LEFT);

    ////////////////////   Right SubVer  ///////////////////
    int oR = LoopNum - 1 - ol;
    int RInTL = InTL + LTauNum;
    RVer = Vertex(RInTL, oR, LoopIndex + 1 + ol, {I, U, S, T}, RENORMALIZED,
                  RIGHT);

    ////////////////////   External Tau  ///////////////////
    map Map(LVer.T.size(), RVer.T.size());
    for (int lt = 0; lt < LVer.T.size(); ++lt)
      for (int rt = 0; rt < RVer.T.size(); ++rt) {
        auto &LVerT = LVer.T[lt];
        auto &RVerT = RVer.T[rt];
        array<int, 4> LegT;
        if (chan == T) {
          LegT = {LVerT[INL], LVerT[OUTL], RVerT[INR], RVerT[OUTR]};
        } else if (chan == U) {
          LegT = {LVerT[INL], RVerT[OUTR], RVerT[INR], LVerT[OUTL]};
        } else if (chan == S) {
          LegT = {LVerT[INL], RVerT[OUTL], LVerT[INR], RVerT[OUTR]};
        }
        int Index = AddToTList(Ver4.T, LegT);
        Map.Set(lt, rt, Index);
      }

    Ver4.Pairs.push_back(pair{LVer, RVer, Map, chan, Sym(chan)});
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