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

// #define TARRAY(LoopNum, InL, OutL, InR, OutR)                                  \
//   {                                                                            \
//     InL, OutL, InR, OutR,                                                      \
//         (OutL - InL) * pow(2 * (LoopNum + 1), 2) +                             \
//             (InR - InL) * 2 * (LoopNum + 1) + (OutR - InL)                     \
//   }

ver4 verDiag::Build(int LoopNum, vector<channel> Channel, vertype Type) {
  ASSERT_ALLWAYS(LoopNum > 0, "LoopNum must be larger than zero!");
  return Vertex(0, LoopNum, Channel, Type);
}

ver4 verDiag::Vertex(int InTL, int LoopNum, vector<channel> Channel,
                     vertype Type) {
  ver4 Ver4;
  Ver4.ID = DiagNum;
  DiagNum++;
  Ver4.LoopNum = LoopNum;
  Ver4.Type = Type;
  if (LoopNum == 0)
    return Ver0(Ver4, InTL, Type);

  Ver4.Channel = Channel;
  for (auto &chan : Channel) {
    if (chan == I)
      Ver4 = ChanI(Ver4, InTL, LoopNum, Type);
    else
      Ver4 = Bubble(Ver4, InTL, LoopNum, chan, Type);
  }
}

ver4 verDiag::Ver0(ver4 Ver4, int InTL, vertype Type) {
  // ver4 Ver4;
  // Ver4.ID = DiagNum;
  // DiagNum++;
  //   Ver4.Channel = 0;
  // Ver4.LoopNum = 0;
  // Ver4.Type = Type;
  ////////////// bare interaction ///////////
  Ver4.T.push_back({InTL, InTL, InTL, InTL});
  Ver4.Weight.push_back(0.0);

  if (Type == RENORMALIZED) {
    //////////// dressed interaction ///////////
    // construct possible OutT pairs
    Ver4.T.push_back({InTL, InTL, InTL + 1, InTL + 1});
    Ver4.T.push_back({InTL, InTL + 1, InTL + 1, InTL});
    Ver4.Weight.push_back(0.0);
    Ver4.Weight.push_back(0.0);
  }

  return Ver4;
}

ver4 verDiag::ChanI(ver4 Ver4, int InTL, int LoopNum, vertype Type) {
  return Ver4;
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

ver4 verDiag::Bubble(ver4 Ver4, int InTL, int LoopNum, channel chan,
                     vertype Type) {
  // ver4 Ver4;
  // Ver4.ID = DiagNum;
  // DiagNum++;
  // Ver4.LoopNum = LoopNum;
  // Ver4.Type = Type;
  // Ver4.Channel = Channel;
  if (chan == I)
    return Ver4;

  for (int i = 0; i < pow(2 * (LoopNum + 1), 3); i++)
    Ver4.Weight.push_back(0.0);

  for (int ol = 0; ol < LoopNum; ol++) {
    int LTauNum = 2 * (ol + 1);

    ver4 LVer, RVer;

    ////////////////////   Left SubVer  ///////////////////
    if (chan == U || chan == T)
      LVer = Vertex(InTL, ol, {I, U, S}, RENORMALIZED);
    else
      LVer = Vertex(InTL, ol, {I, U, T}, RENORMALIZED);

    ////////////////////   Right SubVer  ///////////////////
    int oR = LoopNum - 1 - ol;
    int RInTL = InTL + LTauNum;
    RVer = Vertex(RInTL, oR, {I, U, S, T}, RENORMALIZED);

    ////////////////////   Merge SubVer  ///////////////////
    Ver4.SubVer.push_back(LVer);
    Ver4.SubVer.push_back(RVer);

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

    Ver4.Pairs.push_back(
        pair{Ver4.SubVer.size() - 2, Ver4.SubVer.size() - 1, Map});
  }
  return Ver4;
}

string verDiag::ToString(const ver4 &Ver4) {
  string Info = fmt::format("Root: \n  ID: {0}; OutT: ", Ver4.ID);
  for (auto &t : Ver4.T)
    Info +=
        fmt::format("({0}, {1}, {2}, {3}), ", t[INL], t[OUTL], t[INR], t[OUTR]);
  Info += "\n";
  Info += "SubVer: \n";
  for (auto &sub : Ver4.SubVer) {
    Info += fmt::format("  Sub ID: {0}, OutT: ", sub.ID);
    for (auto &t : sub.T)
      Info += fmt::format("({0}, {1}, {2}, {3}), ", t[INL], t[OUTL], t[INR],
                          t[OUTR]);
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