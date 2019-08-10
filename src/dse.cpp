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
  return Bubble(0, LoopNum, Channel, Type);
}

ver4 verDiag::Ver0(int InTL, vertype Type) {
  ver4 Ver4;
  Ver4.ID = DiagNum;
  DiagNum++;
  //   Ver4.Channel = 0;
  Ver4.LoopNum = 0;
  Ver4.Type = Type;
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

ver4 verDiag::ChanI(int InTL, int LoopNum, vertype Type) {
  if (LoopNum == 0)
    return Ver0(InTL, Type);
}

void AddToTList(vector<array<int, 4>> &TList, const array<int, 4> T) {
  bool Flag = true;
  for (auto &t : TList) {
    ASSERT_ALLWAYS(t[INL] == T[INL],
                   "left Tin must be the same for all subvertex!");
    if (t[OUTL] == T[OUTL] && t[INR] == T[INR] && t[OUTR] == T[OUTR])
      Flag == false;
  }
  if (Flag)
    TList.push_back(T);
}

ver4 verDiag::Bubble(int InTL, int LoopNum, vector<channel> Channel,
                     vertype Type) {
  ver4 Ver4;
  Ver4.ID = DiagNum;
  DiagNum++;
  Ver4.LoopNum = LoopNum;
  Ver4.Type = Type;
  Ver4.Channel = Channel;

  for (int i = 0; i < pow(2 * (LoopNum + 1), 3); i++)
    Ver4.Weight.push_back(0.0);

  for (auto &chan : Channel) {
    // bubble subroutine does not allow I channel
    if (chan == I)
      continue;
    for (int ol = 0; ol < LoopNum; ol++) {
      int LTauNum = 2 * (ol + 1);

      vector<ver4> LVer, RVer;

      ////////////////////   Left SubVer  ///////////////////
      if (ol == 0) {
        LVer.push_back(ChanI(InTL, ol, BARE));
        //   cout << LVer[0]->OutT[0][LEFT] << LVer[0]->OutT[0][RIGHT] << endl;
        //  cout << LVer[0]->ID << endl;
        //  cout << LVer[0]->OutT.size() << endl;
      } else {
        LVer.push_back(ChanI(InTL, ol, RENORMALIZED));
        if (chan == U || chan == T)
          LVer.push_back(Bubble(InTL, ol, {U, S}, RENORMALIZED));
        else
          LVer.push_back(Bubble(InTL, ol, {U, T}, RENORMALIZED));
      }

      ////////////////////   Right SubVer  ///////////////////
      int oR = LoopNum - 1 - ol;
      int RInTL = InTL + LTauNum;
      if (oR == 0) {
        RVer.push_back(ChanI(RInTL, oR, BARE));
      } else {
        RVer.push_back(ChanI(RInTL, oR, RENORMALIZED));
        RVer.push_back(Bubble(RInTL, oR, {U, S, T}, RENORMALIZED));
      }

      ////////////////////   Merge SubVer  ///////////////////
      vector<int> LVerIndex, RVerIndex;
      for (auto &ver : LVer) {
        Ver4.SubVer.push_back(ver);
        LVerIndex.push_back(Ver4.SubVer.size() - 1);
      }
      for (auto &ver : RVer) {
        Ver4.SubVer.push_back(ver);
        RVerIndex.push_back(Ver4.SubVer.size() - 1);
      }

      ////////////////////   External Tau  ///////////////////
      for (auto &i : LVerIndex)
        for (auto &j : RVerIndex) {
          Ver4.Pairs.push_back(pair{i, j, chan, Type});
          auto LVerTList = Ver4.SubVer[i].T;
          auto RVerTList = Ver4.SubVer[j].T;
          for (auto &LVerT : LVerTList)
            for (auto &RVerT : RVerTList) {
              array<int, 4> LegT;
              if (chan == T) {
                LegT = {LVerT[INL], LVerT[OUTL], RVerT[INR], RVerT[OUTR]};
              } else if (chan == U) {
                LegT = {LVerT[INL], RVerT[OUTR], RVerT[INR], LVerT[OUTL]};
              } else if (chan == S) {
                LegT = {LVerT[INL], RVerT[OUTL], LVerT[INR], RVerT[OUTR]};
              }
              AddToTList(Ver4.T, LegT);
            }
        }
    }
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