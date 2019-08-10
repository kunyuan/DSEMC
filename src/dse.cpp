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

#define SETINT(Ver4, InT)                                                      \
  Ver4.InT[LEFT] = InT[LEFT];                                                  \
  Ver4.InT[RIGHT] = InT[RIGHT];

#define SETOUTT(Ver4, DiagNum, OutL, OutR)                                     \
  Ver4.OutT[DiagNum][LEFT] = OutL;                                             \
  Ver4.OutT[DiagNum][RIGHT] = OutR;

#define TIND(LTau, RTau) (LTau * MaxTauNum + RTau)

ver4 verDiag::Build(int LoopNum, vector<channel> Channel, vertype Type) {
  ASSERT_ALLWAYS(LoopNum > 0, "LoopNum must be larger than zero!");
  array<int, 2> InT = {0, 2 * (LoopNum + 1) - 1};
  return Bubble(InT, LoopNum, Channel, Type);
}

ver4 verDiag::Ver0(array<int, 2> InT, int LoopNum, vertype Type) {
  ver4 Ver4;
  Ver4.ID = DiagNum;
  DiagNum++;
  //   Ver4.Channel = 0;
  Ver4.LoopNum = 0;
  Ver4.Weight.fill(0.0);
  SETINT(Ver4, InT);
  ////////////// bare interaction ///////////

  Ver4.Pairs.push_back(pair{-1, -1, IDIR, BARE});

  Ver4.Pairs.push_back(pair{-1, -1, IEX, BARE});

  if (Type == RENORMALIZED) {
    //////////// dressed interaction ///////////
    Ver4.Pairs.push_back(pair{-1, -1, IDIR, RENORMALIZED});

    Ver4.Pairs.push_back(pair{-1, -1, IEX, RENORMALIZED});
  }

  //   cout << Ver4.ID << ", " << InT[LEFT] << ", " << InT[RIGHT] << endl;
  // construct possible OutT pairs
  Ver4.OutT.push_back({InT[LEFT], InT[RIGHT]});
  Ver4.OutT.push_back({InT[RIGHT], InT[LEFT]});

  return Ver4;
}

ver4 verDiag::ChanI(array<int, 2> InT, int LoopNum, vertype Type) {
  if (LoopNum == 0)
    return Ver0(InT, LoopNum, Type);
}

ver4 verDiag::Bubble(array<int, 2> InT, int LoopNum, vector<channel> Channel,
                     vertype Type) {
  ver4 Ver4;
  Ver4.ID = DiagNum;
  DiagNum++;
  Ver4.LoopNum = LoopNum;
  Ver4.Weight.fill(0.0);
  SETINT(Ver4, InT);

  ASSERT_ALLWAYS(InT[RIGHT] - InT[LEFT] == 2 * (LoopNum + 1) - 1,
                 fmt::format("Tau Number must be 2*(loopNum+1), InT: {0}, {1}",
                             InT[LEFT], InT[RIGHT]));

  for (auto &chan : Channel) {
    for (int ol = 0; ol < LoopNum; ol++) {
      int LTauNum = 2 * (ol + 1); //-2 because left and right InT are known
      array<int, 2> LInT = {InT[LEFT], InT[LEFT] + LTauNum - 1};

      vector<ver4> LVer, RVer;

      if (ol == 0) {
        LVer.push_back(ChanI(LInT, ol, BARE));
        //   cout << LVer[0]->OutT[0][LEFT] << LVer[0]->OutT[0][RIGHT] << endl;
        //  cout << LVer[0]->ID << endl;
        //  cout << LVer[0]->OutT.size() << endl;
      } else {
        LVer.push_back(ChanI(LInT, ol, RENORMALIZED));
        if (chan == U || chan == T)
          LVer.push_back(Bubble(LInT, ol, {U, S}, RENORMALIZED));
        else
          LVer.push_back(Bubble(LInT, ol, {U, T}, RENORMALIZED));
      }

      int oR = LoopNum - 1 - ol;
      array<int, 2> RInT = {InT[LEFT] + LTauNum, InT[RIGHT]};

      if (oR == 0) {
        RVer.push_back(ChanI(RInT, oR, BARE));
      } else {
        RVer.push_back(ChanI(RInT, oR, RENORMALIZED));
        RVer.push_back(Bubble(RInT, oR, {U, S, T}, RENORMALIZED));
      }

      //       cout << LVer[0] << ", " << RVer[0] << endl;

      vector<int> LVerIndex, RVerIndex;
      for (auto &ver : LVer) {
        Ver4.SubVer.push_back(ver);
        LVerIndex.push_back(Ver4.SubVer.size() - 1);
      }
      for (auto &ver : RVer) {
        Ver4.SubVer.push_back(ver);
        RVerIndex.push_back(Ver4.SubVer.size() - 1);
      }

      cout << LVer[0].ID << ", " << LVer[0].OutT.size() << endl;

      for (auto &i : LVerIndex)
        for (auto &j : RVerIndex)
          //      cout << LVer[i]->ID << ", " << RVer[j]->ID << endl;
          //      cout << LVer[i]->OutT.size() << ", " << RVer[j]->OutT.size()
          //           << endl;
          Ver4.Pairs.push_back(pair{i, j, chan, Type});
    }
  }

  //   cout << Ver4.Pairs.size() << endl;

  // find all independent tau
  for (auto &pair : Ver4.Pairs) {
    auto LOutTVec = VerPool[pair.LVer].OutT;
    auto ROutTVec = VerPool[pair.RVer].OutT;
    for (auto &LOutT : LOutTVec)
      for (auto &ROutT : ROutTVec) {
        int OutTL = LOutT[LEFT];
        int OutTR = ROutT[RIGHT];
        bool Flag = false;
        for (auto outt : Ver4.OutT)
          if (outt[LEFT] == OutTL && outt[RIGHT] == OutTR)
            Flag = true;
        if (Flag == false)
          Ver4.OutT.push_back({OutTL, OutTR});
      }
  }
  return Ver4;
}

string verDiag::ToString(const ver4 &Ver4) {
  string Info =
      fmt::format("Root: \n  ID: {0}; InT: ({1}, {2}); OutT: ", Ver4.ID,
                  Ver4.InT[LEFT], Ver4.InT[RIGHT]);
  for (auto &outt : Ver4.OutT)
    Info += fmt::format("({0}, {1}), ", outt[LEFT], outt[RIGHT]);
  Info += "\n";
  Info += "SubVer: \n";
  for (auto &sub : Ver4.SubVer) {
    Info += fmt::format("  Sub ID: {0}, OutT: ", sub.ID);
    for (auto &outt : sub.OutT)
      Info += fmt::format("({0}, {1}), ", outt[LEFT], outt[RIGHT]);
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