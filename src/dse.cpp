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
  Ver4.InT[INL] = InT[INL];                                                    \
  Ver4.InT[INR] = InT[INR];

#define SETOUTT(Ver4, DiagNum, OutL, OutR)                                     \
  Ver4.OutT[DiagNum][OUTL] = OutL;                                             \
  Ver4.OutT[DiagNum][OUTR] = OutR;

#define TIND(LTau, RTau) (LTau * MaxTauNum + RTau)

ver4 *verDiag::Ver0(array<int, 2> InT, int LoopNum, vertype Type) {
  VerPool.push_back(ver4());
  ver4 &Ver4 = VerPool.back();
  //   Ver4.Channel = 0;
  Ver4.Weight.fill(0.0);
  SETINT(Ver4, InT);
  Ver4.LoopNum = 0;
  ////////////// bare interaction ///////////

  Ver4.Pairs.push_back(pair{nullptr, nullptr, IDIR, BARE});

  Ver4.Pairs.push_back(pair{nullptr, nullptr, IEX, BARE});

  if (Type == RENORMALIZED) {
    //////////// dressed interaction ///////////
    Ver4.Pairs.push_back(pair{nullptr, nullptr, IDIR, RENORMALIZED});

    Ver4.Pairs.push_back(pair{nullptr, nullptr, IEX, RENORMALIZED});
  }

  // construct possible OutT pairs
  Ver4.OutT.push_back({InT[LEFT], InT[RIGHT]});
  Ver4.OutT.push_back({InT[RIGHT], InT[LEFT]});

  return &VerPool.back();
}

ver4 *verDiag::ChanI(array<int, 2> InT, int LoopNum, vertype Type) {
  if (LoopNum == 0)
    return Ver0(InT, LoopNum, Type);
}

ver4 *verDiag::Bubble(array<int, 2> InT, int LoopNum, vector<channel> Channel,
                      vertype Type) {
  VerPool.push_back(ver4());
  ver4 &Ver4 = VerPool.back();
  Ver4.LoopNum = LoopNum;
  Ver4.Weight.fill(0.0);

  ASSERT_ALLWAYS(InT[RIGHT] - InT[LEFT] == 2 * (LoopNum + 1),
                 "Tau Number must be 2*(loopNum+1)");

  for (auto &chan : Channel) {
    for (int ol = 0; ol < LoopNum - 1; ol++) {
      int LTauNum = 2 * (ol + 1); //-2 because left and right InT are known
      array<int, 2> LInT = {InT[LEFT], InT[LEFT] + LTauNum};

      ver4 *LVer[2] = {nullptr, nullptr};
      ver4 *RVer[2] = {nullptr, nullptr};

      if (ol == 0) {
        LVer[0] = ChanI(LInT, ol, BARE);
      } else {
        LVer[0] = ChanI(LInT, ol, RENORMALIZED);
        if (chan == U || chan == T)
          LVer[1] = Bubble(LInT, ol, {U, S}, RENORMALIZED);
        else
          LVer[1] = Bubble(LInT, ol, {U, T}, RENORMALIZED);
      }

      int oR = LoopNum - 1 - ol;
      array<int, 2> RInT = {InT[LEFT] + LTauNum + 1, InT[RIGHT]};

      if (oR == 0) {
        RVer[0] = ChanI(RInT, oR, BARE);
      } else {
        RVer[0] = ChanI(RInT, oR, RENORMALIZED);
        RVer[1] = Bubble(RInT, oR, {U, S, T}, RENORMALIZED);
      }

      for (int i = 0; i < 2; i++) {
        if (LVer[i] != nullptr)
          Ver4.SubVer.push_back(LVer[i]);
        if (RVer[i] != nullptr)
          Ver4.SubVer.push_back(RVer[i]);
      }

      for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
          if ((LVer[i] != nullptr) && (RVer[j] != nullptr))
            Ver4.Pairs.push_back(pair{LVer[i], RVer[j], chan, Type});
    }
  }

  // find all independent tau
  for (auto &pair : Ver4.Pairs) {
    auto LOutTVec = pair.LVer->OutT;
    auto ROutTVec = pair.RVer->OutT;
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
  return &VerPool.back();
}