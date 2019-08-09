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

ver4 *verDiag::Ver0(array<int, 2> InT, int LoopNum, vertype Type) {
  VerPool.push_back(ver4());
  ver4 &Ver4 = VerPool.back();
  //   Ver4.Channel = 0;
  Ver4.Weight.fill(0.0);
  Ver4.DiagNum = 0;
  SETINT(Ver4, InT);
  Ver4.LoopNum = 0;
  ////////////// bare interaction ///////////

  Ver4.Channel[Ver4.DiagNum] = IDIR;
  Ver4.Type[Ver4.DiagNum] = BARE;
  SETOUTT(Ver4, Ver4.DiagNum, InT[INL], InT[INR]);
  Ver4.DiagNum++;

  Ver4.Channel[Ver4.DiagNum] = IEX;
  Ver4.Type[Ver4.DiagNum] = BARE;
  SETOUTT(Ver4, Ver4.DiagNum, InT[INR], InT[INL]);
  Ver4.DiagNum++;

  if (Type == RENORMALIZED) {
    //////////// dressed interaction ///////////
    Ver4.Channel[Ver4.DiagNum] = IDIR;
    Ver4.Type[Ver4.DiagNum] = RENORMALIZED;
    SETOUTT(Ver4, Ver4.DiagNum, InT[INL], InT[INR]);
    Ver4.DiagNum++;

    Ver4.Channel[Ver4.DiagNum] = IEX;
    Ver4.Type[Ver4.DiagNum] = RENORMALIZED;
    SETOUTT(Ver4, Ver4.DiagNum, InT[INL], InT[INR]);
    Ver4.DiagNum++;
  }
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
  Ver4.DiagNum = 0;
  Ver4.Weight.fill(0.0);

  ASSERT_ALLWAYS(InT[RIGHT] - InT[LEFT] == 2 * (LoopNum + 1),
                 "Tau Number must be 2*(loopNum+1)");

  for (auto &chan : Channel) {
    for (int ol = 0; ol < LoopNum - 1; ol++) {
      int LTauNum = 2 * (ol + 1); //-2 because left and right InT are known
      array<int, 2> LInT = {InT[LEFT], InT[LEFT] + LTauNum};

      if (ol == 0)
        Ver4.LVer.push_back(ChanI(LInT, ol, BARE));
      else
        Ver4.LVer.push_back(ChanI(LInT, ol, RENORMALIZED));

      if (chan == U || chan == T)
        Ver4.LVer.push_back(Bubble(LInT, ol, {U, S}, RENORMALIZED));
      else
        Ver4.LVer.push_back(Bubble(LInT, ol, {U, T}, RENORMALIZED));

      int oR = LoopNum - 1 - ol;
      array<int, 2> RInT = {InT[LEFT] + LTauNum + 1, InT[RIGHT]};

      if (oR == 0)
        Ver4.RVer.push_back(ChanI(RInT, oR, BARE));
      else
        Ver4.RVer.push_back(ChanI(RInT, oR, RENORMALIZED));

      Ver4.RVer.push_back(Bubble(RInT, oR, {U, S, T}, RENORMALIZED));
    }
  }
}