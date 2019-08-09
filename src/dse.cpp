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

#define SETK(Ver4, LegK)                                                       \
  Ver4.LegK[INL] = LegK[INL];                                                  \
  Ver4.LegK[INR] = LegK[INR];                                                  \
  Ver4.LegK[OUTL] = LegK[OUTL];                                                \
  Ver4.LegK[OUTR] = LegK[OUTR];

#define SETINT(Ver4, InT)                                                      \
  Ver4.InT[INL] = InT[INL];                                                    \
  Ver4.InT[INR] = InT[INR];

#define SETOUTT(Ver4, DiagNum, OutL, OutR)                                     \
  Ver4.OutT[DiagNum][OUTL] = OutL;                                             \
  Ver4.OutT[DiagNum][OUTR] = OutR;

ver4 *verDiag::ChanI(array<int, 4> LegK, array<int, 2> InT, int LoopNum,
                     int TauIndex, int LoopIndex, vertype Type) {
  VerPool.push_back(ver4());
  ver4 &Ver4 = VerPool.back();
  Ver4.Channel = 0;
  if (LoopNum == 0) {
    Ver4.LoopNum = 0;
    ////////////// bare interaction ///////////
    Ver4.DiagNum = 0;
    Ver4.Type[Ver4.DiagNum] = BARE;
    SETK(Ver4, LegK);
    SETINT(Ver4, InT);
    SETOUTT(Ver4, Ver4.DiagNum, InT[INL], InT[INR]);
    Ver4.Weight[Ver4.DiagNum] = 0.0;
    Ver4.DiagNum++;

    if (Type == DYNAMIC) {
      //////////// dressed interaction ///////////
      Ver4.Type[Ver4.DiagNum] = DYNAMIC;
      SETOUTT(Ver4, Ver4.DiagNum, InT[INL], InT[INR]);
      Ver4.Weight[Ver4.DiagNum] = 0.0;
      Ver4.DiagNum++;

      Ver4.Type[Ver4.DiagNum] = DYNAMIC;
      SETOUTT(Ver4, Ver4.DiagNum, InT[INL], InT[INR]);
      Ver4.Weight[Ver4.DiagNum] = 0.0;
      Ver4.DiagNum++;
    }
    return &VerPool.back();
  } else {
    // TODO: add envolpe diagrams
    return nullptr;
  }
}

ver4 *verDiag::ChanT(array<int, 4> LegK, array<int, 2> InT, int LoopNum,
                     int TauIndex, int LoopIndex, vertype Type) {
  VerPool.push_back(ver4());
  ver4 &Ver4 = VerPool.back();
  Ver4.LoopNum = LoopNum;
  Ver4.DiagNum = 0;
}