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

ver4 *verDiag::Ver0(array<int, 4> LegK, array<int, 2> InT, bool IsBare) {
  ver4 &Ver4 = Pool.Ver4[Pool.Ver4Index];
  //   struct ver4 {
  //     int Level;
  //     int DiagNum;
  //     int Type;
  //     bool IsProjection;
  //     ver4 *SubVer[2];
  //     momentum *LegK[4];
  //     int LegT[4];
  //     green *Internal[MaxOrder * 2];
  //     double Weight[MaxOrder * MaxOrder * 4];
  //   };
  Ver4.LoopNum = 0;
  ////////////// bare interaction ///////////
  Ver4.DiagNum = 1;
  Ver4.Type[Ver4.DiagNum] = BARE;
  SETK(Ver4, LegK);
  SETINT(Ver4, InT);
  SETOUTT(Ver4, Ver4.DiagNum, InT[INL], InT[INR]);

  if (IsBare == false) {
    //////////// dressed interaction ///////////
    Ver4.DiagNum++;
    SETOUTT(Ver4, Ver4.DiagNum, InT[INL], InT[INR]);
    //     DiWeight = VerQTheta.Interaction(InL, InR, DirTran, Tau, 1);
    //     SETTAU(_GlobalOrder, Level, Index, TauIndex, TauIndex, TauIndex + 1,
    //            TauIndex + 1);
    //     _Weight[Level][Index][0] = DiWeight;
    //     Index += 1;

    //     ExWeight = VerQTheta.Interaction(InL, InR, InR + DirTran - InL, Tau,
    //     1); SETTAU(_GlobalOrder, Level, Index, TauIndex, TauIndex + 1,
    //     TauIndex
    //     + 1,
    //            TauIndex);
    //     _Weight[Level][Index][0] = -ExWeight;
    //     Index += 1;
  }
  return &Ver4;
}