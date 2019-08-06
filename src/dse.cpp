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

dse::ver4 *verDiag::Ver0(array<int, 4> LegK, array<int, 2> InLegT,
                         bool IsBare) {

  //   return VerQTheta.Interaction(InL, InR, DirTran, 0, Var.CurrScale) -
  //   VerQTheta.Interaction(InL, InR, ExTran, 0, Var.CurrScale);
  double Tau = Var.Tau[TauIndex] - Var.Tau[TauIndex + 1];
  int &Index = _DiagIndex[Level];
  ////////////// bare interaction ///////////
  double DiWeight = VerQTheta.Interaction(InL, InR, DirTran, Tau, 0);
  double ExWeight =
      VerQTheta.Interaction(InL, InR, InR + DirTran - InL, Tau, 0);
  SETTAU(_GlobalOrder, Level, Index, TauIndex, TauIndex, TauIndex, TauIndex);
  _Weight[Level][Index][0] = DiWeight - ExWeight;
  // _Weight[Level][DiagIndex][0] = DiWeight;
  Index += 1;

  if (Type != -2) {
    //////////// dressed interaction ///////////
    DiWeight = VerQTheta.Interaction(InL, InR, DirTran, Tau, 1);
    SETTAU(_GlobalOrder, Level, Index, TauIndex, TauIndex, TauIndex + 1,
           TauIndex + 1);
    _Weight[Level][Index][0] = DiWeight;
    Index += 1;

    ExWeight = VerQTheta.Interaction(InL, InR, InR + DirTran - InL, Tau, 1);
    SETTAU(_GlobalOrder, Level, Index, TauIndex, TauIndex + 1, TauIndex + 1,
           TauIndex);
    _Weight[Level][Index][0] = -ExWeight;
    Index += 1;
  }
  return;
}