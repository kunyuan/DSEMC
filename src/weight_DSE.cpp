#include "global.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/vector.h"
#include "weight.h"
#include <array>
#include <iostream>
#include <stack>
#include <string>

using namespace diag;
using namespace std;
using namespace dse;

#define TIND(Shift, LTau, RTau) ((LTau - Shift) * MaxTauNum + RTau - Shift)

double weight::Evaluate(int LoopNum, int ID) {
  if (LoopNum == 0) {
    // normalization
    return VerQTheta.Interaction(Var.LoopMom[1], Var.LoopMom[2], Var.LoopMom[0],
                                 0.0, -2);
  } else {
    ver4 &Root = Ver4Root[LoopNum];
    Var.LoopMom[Root.LegK[OUTL]] = Var.LoopMom[1] - Var.LoopMom[0];
    Var.LoopMom[Root.LegK[OUTR]] = Var.LoopMom[2] + Var.LoopMom[0];

    if (LoopNum == 1) {
      Vertex4(Root);
    }

    double Weight = 0.0;
    for (auto &w : Root.Weight)
      Weight += w;
    // if (Order == 2)
    //   // cout << Weight << endl;
    // cout << count << endl;
    return Weight / pow(2.0 * PI, D * LoopNum);
    // return Weight;
  }
}

void weight::Ver0(ver4 &Ver4) {
  auto &K = Ver4.LegK;
  const momentum &InL = Var.LoopMom[K[INL]];
  const momentum &InR = Var.LoopMom[K[INR]];
  momentum DiQ = InL - Var.LoopMom[K[OUTL]];
  momentum ExQ = InR - Var.LoopMom[K[OUTR]];
  Ver4.Weight[0] = VerQTheta.Interaction(InL, InR, DiQ, 0.0, 0) -
                   VerQTheta.Interaction(InL, InR, ExQ, 0.0, 0);

  if (Ver4.Type != dse::BARE) {
    double Tau = Var.Tau[Ver4.T[0][INR]] - Var.Tau[Ver4.T[0][INL]];
    Ver4.Weight[1] = VerQTheta.Interaction(InL, InR, DiQ, Tau, 1);
    Ver4.Weight[2] = -VerQTheta.Interaction(InL, InR, ExQ, Tau, 1);
  }
  return;
}

void weight::Vertex4(dse::ver4 &Ver4) {
  if (Ver4.LoopNum == 0) {
    Ver0(Ver4);
    return;
  }
  const momentum &InL = Var.LoopMom[Ver4.LegK[INL]];
  const momentum &OutL = Var.LoopMom[Ver4.LegK[OUTL]];
  const momentum &InR = Var.LoopMom[Ver4.LegK[INR]];
  const momentum &OutR = Var.LoopMom[Ver4.LegK[OUTR]];
  const momentum &K1 = Var.LoopMom[Ver4.K1];
  int InTL = Ver4.T[0][INL];
  // set all weight element to be zero
  for (auto &w : Ver4.Weight)
    w = 0.0;

  for (auto &chan : Ver4.Channel) {
    // construct internal momentum
    if (chan == T) {
      Var.LoopMom[Ver4.K2[chan]] = OutL + K1 - InL;
    } else if (chan == U) {
      Var.LoopMom[Ver4.K2[chan]] = OutR + K1 - InL;
    } else if (chan == S) {
      Var.LoopMom[Ver4.K2[chan]] = InL + InR - K1;
    }

    // construct Green's function weight
    const momentum &K2 = Var.LoopMom[Ver4.K2[chan]];
    for (int lt = InTL; lt < InTL + Ver4.TauNum - 2; ++lt)
      for (int rt = InTL + 2; rt < InTL + Ver4.TauNum; ++lt) {
        double dTau = Var.Tau[rt] - Var.Tau[lt];
        Ver4.G1(lt, rt) = Fermi.Green(dTau, K1, UP, 0, Var.CurrScale);
        Ver4.G2[chan](lt, rt) = Fermi.Green(dTau, K2, UP, 0, Var.CurrScale);
      }
  }

  // for vertex4 with one or more loops
  for (auto &pair : Ver4.Pairs) {

    ver4 &LVer = pair.LVer;
    ver4 &RVer = pair.RVer;

    Vertex4(LVer);
    Vertex4(RVer);

    for (int l = 0; l < LVer.T.size(); ++l)
      for (int r = 0; r < RVer.T.size(); ++r) {

        double Weight = pair.SymFactor;
        auto &Int1 = pair.IntT1(l, r);
        auto &Int2 = pair.IntT2(l, r);

        Weight *= Ver4.G1(Int1[IN], Int1[OUT]);
        Weight *= Ver4.G2[pair.Chan](Int2[IN], Int2[OUT]);
        Weight *= LVer.Weight[l] * RVer.Weight[r];

        Ver4.Weight[pair.Map(l, r)] += Weight;
      }
  }
}
