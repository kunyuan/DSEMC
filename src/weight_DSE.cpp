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

    Vertex4(Root);

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
  momentum ExQ = InL - Var.LoopMom[K[OUTR]];
  Ver4.Weight[0] = VerQTheta.Interaction(InL, InR, DiQ, 0.0, 0) -
                   VerQTheta.Interaction(InL, InR, ExQ, 0.0, 0);
  // Ver4.Weight[0] = VerQTheta.Interaction(InL, InR, DiQ, 0.0, 0);
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
  const momentum &K1 = Var.LoopMom[Ver4.K[0]];
  int InTL = Ver4.T[0][INL];
  gMatrix &G1 = Ver4.G[0];
  double Weight;
  // set all weight element to be zero
  for (auto &w : Ver4.Weight)
    w = 0.0;

  for (auto &chan : Ver4.Channel) {
    // construct internal momentum
    momentum &K2 = Var.LoopMom[Ver4.K[chan]];
    gMatrix &G2 = Ver4.G[chan];
    if (chan == T) {
      K2 = OutL + K1 - InL;
    } else if (chan == U) {
      K2 = OutR + K1 - InL;
    } else if (chan == S) {
      K2 = InL + InR - K1;
    } else if (chan == I) {
      // TODO: add envelope diagram
      continue;
    }

    // construct Green's function weight table
    for (int lt = InTL; lt < InTL + Ver4.TauNum - 2; ++lt)
      for (int rt = InTL + 2; rt < InTL + Ver4.TauNum; ++rt) {
        double dTau = Var.Tau[rt] - Var.Tau[lt];
        G1(lt, rt) = Fermi.Green(dTau, K1, UP, 0, Var.CurrScale);
        if (chan == S)
          // LVer to RVer
          G2(lt, rt) = Fermi.Green(dTau, K2, UP, 0, Var.CurrScale);
        else
          // RVer to LVer
          G2(rt, lt) = Fermi.Green(-dTau, K2, UP, 0, Var.CurrScale);
      }
  }

  // for vertex4 with one or more loops
  for (auto &pair : Ver4.Pairs) {

    ver4 &LVer = pair.LVer;
    ver4 &RVer = pair.RVer;

    Vertex4(LVer);
    Vertex4(RVer);

    for (auto &m : pair.Map) {
      Weight = pair.SymFactor;
      Weight *= G1(m.G1T[IN], m.G1T[OUT]);
      Weight *= Ver4.G[pair.Chan](m.G2T[IN], m.G2T[OUT]);
      Weight *= LVer.Weight[m.LVerT] * RVer.Weight[m.RVerT];
      Ver4.Weight[m.T] += Weight;
    }
  }
}
