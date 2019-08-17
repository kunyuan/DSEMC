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

double weight::Evaluate(int LoopNum, int Channel) {
  if (LoopNum == 0) {
    // normalization
    return VerQTheta.Interaction(Var.LoopMom[1], Var.LoopMom[2], Var.LoopMom[0],
                                 0.0, -2);
  } else {
    ver4 &Root = Ver4Root[LoopNum][Channel];
    if (Root.Weight.size() == 0)
      // empty vertex
      return 0.0;

    *Root.LegK[OUTL] = Var.LoopMom[1] - Var.LoopMom[0];
    *Root.LegK[OUTR] = Var.LoopMom[2] + Var.LoopMom[0];

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
  const momentum &InL = *K[INL];
  const momentum &InR = *K[INR];
  momentum DiQ = InL - *K[OUTL];
  momentum ExQ = InL - *K[OUTR];
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
  // cout << Ver4.LoopNum << endl;
  if (Ver4.LoopNum == 0) {
    Ver0(Ver4);
  } else {
    for (auto &w : Ver4.Weight)
      w = 0.0;
    ChanUST(Ver4);
    if (Ver4.LoopNum >= 3)
      ChanI(Ver4);
  }
  return;
}

void weight::ChanUST(dse::ver4 &Ver4) {
  const momentum &InL = *Ver4.LegK[INL];
  const momentum &OutL = *Ver4.LegK[OUTL];
  const momentum &InR = *Ver4.LegK[INR];
  const momentum &OutR = *Ver4.LegK[OUTR];
  const momentum &K1 = *Ver4.K[0];
  int InTL = Ver4.T[0][INL];
  gMatrix &G1 = Ver4.G[0];
  double Weight;
  // set all weight element to be zero
  for (auto &chan : Ver4.Channel) {
    // construct internal momentum
    momentum &K2 = *Ver4.K[chan];
    gMatrix &G2 = Ver4.G[chan];
    if (chan == T) {
      K2 = OutL + K1 - InL;
    } else if (chan == U) {
      K2 = OutR + K1 - InL;
    } else if (chan == S) {
      K2 = InL + InR - K1;
    } else if (chan == I) {
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

void weight::ChanI(dse::ver4 &Ver4) {
  if (Ver4.LoopNum != 3)
    return;
  if (Ver4.Channel[0] != I)
    return;
  envelope &Env = Ver4.Envelope[0];

  const momentum &InL = *Ver4.LegK[INL];
  const momentum &OutL = *Ver4.LegK[OUTL];
  const momentum &InR = *Ver4.LegK[INR];
  const momentum &OutR = *Ver4.LegK[OUTR];

  auto &G = Env.G;

  *G[3].K = *G[0].K + *G[1].K - InL;
  *G[4].K = *G[1].K + *G[2].K - OutL;
  *G[5].K = *G[0].K + InR - *G[2].K;
  *G[6].K = *G[1].K + *G[2].K - OutR;
  *G[7].K = *G[2].K + OutR - *G[0].K;
  *G[8].K = *G[2].K + OutL - *G[0].K;

  for (auto &g : Env.G)
    for (auto &in : g.InT)
      for (auto &out : g.OutT)
        g(in, out) = Fermi.Green(Var.Tau[out] - Var.Tau[in], *(g.K), UP, 0,
                                 Var.CurrScale);

  for (auto &subVer : Env.Ver)
    Vertex4(subVer);

  double Weight = 0.0;
  double ComWeight = 0.0;
  for (auto &map : Env.Map) {
    auto &SubVer = Env.Ver;
    auto &GT = map.GT;
    auto &G = Env.G;
    ComWeight = G[0](GT[0]) * G[1](GT[1]) * G[2](GT[2]) * G[3](GT[3]);
    ComWeight *= SubVer[0].Weight[map.LDVerT];

    Weight = Env.SymFactor[0] * ComWeight;
    Weight *= SubVer[1].Weight[map.LUVerT];
    Weight *= SubVer[3].Weight[map.RDVerT];
    Weight *= SubVer[6].Weight[map.RUVerT];
    Weight *= G[4](GT[4]) * G[5](GT[5]);
    Ver4.Weight[map.T[0]] = Weight;

    Weight = Env.SymFactor[1] * ComWeight;
    Weight *= SubVer[2].Weight[map.LUVerT];
    Weight *= SubVer[3].Weight[map.RDVerT];
    Weight *= SubVer[7].Weight[map.RUVerT];
    Weight *= G[6](GT[6]) * G[5](GT[5]);
    Ver4.Weight[map.T[1]] = Weight;

    Weight = Env.SymFactor[2] * ComWeight;
    Weight *= SubVer[1].Weight[map.LUVerT];
    Weight *= SubVer[4].Weight[map.RDVerT];
    Weight *= SubVer[8].Weight[map.RUVerT];
    Weight *= G[4](GT[4]) * G[7](GT[7]);
    Ver4.Weight[map.T[2]] = Weight;

    Weight = Env.SymFactor[3] * ComWeight;
    Weight *= SubVer[2].Weight[map.LUVerT];
    Weight *= SubVer[5].Weight[map.RDVerT];
    Weight *= SubVer[9].Weight[map.RUVerT];
    Weight *= G[6](GT[6]) * G[8](GT[8]);
    Ver4.Weight[map.T[3]] = Weight;
  }

  return;
}
