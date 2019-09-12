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
    // return VerQTheta.Interaction(Var.LoopMom[1], Var.LoopMom[2],
    // Var.LoopMom[0],
    //                              0.0, -2);
    return 1.0;
  } else {
    // if (Channel != dse::T)
    //   return 0.0;

    ver4 &Root = Ver4Root[LoopNum][Channel];
    if (Root.Weight.size() == 0)
      // empty vertex
      return 0.0;

    if (Channel == dse::S) {
      *Root.LegK[INR] = Var.LoopMom[0] - Var.LoopMom[1];
      *Root.LegK[OUTR] = Var.LoopMom[0] - Var.LoopMom[2];
    } else {
      *Root.LegK[OUTL] = Var.LoopMom[1] - Var.LoopMom[0];
      *Root.LegK[OUTR] = Var.LoopMom[2] + Var.LoopMom[0];
    }

    Vertex4(Root);

    double Weight = 0.0;
    for (auto &w : Root.Weight)
      Weight += w;
    // if (LoopNum == 3 && Channel == dse::I) {
    //   cout << "loopnum: " << Root.LoopNum << endl;
    //   cout << "channel: " << Root.Channel[0] << endl;
    //   cout << Weight << endl;
    // }
    // cout << count << endl;
    return Weight / pow(2.0 * PI, D * LoopNum);
    // return Weight;
  }
}

void weight::Ver0(ver4 &Ver4) {
  array<momentum *, 4> &K = Ver4.LegK;
  momentum DiQ = *K[INL] - *K[OUTL];
  momentum ExQ = *K[INL] - *K[OUTR];
  Ver4.Weight[0] = VerQTheta.Interaction(K, DiQ, 0.0, 0) -
                   VerQTheta.Interaction(K, ExQ, 0.0, 0);
  // Ver4.Weight[0] = 1.0 / Para.Beta;
  if (Ver4.RexpandBare) {
    // cout << Ver4.T[0][INR] << ", " << Ver4.T[0][INL] << endl;
    // double Tau = Var.Tau[Ver4.T[1][INR]] - Var.Tau[Ver4.T[1][INL]];
    // cout << Ver4.T[1][INR] << ", " << Ver4.T[1][INL] << "; " <<
    // Ver4.T[2][INR]
    //      << ", " << Ver4.T[2][INL] << endl;
    Ver4.Weight[0] = +VerQTheta.Interaction(K, DiQ, 0.0, 1) -
                     VerQTheta.Interaction(K, ExQ, 0.0, 1);
    // Ver4.Weight[1] = 0.0;
    // Ver4.Weight[2] = 0.0;

    // Ver4.Weight[1] = +VerQTheta.Interaction(K, DiQ, Tau, 1);
    // Ver4.Weight[2] = -VerQTheta.Interaction(K, ExQ, Tau, 1);
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
  double Weight = 0.0;
  double Ratio;
  momentum Transfer;
  array<momentum *, 4> &LegK0 = Ver4.LegK;
  // if (Para.Counter == 181406 && Ver4.ID == 2 && Ver4.LoopNum == 1) {
  //   cout << "here" << endl;
  // }

  for (auto &bubble : Ver4.Bubble) {
    auto &G = bubble.G;
    const momentum &K0 = *G[0].K;
    int InTL = bubble.InTL;

    for (auto &chan : bubble.Channel) {
      bubble.ProjFactor[chan] = 1;
      array<momentum *, 4> &LegK = bubble.LegK[chan];

      if (bubble.IsProjected) {
        if (chan == S) {
          continue;
        } else if (chan == T) {
          Transfer = *LegK0[INL] - *LegK0[OUTL];
          double Q = Transfer.norm();
          if (Q < 0.05 * Para.Kf) {
            Ratio = Para.Kf / (*LegK0[INL]).norm();
            *LegK[INL] = *LegK0[INL] * Ratio;
            Ratio = Para.Kf / (*LegK0[INR]).norm();
            *LegK[INR] = *LegK0[INR] * Ratio;

            *LegK[OUTL] = *LegK[INL] - Transfer;
            *LegK[OUTR] = *LegK[INR] + Transfer;
            // } else if (Q > 1.8 * Para.Kf && Q < 2.2 * Para.Kf) {
            //   if (((*LegK0[INL]).norm() > 0.8 * Para.Kf &&
            //        (*LegK0[INL]).norm() < 1.2 * Para.Kf) &&
            //       ((*LegK0[INR]).norm() > 0.8 * Para.Kf &&
            //        (*LegK0[INR]).norm() < 1.2 * Para.Kf) &&
            //       ((*LegK0[OUTL]).norm() > 0.8 * Para.Kf &&
            //        (*LegK0[OUTL]).norm() < 1.2 * Para.Kf) &&
            //       ((*LegK0[OUTR]).norm() > 0.8 * Para.Kf &&
            //        (*LegK0[OUTR]).norm() < 1.2 * Para.Kf)) {

            //     Ratio = 2.0 * Para.Kf / Q;
            //     Transfer = Transfer * Ratio;
            //     *LegK[INL] = Transfer * 0.5;
            //     *LegK[INR] = Transfer * (-0.5);
            //     *LegK[OUTL] = *LegK[INR];
            //     *LegK[OUTR] = *LegK[INL];
            //   } else {
            //     bubble.ProjFactor[T] = 0.0;
            //   }
          } else {
            bubble.ProjFactor[T] = 0.0;
            continue;
          }
        } else {
          Transfer = *LegK0[INL] - *LegK0[OUTR];
          double Q = Transfer.norm();
          if (Q < 0.05 * Para.Kf) {
            Ratio = Para.Kf / (*LegK0[INL]).norm();
            *LegK[INL] = *LegK0[INL] * Ratio;
            Ratio = Para.Kf / (*LegK0[INR]).norm();
            *LegK[INR] = *LegK0[INR] * Ratio;

            *LegK[OUTL] = *LegK[INR] + Transfer;
            *LegK[OUTR] = *LegK[INL] - Transfer;
            // } else if (Q > 1.8 * Para.Kf && Q < 2.2 * Para.Kf) {

            //   if (((*LegK0[INL]).norm() > 0.8 * Para.Kf &&
            //        (*LegK0[INL]).norm() < 1.2 * Para.Kf) &&
            //       ((*LegK0[INR]).norm() > 0.8 * Para.Kf &&
            //        (*LegK0[INR]).norm() < 1.2 * Para.Kf) &&
            //       ((*LegK0[OUTL]).norm() > 0.8 * Para.Kf &&
            //        (*LegK0[OUTL]).norm() < 1.2 * Para.Kf) &&
            //       ((*LegK0[OUTR]).norm() > 0.8 * Para.Kf &&
            //        (*LegK0[OUTR]).norm() < 1.2 * Para.Kf)) {

            //     Ratio = 2.0 * Para.Kf / Q;
            //     Transfer = Transfer * Ratio;
            //     *LegK[INL] = Transfer * 0.5;
            //     *LegK[INR] = Transfer * (-0.5);
            //     *LegK[OUTL] = *LegK[INL];
            //     *LegK[OUTR] = *LegK[INR];
            //   } else {
            //     bubble.ProjFactor[U] = 0.0;
            //   }
          } else {
            bubble.ProjFactor[U] = 0.0;
            continue;
          }
        }
      }

      if (chan == T)
        *G[T].K = *LegK[OUTL] + K0 - *LegK[INL];
      else if (chan == U)
        *G[U].K = *LegK[OUTR] + K0 - *LegK[INL];
      else if (chan == S)
        *G[S].K = *LegK[INL] + *LegK[INR] - K0;
    }

    for (int lt = InTL; lt < InTL + Ver4.TauNum - 1; ++lt)
      for (int rt = InTL + 1; rt < InTL + Ver4.TauNum; ++rt) {
        double dTau = Var.Tau[rt] - Var.Tau[lt];
        G[0](lt, rt) = Fermi.Green(dTau, K0, UP, 0, Var.CurrScale);
        for (auto &chan : bubble.Channel) {
          if (abs(bubble.ProjFactor[chan]) > EPS)
            if (chan == S)
              // LVer to RVer
              G[S](lt, rt) = Fermi.Green(dTau, *G[S].K, UP, 0, Var.CurrScale);
            else
              // RVer to LVer
              G[chan](rt, lt) =
                  Fermi.Green(-dTau, *G[chan].K, UP, 0, Var.CurrScale);
        }
      }

    // for vertex4 with one or more loops
    for (auto &pair : bubble.Pair) {
      if (abs(bubble.ProjFactor[pair.Channel]) < EPS)
        continue;
      ver4 &LVer = pair.LVer;
      ver4 &RVer = pair.RVer;
      Vertex4(LVer);
      Vertex4(RVer);

      for (auto &map : pair.Map) {
        Weight = pair.SymFactor * bubble.ProjFactor[pair.Channel];
        Weight *= G[0](map.G0T) * G[pair.Channel](map.GT);
        Weight *= LVer.Weight[map.LVerTidx] * RVer.Weight[map.RVerTidx];
        Ver4.Weight[map.Tidx] += Weight;
      }
    }
  }
}

void weight::ChanI(dse::ver4 &Ver4) {
  if (Ver4.LoopNum != 3)
    return;
  for (auto &Env : Ver4.Envelope) {
    const momentum &InL = *Env.LegK[INL];
    const momentum &OutL = *Env.LegK[OUTL];
    const momentum &InR = *Env.LegK[INR];
    const momentum &OutR = *Env.LegK[OUTR];

    auto &G = Env.G;

    *G[3].K = *G[0].K + *G[1].K - InL;
    *G[4].K = *G[1].K + *G[2].K - OutL;
    *G[5].K = *G[0].K + InR - *G[2].K;
    *G[6].K = *G[1].K + *G[2].K - OutR;
    *G[7].K = *G[2].K + OutR - *G[0].K;
    *G[8].K = *G[2].K + OutL - *G[0].K;

    for (auto &g : Env.G)
      g.Weight = Fermi.Green(Var.Tau[g.OutT] - Var.Tau[g.InT], *(g.K), UP, 0,
                             Var.CurrScale);

    for (auto &subVer : Env.Ver)
      Vertex4(subVer);

    double Weight = 0.0;
    double ComWeight = 0.0;
    for (auto &map : Env.Map) {
      auto &SubVer = Env.Ver;
      auto &GT = map.GT;
      auto &G = Env.G;
      ComWeight = G[0].Weight * G[1].Weight * G[2].Weight * G[3].Weight;
      // cout << "G: " << ComWeight << endl;
      ComWeight *= SubVer[0].Weight[map.LDVerTidx];
      // cout << "Ver: " << SubVer[0].Weight[map.LDVerT] << endl;
      // cout << "T: " << map.LDVerT << endl;

      Weight = Env.SymFactor[0] * ComWeight;
      Weight *= SubVer[1].Weight[map.LUVerTidx];
      Weight *= SubVer[3].Weight[map.RDVerTidx];
      Weight *= SubVer[6].Weight[map.RUVerTidx];
      Weight *= G[4].Weight * G[5].Weight;
      Ver4.Weight[map.Tidx[0]] += Weight;

      Weight = Env.SymFactor[1] * ComWeight;
      Weight *= SubVer[2].Weight[map.LUVerTidx];
      Weight *= SubVer[3].Weight[map.RDVerTidx];
      Weight *= SubVer[7].Weight[map.RUVerTidx];
      Weight *= G[6].Weight * G[5].Weight;
      Ver4.Weight[map.Tidx[1]] += Weight;
      // cout << Weight << endl;

      Weight = Env.SymFactor[2] * ComWeight;
      Weight *= SubVer[1].Weight[map.LUVerTidx];
      Weight *= SubVer[4].Weight[map.RDVerTidx];
      Weight *= SubVer[8].Weight[map.RUVerTidx];
      Weight *= G[4].Weight * G[7].Weight;
      Ver4.Weight[map.Tidx[2]] += Weight;
      // cout << Weight << endl;

      Weight = Env.SymFactor[3] * ComWeight;
      Weight *= SubVer[2].Weight[map.LUVerTidx];
      Weight *= SubVer[5].Weight[map.RDVerTidx];
      Weight *= SubVer[9].Weight[map.RUVerTidx];
      Weight *= G[6].Weight * G[8].Weight;
      Ver4.Weight[map.Tidx[3]] += Weight;
      // cout << Weight << endl;

      // if (map.LDVerT == 0 && map.LUVerT == 0 && map.RDVerT == 0 &&
      //     map.RUVerT == 0) {
      // cout << "Com: " << ComWeight << endl;
      // cout << "G[4]: " << G[4](GT[4]) << endl;
      // cout << "G[5]: " << G[5](GT[5]) << endl;
      // cout << SubVer[1].Weight[map.LUVerT] << endl;
      // cout << SubVer[3].Weight[map.RDVerT] << endl;
      // cout << SubVer[6].Weight[map.RUVerT] << endl;
      // cout << "First: " << Weight << endl;
      // }
    }
  }

  return;
}
