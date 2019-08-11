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

    double Weight = 0.0;
    if (LoopNum == 1) {
      Vertex4(Root);
    }
    // if (Order == 2)
    //   // cout << Weight << endl;
    // cout << count << endl;
    return Weight / pow(2.0 * PI, D * LoopNum);
    // return Weight;
  }
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

  for (auto &chan : Ver4.Channel) {
    if (chan == T) {
      Var.LoopMom[Ver4.K2t] = OutL + K1 - InL;
    } else if (chan == U) {
      Var.LoopMom[Ver4.K2u] = OutR + K1 - InL;
    } else if (chan == S) {
      Var.LoopMom[Ver4.K2u] = InL + InR - K1;
    }

    // for vertex4 with one or more loops
    for (auto &pair : Ver4.Pairs) {
      ver4 &LVer = pair.LVer;
      Bubble(LVer);

      ver4 &RVer = pair.RVer;
      Bubble(RVer);
    }
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

// int weight::Bubble(
//     const momentum &InL, const momentum &InR, const momentum &DirTran,
//     int Order, int TauIndex, int LoopIndex, int Level,
//     bool *Channel, // three flags, calculate t, u, s or not
//     int VerType,   // -1: normal, 0: left(to project), 1: right(to diff)
//     int LVerOrder, // order of left vertex
//     bool IsPenguin) {

//   // calculate renormalized diagrams
//   for (int OL = 0; OL < Order; OL++) {
//     if (LVerOrder >= 0 && OL != LVerOrder)
//       continue;
//     if (IsPenguin && OL < 1)
//       continue;

//     if (VerType == -1 || VerType == 1) {
//       // for normal vertex or projected vertex, just return
//       OneLoop(InL, InR, DirTran, Order, OL, TauIndex, LoopIndex, Level,
//       Channel,
//               false, // do not project
//               IsPenguin);
//       // if (Order == 2 && Level == 0) {
//       //   cout << VerType << ", index=" << _DiagIndex[Level]
//       //        << ", level=" << Level << " OL:" << OL << endl;
//       // }
//     }

//     /////// comment this block to calculate original diagrams  //////
//     if (VerType == 0 || VerType == 1) {
//       // do projection
//       OneLoop(InL, InR, DirTran, Order, OL, TauIndex, LoopIndex, Level,
//       Channel,
//               true, // do projection
//               IsPenguin);
//       // if (Order == 2 && Level == 0) {
//       //   cout << VerType << ", index=" << _DiagIndex[Level]
//       //        << ", level=" << Level << " OL:" << OL << endl;
//       // }
//     }
//   }
//   return _DiagIndex[Level];
// }

// int weight::Ver4Loop0(const momentum &InL, const momentum &InR,
//                       const momentum &DirTran, int TauIndex, int LoopIndex,
//                       int Level, int Type) {
//   //   return VerQTheta.Interaction(InL, InR, DirTran, 0, Var.CurrScale) -
//   //   VerQTheta.Interaction(InL, InR, ExTran, 0, Var.CurrScale);
//   double Tau = Var.Tau[TauIndex] - Var.Tau[TauIndex + 1];
//   int &Index = _DiagIndex[Level];
//   ////////////// bare interaction ///////////
//   double DiWeight = VerQTheta.Interaction(InL, InR, DirTran, Tau, 0);
//   double ExWeight =
//       VerQTheta.Interaction(InL, InR, InR + DirTran - InL, Tau, 0);
//   SETTAU(_GlobalOrder, Level, Index, TauIndex, TauIndex, TauIndex,
//   TauIndex); _Weight[Level][Index][0] = DiWeight - ExWeight;
//   // _Weight[Level][DiagIndex][0] = DiWeight;
//   Index += 1;

//   if (Type != -2) {
//     //////////// dressed interaction ///////////
//     DiWeight = VerQTheta.Interaction(InL, InR, DirTran, Tau, 1);
//     SETTAU(_GlobalOrder, Level, Index, TauIndex, TauIndex, TauIndex + 1,
//            TauIndex + 1);
//     _Weight[Level][Index][0] = DiWeight;
//     Index += 1;

//     ExWeight = VerQTheta.Interaction(InL, InR, InR + DirTran - InL, Tau,
//     1); SETTAU(_GlobalOrder, Level, Index, TauIndex, TauIndex + 1, TauIndex
//     + 1,
//            TauIndex);
//     _Weight[Level][Index][0] = -ExWeight;
//     Index += 1;
//   }

//   return Index;
// }

// int weight::OneLoop(const momentum &InL, const momentum &InR,
//                     const momentum &DirTran, int Order, int LVerOrder,
//                     int TauIndex, int LoopIndex, int Level,
//                     bool *Channel, // three flags, calculate t, u, s or not
//                     bool IsProjected, bool IsPenguin) {

//   double VerWeight;
//   double GR2L, GL2R;
//   double SymFactor = 1.0;
//   int nLevel = Level + 1;
//   int &Index = _DiagIndex[Level];

//   if (Order < 1)
//     return Index;

//   if (IsPenguin && (LVerOrder < 1 || Order < 2))
//     return Index;

//   // do projection
//   double ProjSign = 1.0;
//   if (IsProjected)
//     ProjSign = -1.0; // projection always comes with a minus sign

//   momentum Internal = Var.LoopMom[LoopIndex + LVerOrder];
//   momentum Internal2, VerLInL, VerLInR, VerLDiTran, VerRInL, VerRInR,
//       VerRDiTran;

//   int LTauIndex = TauIndex;
//   int RTauIndex = TauIndex + (LVerOrder + 1) * 2;

//   for (int chan = 0; chan < 3; chan++) {
//     if (!Channel[chan])
//       continue;
//     if (chan == 0) {
//       // t diagram
//       if (IsProjected) {
//         VerLInL = InL * (Para.Kf / InL.norm());
//         VerRInR = InR * (Para.Kf / InR.norm());
//       } else {
//         VerLInL = InL;
//         VerRInR = InR;
//       }
//       Internal2 = Internal - DirTran;
//       VerLInR = Internal2;
//       VerLDiTran = DirTran;

//       VerRInL = Internal;
//       VerRDiTran = DirTran;

//       SymFactor = -1.0;
//     } else if (chan == 1) {
//       // u diagram
//       if (IsProjected) {
//         VerLInL = InL * (Para.Kf / InL.norm());
//         VerRInR = InR * (Para.Kf / InR.norm());
//       } else {
//         VerLInL = InL;
//         VerRInR = InR;
//       }
//       Internal2 = Internal + DirTran - InL + InR;
//       // after the projection, the exchange transfer momentum
//       // should remain the same
//       VerLInR = Internal2;
//       VerLDiTran = Internal - Internal2;
//       // after the projection, the direct transfer momentum on left and right
//       // vertex also remain the same !!!

//       VerRInL = Internal;
//       VerRDiTran = VerLDiTran;
//       SymFactor = 1.0;
//     } else if (chan == 2) {
//       // projection is non-zero only for t and u channel
//       if (IsProjected)
//         continue;
//       // if (IsProjected) {
//       //   VerLInL = InL * (Para.Kf / InL.norm());
//       //   VerRInR = InR * (Para.Kf / InR.norm());
//       // } else {
//       //   VerLInL = InL;
//       //   VerRInR = InR;
//       // }

//       // s diagram
//       Internal2 = InL + InR - Internal;
//       VerLInL = InL;
//       VerLInR = InR;
//       VerLDiTran = VerLInL - Internal2;

//       VerRInL = Internal2;
//       VerRInR = Internal;
//       VerRDiTran = DirTran + InR - Internal;
//       SymFactor = 0.5;
//     }

//     //====================  DIRECT  Diagram =============================
//     // left vertex
//     bool *nChannel = ALL;
//     int LVerType = LEFT;
//     int LIndex = _DiagIndex[nLevel];

//     if (IsPenguin) {
//       LVerType = -1; // normal left vertex;
//       if (chan == 0 || chan == 1)
//         nChannel = US;
//       else
//         nChannel = UT;
//       Bubble(VerLInL, VerLInR, VerLDiTran, LVerOrder, LTauIndex, LoopIndex,
//              nLevel, nChannel,
//              LVerType, // VerType
//              -1,       // LVerOrder
//              false     // not penguin
//       );
//     } else {
//       // nChannel = T;
//       // if (Level == 0)
//       //   LVerType = -2;
//       Vertex4(VerLInL, VerLInR, VerLDiTran, LVerOrder, LTauIndex, LoopIndex,
//               nLevel, nChannel,
//               LVerType, // VerType
//               -1        // LVerOrder
//       );
//     }
//     int LDiagIndex = _DiagIndex[nLevel];

//     // nChannel = T;
//     // right vertex
//     int RIndex = LDiagIndex;
//     Vertex4(VerRInL, VerRInR, VerRDiTran, Order - 1 - LVerOrder, RTauIndex,
//             LoopIndex + 1 + LVerOrder, nLevel, ALL,
//             // nChannel,
//             RIGHT, // VerType
//             // -1, // VerType
//             -1 // LVerOrder
//     );

//     // if (LVerLoopNum == 0 && LoopNum == 2) {
//     //   nChannel = U;
//     //   // right vertex
//     //   nDiagIndex =
//     //       Vertex4(VerRInL, VerRInR, VerRDiTran, LoopNum - 1 - LVerLoopNum,
//     //               RTauIndex, LoopIndex + 1 + LVerLoopNum, nDiagIndex,
//     //               nLevel,
//     //               // ALL,
//     //               nChannel,
//     //               RIGHT, // VerType
//     //               -1     // LVerOrder
//     //       );
//     // }
//     int RDiagIndex = _DiagIndex[nLevel];

//     ////////////// construct  G table  /////////////////////////////////////
//     // for (int tL = LTauIndex; tL < RTauIndex; tL++)
//     //   for (int tR = RTauIndex; tR < TauIndex + (Order + 1) * 2; tR++) {
//     //     _GL2R[tL][tR] = Fermi.Green(Var.Tau[tR] - Var.Tau[tL], Internal,
//     //     UP, 0,
//     //                                 Var.CurrScale);
//     //     if (chan == 2)
//     //       _GR2L[tL][tR] = Fermi.Green(Var.Tau[tR] - Var.Tau[tL],
//     Internal2,
//     //       UP,
//     //                                   0, Var.CurrScale);
//     //     else
//     //       _GR2L[tR][tL] = Fermi.Green(Var.Tau[tL] - Var.Tau[tR],
//     Internal2,
//     //       UP,
//     //                                   0, Var.CurrScale);
//     //     // cout << tR << "-" << tL << ", " << _GR2L[tR][tL] << endl;
//     //   }

//     if (chan == 2) {
//       for (int tL = LTauIndex; tL < RTauIndex; tL++) {
//         _GR2L[tL][RTauIndex] = Fermi.Green(Var.Tau[RTauIndex] - Var.Tau[tL],
//                                            Internal2, UP, 0, Var.CurrScale);
//         for (int tR = RTauIndex; tR < TauIndex + (Order + 1) * 2; tR++)
//           _GL2R[tL][tR] = Fermi.Green(Var.Tau[tR] - Var.Tau[tL], Internal,
//           UP,
//                                       0, Var.CurrScale);
//       }
//     } else {
//       for (int tL = LTauIndex; tL < RTauIndex; tL++) {
//         _GL2R[tL][RTauIndex] = Fermi.Green(Var.Tau[RTauIndex] - Var.Tau[tL],
//                                            Internal, UP, 0, Var.CurrScale);
//         for (int tR = RTauIndex; tR < TauIndex + (Order + 1) * 2; tR++)
//           _GR2L[tR][tL] = Fermi.Green(Var.Tau[tL] - Var.Tau[tR], Internal2,
//           UP,
//                                       0, Var.CurrScale);
//       }
//     }

//     int *_ExtTauC = _ExtTau[_GlobalOrder][Level][Index];
//     for (int l = LIndex; l < LDiagIndex; l++) {
//       int *_ExtTauL = _ExtTau[_GlobalOrder][nLevel][l];
//       for (int r = RIndex; r < RDiagIndex; r++) {
//         int *_ExtTauR = _ExtTau[_GlobalOrder][nLevel][r];

//         if (chan == 0) {
//           SETTAU(_GlobalOrder, Level, Index, _ExtTauL[INL], _ExtTauL[OUTL],
//                  _ExtTauR[INR], _ExtTauR[OUTR]);

//           GR2L = _GR2L[_ExtTauR[OUTL]][_ExtTauL[INR]];
//           GL2R = _GL2R[_ExtTauL[OUTR]][_ExtTauR[INL]];
//           // cout << "0me: " << GR2L << ", " << GL2R << endl;
//         } else if (chan == 1) {
//           SETTAU(_GlobalOrder, Level, Index, _ExtTauL[INL], _ExtTauR[OUTR],
//                  _ExtTauR[INR], _ExtTauL[OUTL]);

//           GR2L = _GR2L[_ExtTauR[OUTL]][_ExtTauL[INR]];
//           GL2R = _GL2R[_ExtTauL[OUTR]][_ExtTauR[INL]];
//           // cout << "1me: " << GR2L << ", " << GL2R << endl;
//         } else if (chan == 2) {
//           SETTAU(_GlobalOrder, Level, Index, _ExtTauL[INL], _ExtTauR[OUTL],
//                  _ExtTauL[INR], _ExtTauR[OUTR]);

//           GR2L = _GR2L[_ExtTauL[OUTL]][_ExtTauR[INL]];
//           GL2R = _GL2R[_ExtTauL[OUTR]][_ExtTauR[INR]];
//           // cout << "2me: " << GR2L << ", " << GL2R << endl;
//         }

//         if (IsProjected)
//           if (chan == 0) {
//             SETTAU(_GlobalOrder, Level, Index, _ExtTauC[INL], _ExtTauC[INL],
//                    _ExtTauC[INR], _ExtTauC[INR]);
//           } else if (chan == 1) {
//             SETTAU(_GlobalOrder, Level, Index, _ExtTauC[INL], _ExtTauC[INR],
//                    _ExtTauC[INR], _ExtTauC[INL]);
//           }

//         _Weight[Level][Index][0] = _Weight[nLevel][l][0] *
//                                    _Weight[nLevel][r][0] * GL2R * GR2L *
//                                    SymFactor * ProjSign;

//         // if (_GlobalOrder == 2 && l == LIndex && r == RIndex) {
//         //   cout << "Lindex: " << LIndex << " Weight:" <<
//         //   _Weight[Level][Index][0]
//         //        << endl;
//         //   cout << Var.Tau[_ExtTauL[OUTR]] << ", " <<
//         Var.Tau[_ExtTauR[INL]]
//         //        << ", Mom=" << Internal.norm() << ": GL2R=" << GL2R <<
//         endl;
//         //   // cout << Fermi.Green(Var.Tau[_ExtTauR[INL]] -
//         //   // Var.Tau[_ExtTauL[OUTR]],
//         //   //                     Internal, UP, 0, Var.CurrScale)
//         //   //      << endl;
//         //   cout << _ExtTauR[OUTL] << ", " << _ExtTauL[INR]
//         //        << ", Mom=" << Internal2.norm() << ": GR2L=" << GR2L
//         //        << ", chan=" << chan << endl;
//         //   cout << Fermi.Green(Var.Tau[_ExtTauL[INR]] -
//         //   Var.Tau[_ExtTauR[OUTL]],
//         //                       Internal2, UP, 0, Var.CurrScale)
//         //        << endl;
//         //   cout << _Weight[nLevel][l][0] << endl;
//         //   cout << _Weight[nLevel][r][0] << endl;
//         //   // cout << GL2R << " vs " << GR2L << endl;
//         //   cout << "end" << endl;
//         //   cout << endl;
//         // }

//         Index++;
//       }
//     }
//   }
//   return Index;
// }