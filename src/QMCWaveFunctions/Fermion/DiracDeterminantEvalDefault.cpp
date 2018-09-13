#include "QMCWaveFunctions/Fermion/DiracDeterminantEval.h"
#include "Configuration.h"

namespace qmcplusplus
{
using QMCT = QMCTraits;
  
// QMCT::ValueType DiracDeterminantEvalDefault::ratio(ParticleSet& P, int iat)
// {
//   UpdateMode=ORB_PBYP_RATIO;
//   WorkingIndex = iat-FirstIndex;
//   SPOVTimer.start();
//   Phi->evaluate(P, iat, psiV);
//   SPOVTimer.stop();
//   RatioTimer.start();
//   curRatio=simd::dot(psiM[WorkingIndex],psiV.data(),NumOrbitals);
//   //curRatio = DetRatioByRow(psiM, psiV,WorkingIndex);
//   RatioTimer.stop();
//   return curRatio;
// }

}
