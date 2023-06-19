//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 Jeongnim Kim and QMCPACK developers.
//
// File developed by: William F Godoy, godoywf@ornl.gov, Oak Ridge National Laboratory
//
// File created by: William F Godoy, godoywf@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "AGPDeterminant.h"
#include "Numerics/DeterminantOperators.h"

namespace qmcplusplus
{

template<class T>
T AGPDeterminant::do_ratioT(ParticleSet& P, int iat)
{
  UpdateMode = ORB_PBYP_RATIO;
  //GeminalBasis->evaluate(P,iat);
  GeminalBasis->evaluateForPtclMove(P, iat); //@@
  //To dune with gemv
  //BLAS::gemv(Lambda.rows(),Lambda.cols(), Lambda.data(), GeminalBasis->y(0), phiT[iat]);
  //const ValueType* restrict y_ptr=GeminalBasis->y(0);
  const BasisSetType::ValueType* restrict y_ptr = GeminalBasis->Phi.data(); //@@
  if (iat < Nup)
  {
    for (int d = 0, jat = Nup; d < Ndown; d++, jat++)
      psiU[d] = simd::dot(y_ptr, phiT[jat], BasisSize);
    //psiU[d]=BLAS::dot(BasisSize,y_ptr,phiT[jat]);
    //unpaired block Ndown x unpaired
    for (int d = Ndown, unpaired = 0; d < Nup; d++, unpaired++)
      //psiU[d] = BLAS::dot(BasisSize,LambdaUP[unpaired],y_ptr);
      psiU[d] = simd::dot(LambdaUP[unpaired], y_ptr, BasisSize);
    //curRatio=DetRatio(psiM, psiU.data(),iat);
    curRatio = DetRatioByRow(psiM, psiU, iat);
  }
  else
  {
    for (int u = 0; u < Nup; u++)
      //psiD[u]=BLAS::dot(BasisSize,y_ptr,phiT[u]);
      psiD[u] = simd::dot(y_ptr, phiT[u], BasisSize);
    //curRatio=DetRatioTranspose(psiM, psiD.data(),iat-Nup);
    curRatio = DetRatioByColumn(psiM, psiD, iat - Nup);
  }
  return curRatio;
}

} // namespace qmcplusplus