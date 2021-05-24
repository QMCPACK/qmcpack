//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Bryan Clark, bclark@Princeton.edu, Princeton University
//////////////////////////////////////////////////////////////////////////////////////


/**@file
 * @brief Implement build functions: Function bodies are too big to be in a header file
 */
#include "QMCWaveFunctions/Fermion/MultiDiracDeterminant.h"
#include "Numerics/MatrixOperators.h"

namespace qmcplusplus
{
/** shared function used by BuildDotProductsAndCalculateRatios */
void MultiDiracDeterminant::BuildDotProductsAndCalculateRatios_impl(int ref,
                                                                    ValueType det0,
                                                                    ValueType* restrict ratios,
                                                                    const ValueMatrix_t& psiinv,
                                                                    const ValueMatrix_t& psi,
                                                                    ValueMatrix_t& dotProducts,
                                                                    const std::vector<int>& data,
                                                                    const std::vector<std::pair<int, int>>& pairs,
                                                                    const std::vector<RealType>& sign)
{
  buildTableTimer.start();
  const size_t num    = psi.extent(1);
  const size_t npairs = pairs.size();
  //MatrixOperators::product_ABt(psiinv,psi,dotProducts);
  const std::pair<int, int>* restrict p = pairs.data();
  for (size_t i = 0; i < npairs; ++i)
  {
    const int I       = p[i].first;
    const int J       = p[i].second;
    dotProducts(I, J) = simd::dot(psiinv[I], psi[J], num);
  }
  buildTableTimer.stop();
  readMatTimer.start();
  std::vector<int>::const_iterator it2 = data.begin();
  const size_t nitems                  = sign.size();
  // explore Inclusive Scan for OpenMP
  for (size_t count = 0; count < nitems; ++count)
  {
    const size_t n = *it2;
    //ratios[count]=(count!=ref)?sign[count]*det0*CalculateRatioFromMatrixElements(n,dotProducts,it2+1):det0;
    if (count != ref)
      ratios[count] = sign[count] * det0 * CalculateRatioFromMatrixElements(n, dotProducts, it2 + 1);
    it2 += 3 * n + 1;
  }
  ratios[ref] = det0;
  readMatTimer.stop();
}

void MultiDiracDeterminant::BuildDotProductsAndCalculateRatios(int ref,
                                                               int iat,
                                                               ValueVector_t& ratios,
                                                               const ValueMatrix_t& psiinv,
                                                               const ValueMatrix_t& psi,
                                                               ValueMatrix_t& dotProducts,
                                                               const std::vector<int>& data,
                                                               const std::vector<std::pair<int, int>>& pairs,
                                                               const std::vector<RealType>& sign)
{
  BuildDotProductsAndCalculateRatios_impl(ref, ratios[ref], ratios.data(), psiinv, psi, dotProducts, data, pairs, sign);
#if 0
    buildTableTimer.start();
    ValueType det0 = ratios[ref];
    int num=psi.extent(1);
    std::vector<std::pair<int,int> >::iterator it(pairs.begin()), last(pairs.end());
    while(it != last)
    {
      dotProducts((*it).first,(*it).second) = simd::dot(psiinv[(*it).first],psi[(*it).second],num);
      it++;
    }
    std::vector<int>::iterator it2 = data.begin();
    int count= 0;  // number of determinants processed
    while(it2 != data.end())
    {
      const int n = *it2; // number of excitations
      if(count == ref)
      {
        it2+=3*n+1;  // number of integers used to encode the current excitation
        count++;
        continue;
      }
      ratios[count] = sign[count]*det0*CalculateRatioFromMatrixElements(n,dotProducts,it2+1);
      count++;
      it2+=3*n+1;
    }
    readMatTimer.stop();
#endif
}

void MultiDiracDeterminant::BuildDotProductsAndCalculateRatios(int ref,
                                                               int iat,
                                                               GradMatrix_t& ratios,
                                                               ValueMatrix_t& psiinv,
                                                               ValueMatrix_t& psi,
                                                               ValueMatrix_t& dotProducts,
                                                               std::vector<int>& data,
                                                               std::vector<std::pair<int, int>>& pairs,
                                                               std::vector<RealType>& sign,
                                                               int dx)
{
  const ValueType det0 = ratios(ref, iat)[dx];
  BuildDotProductsAndCalculateRatios_impl(ref, det0, WorkSpace.data(), psiinv, psi, dotProducts, data, pairs, sign);
  for (size_t count = 0; count < getNumDets(); ++count)
    ratios(count, iat)[dx] = WorkSpace[count];
#if 0
    ValueType det0 = ratios(ref,iat)[dx];
    buildTableGradTimer.start();
    int num=psi.extent(1);
    std::vector<std::pair<int,int> >::iterator it(pairs.begin()), last(pairs.end());
    while(it != last)
    {
      dotProducts((*it).first,(*it).second) = simd::dot(psiinv[(*it).first],psi[(*it).second],num);
      it++;
    }
    buildTableGradTimer.stop();
    readMatGradTimer.start();
    std::vector<int>::iterator it2 = data.begin();
    int count= 0;  // number of determinants processed
    while(it2 != data.end())
    {
      int n = *it2; // number of excitations
      if(count == ref)
      {
        it2+=3*n+1;  // number of integers used to encode the current excitation
        count++;
        continue;
      }
      ratios(count,iat)[dx] = sign[count]*det0*CalculateRatioFromMatrixElements(n,dotProducts,it2+1);
      count++;
      it2+=3*n+1;
    }
    readMatGradTimer.stop();
#endif
}

void MultiDiracDeterminant::BuildDotProductsAndCalculateRatios(int ref,
                                                               int iat,
                                                               ValueMatrix_t& ratios,
                                                               ValueMatrix_t& psiinv,
                                                               ValueMatrix_t& psi,
                                                               ValueMatrix_t& dotProducts,
                                                               std::vector<int>& data,
                                                               std::vector<std::pair<int, int>>& pairs,
                                                               std::vector<RealType>& sign)
{
  const ValueType det0 = ratios(ref, iat);
  BuildDotProductsAndCalculateRatios_impl(ref, det0, WorkSpace.data(), psiinv, psi, dotProducts, data, pairs, sign);
  //splatt
  for (size_t count = 0; count < getNumDets(); ++count)
    ratios(count, iat) = WorkSpace[count];
#if 0
    ValueType det0 = ratios(ref,iat);
    int num=psi.extent(1);
    std::vector<std::pair<int,int> >::iterator it(pairs.begin()), last(pairs.end());
    while(it != last)
    {
      dotProducts((*it).first,(*it).second) = simd::dot(psiinv[(*it).first],psi[(*it).second],num);
      it++;
    }
    std::vector<int>::iterator it2 = data.begin();
    int count= 0;  // number of determinants processed
    while(it2 != data.end())
    {
      int n = *it2; // number of excitations
      if(count == ref)
      {
        it2+=3*n+1;  // number of integers used to encode the current excitation
        count++;
        continue;
      }
      ratios(count,iat) = sign[count]*det0*CalculateRatioFromMatrixElements(n,dotProducts,it2+1);
      count++;
      it2+=3*n+1;
    }
#endif
}

void MultiDiracDeterminant::evaluateDetsForPtclMove(const ParticleSet& P, int iat, int refPtcl)
{
  UpdateMode = ORB_PBYP_RATIO;
  RatioTimer.start();
  evalOrbTimer.start();
  Phi->evaluateValue(P, iat, psiV);
  evalOrbTimer.stop();
  const int WorkingIndex = (refPtcl < 0 ? iat : refPtcl) - FirstIndex;
  assert(WorkingIndex >= 0 && WorkingIndex < LastIndex - FirstIndex);
  if (NumPtcls == 1)
  {
    std::vector<ci_configuration2>::iterator it(ciConfigList->begin());
    std::vector<ci_configuration2>::iterator last(ciConfigList->end());
    ValueVector_t::iterator det(new_detValues.begin());
    while (it != last)
    {
      size_t orb = (it++)->occup[0];
      *(det++)   = psiV[orb];
    }
  }
  else
  {
    const auto& confgList = *ciConfigList;
    //std::vector<int>::iterator it(confgList[ReferenceDeterminant].occup.begin());
    auto it(confgList[ReferenceDeterminant].occup.begin());
    // mmorales: the only reason this is here is because
    // NonlocalECP do not necessarily call rejectMove after
    // calling ratio(), and even if the move is rejected
    // this matrix needs to be restored
    // If we always restore after ratio, then this is not needed
    // For efficiency reasons, I don't do this for ratioGrad or ratio(P,dG,dL)
    ExtraStuffTimer.start();
    psiMinv_temp = psiMinv;
    for (size_t i = 0; i < NumPtcls; i++)
      psiV_temp[i] = psiV[*(it++)];
    ValueType ratioRef                  = DetRatioByColumn(psiMinv_temp, psiV_temp, WorkingIndex);
    new_detValues[ReferenceDeterminant] = ratioRef * detValues[ReferenceDeterminant];
    InverseUpdateByColumn(psiMinv_temp, psiV_temp, workV1, workV2, WorkingIndex, ratioRef);
    for (size_t i = 0; i < NumOrbitals; i++)
      TpsiM(i, WorkingIndex) = psiV[i];
    ExtraStuffTimer.stop();
    BuildDotProductsAndCalculateRatios(ReferenceDeterminant, WorkingIndex, new_detValues, psiMinv_temp, TpsiM,
                                       dotProducts, *detData, *uniquePairs, *DetSigns);
    // check comment above
    for (size_t i = 0; i < NumOrbitals; i++)
      TpsiM(i, WorkingIndex) = psiM(WorkingIndex, i);
  }
  RatioTimer.stop();
}

void MultiDiracDeterminant::evaluateDetsAndGradsForPtclMove(const ParticleSet& P, int iat)
{
  UpdateMode = ORB_PBYP_PARTIAL;
  evalOrb1Timer.start();
  Phi->evaluateVGL(P, iat, psiV, dpsiV, d2psiV);
  evalOrb1Timer.stop();
  const int WorkingIndex = iat - FirstIndex;
  assert(WorkingIndex >= 0 && WorkingIndex < LastIndex - FirstIndex);
  if (NumPtcls == 1)
  {
    std::vector<ci_configuration2>::iterator it(ciConfigList->begin());
    std::vector<ci_configuration2>::iterator last(ciConfigList->end());
    ValueVector_t::iterator det(new_detValues.begin());
    GradMatrix_t::iterator grad(new_grads.begin());
    while (it != last)
    {
      size_t orb = (it++)->occup[0];
      *(det++)   = psiV[orb];
      *(grad++)  = dpsiV[orb];
    }
  }
  else
  {
    ExtraStuffTimer.start();
    //mmorales: check comment above
    psiMinv_temp          = psiMinv;
    const auto& confgList = *ciConfigList;
    //std::vector<int>::iterator it(confgList[ReferenceDeterminant].occup.begin());
    auto it(confgList[ReferenceDeterminant].occup.begin());
    GradType ratioGradRef;
    for (size_t i = 0; i < NumPtcls; i++)
    {
      psiV_temp[i] = psiV[*it];
      ratioGradRef += psiMinv_temp(i, WorkingIndex) * dpsiV[*it];
      it++;
    }
    ValueType ratioRef                            = DetRatioByColumn(psiMinv_temp, psiV_temp, WorkingIndex);
    new_grads(ReferenceDeterminant, WorkingIndex) = ratioGradRef * detValues[ReferenceDeterminant];
    new_detValues[ReferenceDeterminant]           = ratioRef * detValues[ReferenceDeterminant];
    InverseUpdateByColumn(psiMinv_temp, psiV_temp, workV1, workV2, WorkingIndex, ratioRef);
    for (size_t i = 0; i < NumOrbitals; i++)
      TpsiM(i, WorkingIndex) = psiV[i];
    ExtraStuffTimer.stop();
    BuildDotProductsAndCalculateRatios(ReferenceDeterminant, WorkingIndex, new_detValues, psiMinv_temp, TpsiM,
                                       dotProducts, *detData, *uniquePairs, *DetSigns);
    for (size_t idim = 0; idim < OHMMS_DIM; idim++)
    {
      ExtraStuffTimer.start();
      //dpsiMinv = psiMinv_temp;
      dpsiMinv = psiMinv;
      it       = confgList[ReferenceDeterminant].occup.begin();
      for (size_t i = 0; i < NumPtcls; i++)
        psiV_temp[i] = dpsiV[*(it++)][idim];
      InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, WorkingIndex, ratioGradRef[idim]);
      for (size_t i = 0; i < NumOrbitals; i++)
        TpsiM(i, WorkingIndex) = dpsiV[i][idim];
      ExtraStuffTimer.stop();
      BuildDotProductsAndCalculateRatios(ReferenceDeterminant, WorkingIndex, new_grads, dpsiMinv, TpsiM, dotProducts,
                                         *detData, *uniquePairs, *DetSigns, idim);
    }
    // check comment above
    for (int i = 0; i < NumOrbitals; i++)
      TpsiM(i, WorkingIndex) = psiM(WorkingIndex, i);
  }
}

void MultiDiracDeterminant::evaluateDetsAndGradsForPtclMoveWithSpin(const ParticleSet& P, int iat)
{
  UpdateMode = ORB_PBYP_PARTIAL;
  evalOrb1Timer.start();
  Phi->evaluateVGL_spin(P, iat, psiV, dpsiV, d2psiV, dspin_psiV);
  evalOrb1Timer.stop();
  const int WorkingIndex = iat - FirstIndex;
  assert(WorkingIndex >= 0 && WorkingIndex < LastIndex - FirstIndex);
  if (NumPtcls == 1)
  {
    std::vector<ci_configuration2>::iterator it(ciConfigList->begin());
    std::vector<ci_configuration2>::iterator last(ciConfigList->end());
    ValueVector_t::iterator det(new_detValues.begin());
    GradMatrix_t::iterator grad(new_grads.begin());
    ValueMatrix_t::iterator spingrad(new_spingrads.begin());
    while (it != last)
    {
      size_t orb    = (it++)->occup[0];
      *(det++)      = psiV[orb];
      *(grad++)     = dpsiV[orb];
      *(spingrad++) = dspin_psiV[orb];
    }
  }
  else
  {
    ExtraStuffTimer.start();
    //mmorales: check comment above
    psiMinv_temp          = psiMinv;
    const auto& confgList = *ciConfigList;
    //std::vector<int>::iterator it(confgList[ReferenceDeterminant].occup.begin());
    auto it(confgList[ReferenceDeterminant].occup.begin());
    GradType ratioGradRef;
    ValueType ratioSpinGradRef;
    for (size_t i = 0; i < NumPtcls; i++)
    {
      psiV_temp[i] = psiV[*it];
      ratioGradRef += psiMinv_temp(i, WorkingIndex) * dpsiV[*it];
      ratioSpinGradRef += psiMinv_temp(i, WorkingIndex) * dspin_psiV[*it];
      it++;
    }
    ValueType ratioRef                                = DetRatioByColumn(psiMinv_temp, psiV_temp, WorkingIndex);
    new_grads(ReferenceDeterminant, WorkingIndex)     = ratioGradRef * detValues[ReferenceDeterminant];
    new_spingrads(ReferenceDeterminant, WorkingIndex) = ratioSpinGradRef * detValues[ReferenceDeterminant];
    new_detValues[ReferenceDeterminant]               = ratioRef * detValues[ReferenceDeterminant];
    InverseUpdateByColumn(psiMinv_temp, psiV_temp, workV1, workV2, WorkingIndex, ratioRef);
    for (size_t i = 0; i < NumOrbitals; i++)
      TpsiM(i, WorkingIndex) = psiV[i];
    ExtraStuffTimer.stop();
    BuildDotProductsAndCalculateRatios(ReferenceDeterminant, WorkingIndex, new_detValues, psiMinv_temp, TpsiM,
                                       dotProducts, *detData, *uniquePairs, *DetSigns);
    for (size_t idim = 0; idim < OHMMS_DIM; idim++)
    {
      ExtraStuffTimer.start();
      //dpsiMinv = psiMinv_temp;
      dpsiMinv = psiMinv;
      it       = confgList[ReferenceDeterminant].occup.begin();
      for (size_t i = 0; i < NumPtcls; i++)
        psiV_temp[i] = dpsiV[*(it++)][idim];
      InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, WorkingIndex, ratioGradRef[idim]);
      for (size_t i = 0; i < NumOrbitals; i++)
        TpsiM(i, WorkingIndex) = dpsiV[i][idim];
      ExtraStuffTimer.stop();
      BuildDotProductsAndCalculateRatios(ReferenceDeterminant, WorkingIndex, new_grads, dpsiMinv, TpsiM, dotProducts,
                                         *detData, *uniquePairs, *DetSigns, idim);
    }
    ExtraStuffTimer.start();
    dpsiMinv = psiMinv;
    it       = confgList[ReferenceDeterminant].occup.begin();
    for (size_t i = 0; i < NumPtcls; i++)
      psiV_temp[i] = dspin_psiV[*(it++)];
    InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, WorkingIndex, ratioSpinGradRef);
    for (size_t i = 0; i < NumOrbitals; i++)
      TpsiM(i, WorkingIndex) = dspin_psiV[i];
    ExtraStuffTimer.stop();
    BuildDotProductsAndCalculateRatios(ReferenceDeterminant, WorkingIndex, new_spingrads, dpsiMinv, TpsiM, dotProducts,
                                       *detData, *uniquePairs, *DetSigns);
    // check comment above
    for (int i = 0; i < NumOrbitals; i++)
      TpsiM(i, WorkingIndex) = psiM(WorkingIndex, i);
  }
}

void MultiDiracDeterminant::evaluateGrads(ParticleSet& P, int iat)
{
  const int WorkingIndex = iat - FirstIndex;
  assert(WorkingIndex >= 0 && WorkingIndex < LastIndex - FirstIndex);

  const auto& confgList = *ciConfigList;
  auto it               = confgList[0].occup.begin(); //just to avoid using the type
  //std::vector<size_t>::iterator it;
  if (NumPtcls == 1)
  {
    std::vector<ci_configuration2>::const_iterator it(confgList.begin());
    std::vector<ci_configuration2>::const_iterator last(confgList.end());
    GradMatrix_t::iterator grad(grads.begin());
    while (it != last)
    {
      size_t orb = (it++)->occup[0];
      *(grad++)  = dpsiM(0, orb);
    }
  }
  else
  {
    for (size_t idim = 0; idim < OHMMS_DIM; idim++)
    {
      //dpsiMinv = psiMinv_temp;
      dpsiMinv         = psiMinv;
      it               = confgList[ReferenceDeterminant].occup.begin();
      ValueType ratioG = 0.0;
      for (size_t i = 0; i < NumPtcls; i++)
      {
        psiV_temp[i] = dpsiM(WorkingIndex, *it)[idim];
        ratioG += psiMinv(i, WorkingIndex) * dpsiM(WorkingIndex, *it)[idim];
        it++;
      }
      grads(ReferenceDeterminant, WorkingIndex)[idim] = ratioG * detValues[ReferenceDeterminant];
      InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, WorkingIndex, ratioG);
      for (size_t i = 0; i < NumOrbitals; i++)
        TpsiM(i, WorkingIndex) = dpsiM(WorkingIndex, i)[idim];
      BuildDotProductsAndCalculateRatios(ReferenceDeterminant, WorkingIndex, grads, dpsiMinv, TpsiM, dotProducts,
                                         *detData, *uniquePairs, *DetSigns, idim);
    }
    // check comment above
    for (size_t i = 0; i < NumOrbitals; i++)
      TpsiM(i, WorkingIndex) = psiM(WorkingIndex, i);
  }
}

void MultiDiracDeterminant::evaluateGradsWithSpin(ParticleSet& P, int iat)
{
  const int WorkingIndex = iat - FirstIndex;
  assert(WorkingIndex >= 0 && WorkingIndex < LastIndex - FirstIndex);

  const auto& confgList = *ciConfigList;
  auto it               = confgList[0].occup.begin(); //just to avoid using the type
  //std::vector<size_t>::iterator it;
  if (NumPtcls == 1)
  {
    std::vector<ci_configuration2>::const_iterator it(confgList.begin());
    std::vector<ci_configuration2>::const_iterator last(confgList.end());
    GradMatrix_t::iterator grad(grads.begin());
    ValueMatrix_t::iterator spingrad(spingrads.begin());
    while (it != last)
    {
      size_t orb    = (it++)->occup[0];
      *(grad++)     = dpsiM(0, orb);
      *(spingrad++) = dspin_psiM(0, orb);
    }
  }
  else
  {
    for (size_t idim = 0; idim < OHMMS_DIM; idim++)
    {
      //dpsiMinv = psiMinv_temp;
      dpsiMinv         = psiMinv;
      it               = confgList[ReferenceDeterminant].occup.begin();
      ValueType ratioG = 0.0;
      for (size_t i = 0; i < NumPtcls; i++)
      {
        psiV_temp[i] = dpsiM(WorkingIndex, *it)[idim];
        ratioG += psiMinv(i, WorkingIndex) * dpsiM(WorkingIndex, *it)[idim];
        it++;
      }
      grads(ReferenceDeterminant, WorkingIndex)[idim] = ratioG * detValues[ReferenceDeterminant];
      InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, WorkingIndex, ratioG);
      for (size_t i = 0; i < NumOrbitals; i++)
        TpsiM(i, WorkingIndex) = dpsiM(WorkingIndex, i)[idim];
      BuildDotProductsAndCalculateRatios(ReferenceDeterminant, WorkingIndex, grads, dpsiMinv, TpsiM, dotProducts,
                                         *detData, *uniquePairs, *DetSigns, idim);
    }
    dpsiMinv          = psiMinv;
    it                = confgList[ReferenceDeterminant].occup.begin();
    ValueType ratioSG = 0.0;
    for (size_t i = 0; i < NumPtcls; i++)
    {
      psiV_temp[i] = dspin_psiM(WorkingIndex, *it);
      ratioSG += psiMinv(i, WorkingIndex) * dspin_psiM(WorkingIndex, *it);
      it++;
    }
    spingrads(ReferenceDeterminant, WorkingIndex) = ratioSG * detValues[ReferenceDeterminant];
    InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, WorkingIndex, ratioSG);
    for (size_t i = 0; i < NumOrbitals; i++)
      TpsiM(i, WorkingIndex) = dspin_psiM(WorkingIndex, i);
    BuildDotProductsAndCalculateRatios(ReferenceDeterminant, WorkingIndex, spingrads, dpsiMinv, TpsiM, dotProducts,
                                       *detData, *uniquePairs, *DetSigns);
    // check comment above
    for (size_t i = 0; i < NumOrbitals; i++)
      TpsiM(i, WorkingIndex) = psiM(WorkingIndex, i);
  }
}

void MultiDiracDeterminant::evaluateAllForPtclMove(const ParticleSet& P, int iat)
{
  UpdateMode = ORB_PBYP_ALL;
  Phi->evaluateVGL(P, iat, psiV, dpsiV, d2psiV);
  const int WorkingIndex = iat - FirstIndex;
  assert(WorkingIndex >= 0 && WorkingIndex < LastIndex - FirstIndex);
  if (NumPtcls == 1)
  {
    std::vector<ci_configuration2>::const_iterator it(ciConfigList->begin());
    std::vector<ci_configuration2>::const_iterator last(ciConfigList->end());
    ValueVector_t::iterator det(new_detValues.begin());
    ValueMatrix_t::iterator lap(new_lapls.begin());
    GradMatrix_t::iterator grad(new_grads.begin());
    while (it != last)
    {
      size_t orb = (it++)->occup[0];
      *(det++)   = psiV[orb];
      *(lap++)   = d2psiV[orb];
      *(grad++)  = dpsiV[orb];
    }
  }
  else
  {
    //mmorales: check comment above
    psiMinv_temp = psiMinv;

    const auto& confgList = *ciConfigList;

    //std::vector<int>::iterator it(confgList[ReferenceDeterminant].occup.begin());
    auto it(confgList[ReferenceDeterminant].occup.begin());
    //GradType ratioGradRef;
    //ValueType ratioLapl = 0.0;
    for (size_t i = 0; i < NumPtcls; i++)
    {
      psiV_temp[i] = psiV[*it];
      //ratioGradRef += psiMinv_temp(i,WorkingIndex)*dpsiV(*it);
      //ratioLapl += psiMinv_temp(i,WorkingIndex)*d2psiV(*it);
      it++;
    }
    ValueType det0, ratioRef = DetRatioByColumn(psiMinv_temp, psiV_temp, WorkingIndex);
    //new_lapls(ReferenceDeterminant,WorkingIndex) = ratioLapl*detValues[ReferenceDeterminant];
    //new_grads(ReferenceDeterminant,WorkingIndex) = ratioGradRef*detValues[ReferenceDeterminant];
    det0                                = ratioRef * detValues[ReferenceDeterminant];
    new_detValues[ReferenceDeterminant] = det0;
    InverseUpdateByColumn(psiMinv_temp, psiV_temp, workV1, workV2, WorkingIndex, ratioRef);
    for (size_t i = 0; i < NumOrbitals; i++)
      TpsiM(i, WorkingIndex) = psiV[i];
    BuildDotProductsAndCalculateRatios(ReferenceDeterminant, WorkingIndex, new_detValues, psiMinv_temp, TpsiM,
                                       dotProducts, *detData, *uniquePairs, *DetSigns);
    for (size_t jat = 0; jat < NumPtcls; jat++)
    {
      it = confgList[ReferenceDeterminant].occup.begin();
      GradType gradRatio; // = 0.0;
      ValueType ratioLapl = 0.0;
      if (jat == WorkingIndex)
      {
        for (size_t i = 0; i < NumPtcls; i++)
        {
          gradRatio += psiMinv_temp(i, jat) * dpsiV[*it];
          ratioLapl += psiMinv_temp(i, jat) * d2psiV[*it];
          it++;
        }
        new_grads(ReferenceDeterminant, jat) = det0 * gradRatio;
        new_lapls(ReferenceDeterminant, jat) = det0 * ratioLapl;
        for (size_t idim = 0; idim < OHMMS_DIM; idim++)
        {
          dpsiMinv = psiMinv_temp;
          it       = confgList[ReferenceDeterminant].occup.begin();
          for (size_t i = 0; i < NumPtcls; i++)
            psiV_temp[i] = dpsiV[*(it++)][idim];
          InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, jat, gradRatio[idim]);
          for (size_t i = 0; i < NumOrbitals; i++)
            TpsiM(i, jat) = dpsiV[i][idim];
          BuildDotProductsAndCalculateRatios(ReferenceDeterminant, jat, new_grads, dpsiMinv, TpsiM, dotProducts,
                                             *detData, *uniquePairs, *DetSigns, idim);
        }
        dpsiMinv = psiMinv_temp;
        it       = confgList[ReferenceDeterminant].occup.begin();
        for (size_t i = 0; i < NumPtcls; i++)
          psiV_temp[i] = d2psiV[*(it++)];
        InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, jat, ratioLapl);
        for (size_t i = 0; i < NumOrbitals; i++)
          TpsiM(i, jat) = d2psiV[i];
        BuildDotProductsAndCalculateRatios(ReferenceDeterminant, jat, new_lapls, dpsiMinv, TpsiM, dotProducts, *detData,
                                           *uniquePairs, *DetSigns);
        for (size_t i = 0; i < NumOrbitals; i++)
          TpsiM(i, jat) = psiV[i];
      }
      else
      {
        for (size_t i = 0; i < NumPtcls; i++)
        {
          gradRatio += psiMinv_temp(i, jat) * dpsiM(jat, *it);
          ratioLapl += psiMinv_temp(i, jat) * d2psiM(jat, *it);
          it++;
        }
        new_grads(ReferenceDeterminant, jat) = det0 * gradRatio;
        new_lapls(ReferenceDeterminant, jat) = det0 * ratioLapl;
        for (size_t idim = 0; idim < OHMMS_DIM; idim++)
        {
          dpsiMinv = psiMinv_temp;
          it       = confgList[ReferenceDeterminant].occup.begin();
          for (size_t i = 0; i < NumPtcls; i++)
            psiV_temp[i] = dpsiM(jat, *(it++))[idim];
          InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, jat, gradRatio[idim]);
          for (size_t i = 0; i < NumOrbitals; i++)
            TpsiM(i, jat) = dpsiM(jat, i)[idim];
          BuildDotProductsAndCalculateRatios(ReferenceDeterminant, jat, new_grads, dpsiMinv, TpsiM, dotProducts,
                                             *detData, *uniquePairs, *DetSigns, idim);
        }
        dpsiMinv = psiMinv_temp;
        it       = confgList[ReferenceDeterminant].occup.begin();
        for (size_t i = 0; i < NumPtcls; i++)
          psiV_temp[i] = d2psiM(jat, *(it++));
        InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, jat, ratioLapl);
        for (size_t i = 0; i < NumOrbitals; i++)
          TpsiM(i, jat) = d2psiM(jat, i);
        BuildDotProductsAndCalculateRatios(ReferenceDeterminant, jat, new_lapls, dpsiMinv, TpsiM, dotProducts, *detData,
                                           *uniquePairs, *DetSigns);
        for (size_t i = 0; i < NumOrbitals; i++)
          TpsiM(i, jat) = psiM(jat, i);
      }
    } // jat
    // check comment above
    for (size_t i = 0; i < NumOrbitals; i++)
      TpsiM(i, WorkingIndex) = psiM(WorkingIndex, i);
  }
}


} // namespace qmcplusplus
