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
                                                                    const ValueMatrix& psiinv,
                                                                    const ValueMatrix& psi,
                                                                    ValueMatrix& dotProducts,
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
                                                               ValueVector& ratios,
                                                               const ValueMatrix& psiinv,
                                                               const ValueMatrix& psi,
                                                               ValueMatrix& dotProducts,
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
                                                               GradMatrix& ratios,
                                                               ValueMatrix& psiinv,
                                                               ValueMatrix& psi,
                                                               ValueMatrix& dotProducts,
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
                                                               ValueMatrix& ratios,
                                                               ValueMatrix& psiinv,
                                                               ValueMatrix& psi,
                                                               ValueMatrix& dotProducts,
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

void MultiDiracDeterminant::mw_evaluateDetsForPtclMove(const RefVectorWithLeader<MultiDiracDeterminant>& det_list,
                                                       const RefVectorWithLeader<ParticleSet>& P_list,
                                                       int iat)
{
  const int nw                      = det_list.size();
  MultiDiracDeterminant& det_leader = det_list.getLeader();
  det_leader.RatioTimer.start();

  det_leader.UpdateMode = ORB_PBYP_RATIO;
  /*  FOR YE: THIS IS NOT compiling...
  std::vector<RefVector<Vector<ValueVector>>> psiV_list;
  for (size_t iw=0;iw<nw;iw++)
  {
    MultiDiracDeterminant& det= (det_list[iw]);
    psiV_list.push_back(det.psiV);
  }
  det_leader.evalOrbTimer.start()
  det_leader.Phi->mw_evaluateValue(P_list[iw], iat, psiV_list);
*/


  for (size_t iw = 0; iw < nw; iw++)
  {
    MultiDiracDeterminant& det = (det_list[iw]);
    det.UpdateMode             = ORB_PBYP_RATIO;
    det.evalOrbTimer.start();
    det.Phi->evaluateValue(P_list[iw], iat, det.psiV);
    det.evalOrbTimer.stop();
    const int WorkingIndex = iat - det.FirstIndex;
    const auto& confgList  = *det.ciConfigList;
    ///std::vector<int>::iterator it(confgList[ReferenceDeterminant].occup.begin());
    auto it(confgList[det.ReferenceDeterminant].occup.begin());
    // mmorales: the only reason this is here is because
    // NonlocalECP do not necessarily call rejectMove after
    // calling ratio(), and even if the move is rejected
    // this matrix needs to be restored
    // If we always restore after ratio, then this is not needed
    // For efficiency reasons, I don't do this for ratioGrad or ratio(P,dG,dL)
    det.ExtraStuffTimer.start();
    det.psiMinv_temp = det.psiMinv;
    for (size_t i = 0; i < det_leader.NumPtcls; i++)
      det.psiV_temp[i] = det.psiV[*(it++)];

    //template<typename MatA, typename VecB>
    //inline typename MatA::value_type DetRatioByColumn(const MatA& Minv, const VecB& newv, int colchanged)
    //{
    //  //use BLAS dot since the stride is not uniform
    //  //  return simd::dot(Minv.cols(), Minv.data() + colchanged, Minv.cols(), newv.data(), 1);
    //  //  }
    //  //
    det.curRatio                                     = DetRatioByColumn(det.psiMinv_temp, det.psiV_temp, WorkingIndex);
    det.new_ratios_to_ref_[det.ReferenceDeterminant] = ValueType(1);
    InverseUpdateByColumn(det.psiMinv_temp, det.psiV_temp, det.workV1, det.workV2, WorkingIndex, det.curRatio);
    for (size_t i = 0; i < det_leader.NumOrbitals; i++)
      det.TpsiM(i, WorkingIndex) = det.psiV[i];
    det.ExtraStuffTimer.stop();
    det.BuildDotProductsAndCalculateRatios(det.ReferenceDeterminant, det.new_ratios_to_ref_, det.psiMinv_temp,
                                           det.TpsiM, det.dotProducts, *det.detData, *det.uniquePairs, *det.DetSigns);

    for (size_t i = 0; i < det_leader.NumOrbitals; i++)
      det.TpsiM(i, WorkingIndex) = det.psiM(WorkingIndex, i);
  }
  det_leader.RatioTimer.stop();
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
  auto ratio_old_ref_det                   = DetRatioByColumn(psiMinv_temp, psiV_temp, WorkingIndex);
  curRatio                                 = ratio_old_ref_det;
  new_ratios_to_ref_[ReferenceDeterminant] = ValueType(1);
  InverseUpdateByColumn(psiMinv_temp, psiV_temp, workV1, workV2, WorkingIndex, ratio_old_ref_det);
  for (size_t i = 0; i < NumOrbitals; i++)
    TpsiM(i, WorkingIndex) = psiV[i];
  ExtraStuffTimer.stop();
  BuildDotProductsAndCalculateRatios(ReferenceDeterminant, new_ratios_to_ref_, psiMinv_temp, TpsiM, dotProducts,
                                     *detData, *uniquePairs, *DetSigns);
  // check comment above
  for (size_t i = 0; i < NumOrbitals; i++)
    TpsiM(i, WorkingIndex) = psiM(WorkingIndex, i);
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
  curRatio                                      = DetRatioByColumn(psiMinv_temp, psiV_temp, WorkingIndex);
  new_grads(ReferenceDeterminant, WorkingIndex) = ratioGradRef / curRatio;
  new_ratios_to_ref_[ReferenceDeterminant]      = 1.0;
  InverseUpdateByColumn(psiMinv_temp, psiV_temp, workV1, workV2, WorkingIndex, curRatio);
  for (size_t i = 0; i < NumOrbitals; i++)
    TpsiM(i, WorkingIndex) = psiV[i];
  ExtraStuffTimer.stop();
  BuildDotProductsAndCalculateRatios(ReferenceDeterminant, new_ratios_to_ref_, psiMinv_temp, TpsiM, dotProducts,
                                     *detData, *uniquePairs, *DetSigns);
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

void MultiDiracDeterminant::evaluateDetsAndGradsForPtclMoveWithSpin(const ParticleSet& P, int iat)
{
  assert(P.isSpinor() == is_spinor_);
  UpdateMode = ORB_PBYP_PARTIAL;
  evalOrb1Timer.start();
  Phi->evaluateVGL_spin(P, iat, psiV, dpsiV, d2psiV, dspin_psiV);
  evalOrb1Timer.stop();
  const int WorkingIndex = iat - FirstIndex;
  assert(WorkingIndex >= 0 && WorkingIndex < LastIndex - FirstIndex);
  ExtraStuffTimer.start();
  //mmorales: check comment above
  psiMinv_temp          = psiMinv;
  const auto& confgList = *ciConfigList;
  //std::vector<int>::iterator it(confgList[ReferenceDeterminant].occup.begin());
  auto it(confgList[ReferenceDeterminant].occup.begin());
  GradType ratioGradRef;
  ValueType ratioSpinGradRef = 0.0;
  for (size_t i = 0; i < NumPtcls; i++)
  {
    psiV_temp[i] = psiV[*it];
    ratioGradRef += psiMinv_temp(i, WorkingIndex) * dpsiV[*it];
    ratioSpinGradRef += psiMinv_temp(i, WorkingIndex) * dspin_psiV[*it];
    it++;
  }
  curRatio                                          = DetRatioByColumn(psiMinv_temp, psiV_temp, WorkingIndex);
  new_grads(ReferenceDeterminant, WorkingIndex)     = ratioGradRef / curRatio;
  new_spingrads(ReferenceDeterminant, WorkingIndex) = ratioSpinGradRef / curRatio;
  new_ratios_to_ref_[ReferenceDeterminant]          = 1.0;
  InverseUpdateByColumn(psiMinv_temp, psiV_temp, workV1, workV2, WorkingIndex, curRatio);
  for (size_t i = 0; i < NumOrbitals; i++)
    TpsiM(i, WorkingIndex) = psiV[i];
  ExtraStuffTimer.stop();
  BuildDotProductsAndCalculateRatios(ReferenceDeterminant, new_ratios_to_ref_, psiMinv_temp, TpsiM, dotProducts,
                                     *detData, *uniquePairs, *DetSigns);
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
  //Now compute the spin gradient, same procedure as normal gradient components above
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

void MultiDiracDeterminant::mw_evaluateDetsAndGradsForPtclMove(
    const RefVectorWithLeader<MultiDiracDeterminant>& det_list,
    const RefVectorWithLeader<ParticleSet>& P_list,
    int iat)
{
  const int nw                      = det_list.size();
  MultiDiracDeterminant& det_leader = det_list.getLeader();

  det_leader.UpdateMode = ORB_PBYP_PARTIAL;
  //det_leader.evalOrb1Timer.start();
  //mw_evaluateVGL(P_list[iw], iat, psiV_list, dpsiV_list, d2psiV_list);
  //det_leader.evalOrb1Timer.stop();

  for (size_t iw = 0; iw < nw; iw++)
  {
    MultiDiracDeterminant& det = (det_list[iw]);
    det.UpdateMode             = ORB_PBYP_PARTIAL;
    det.evalOrb1Timer.start();
    det.Phi->evaluateVGL(P_list[iw], iat, det.psiV, det.dpsiV, det.d2psiV);
    det.evalOrb1Timer.stop();
    const int WorkingIndex = iat - det.FirstIndex;
    det.ExtraStuffTimer.start();
    det.psiMinv_temp      = det.psiMinv;
    const auto& confgList = *det.ciConfigList;
    auto it(confgList[det.ReferenceDeterminant].occup.begin());
    GradType ratioGradRef;
    for (size_t i = 0; i < det_leader.NumPtcls; i++)
    {
      det.psiV_temp[i] = det.psiV[*it];
      ratioGradRef += det.psiMinv_temp(i, WorkingIndex) * det.dpsiV[*it];
      it++;
    }
    det.curRatio = DetRatioByColumn(det.psiMinv_temp, det.psiV_temp, WorkingIndex);
    det.new_grads(det.ReferenceDeterminant, WorkingIndex) = ratioGradRef / det.curRatio;
    det.new_ratios_to_ref_[det.ReferenceDeterminant]      = ValueType(1);
    InverseUpdateByColumn(det.psiMinv_temp, det.psiV_temp, det.workV1, det.workV2, WorkingIndex, det.curRatio);
    for (size_t i = 0; i < det_leader.NumOrbitals; i++)
      det.TpsiM(i, WorkingIndex) = det.psiV[i];
    det.ExtraStuffTimer.stop();
    det.BuildDotProductsAndCalculateRatios(det.ReferenceDeterminant, det.new_ratios_to_ref_, det.psiMinv_temp,
                                           det.TpsiM, det.dotProducts, *det.detData, *det.uniquePairs, *det.DetSigns);
    for (size_t idim = 0; idim < OHMMS_DIM; idim++)
    {
      det.ExtraStuffTimer.start();
      det.dpsiMinv = det.psiMinv;
      it           = confgList[det.ReferenceDeterminant].occup.begin();
      for (size_t i = 0; i < det_leader.NumPtcls; i++)
        det.psiV_temp[i] = det.dpsiV[*(it++)][idim];
      InverseUpdateByColumn(det.dpsiMinv, det.psiV_temp, det.workV1, det.workV2, WorkingIndex, ratioGradRef[idim]);
      for (size_t i = 0; i < det_leader.NumOrbitals; i++)
        det.TpsiM(i, WorkingIndex) = det.dpsiV[i][idim];
      det.ExtraStuffTimer.stop();
      det.BuildDotProductsAndCalculateRatios(det.ReferenceDeterminant, WorkingIndex, det.new_grads, det.dpsiMinv,
                                             det.TpsiM, det.dotProducts, *det.detData, *det.uniquePairs, *det.DetSigns,
                                             idim);
    }
    for (size_t i = 0; i < det_leader.NumOrbitals; i++)
      det.TpsiM(i, WorkingIndex) = det.psiM(WorkingIndex, i);
  }
}

void MultiDiracDeterminant::evaluateGrads(ParticleSet& P, int iat)
{
  const int WorkingIndex = iat - FirstIndex;
  assert(WorkingIndex >= 0 && WorkingIndex < LastIndex - FirstIndex);

  const auto& confgList = *ciConfigList;
  for (size_t idim = 0; idim < OHMMS_DIM; idim++)
  {
    //dpsiMinv = psiMinv_temp;
    dpsiMinv         = psiMinv;
    auto it          = confgList[ReferenceDeterminant].occup.begin();
    ValueType ratioG = 0.0;
    for (size_t i = 0; i < NumPtcls; i++)
    {
      psiV_temp[i] = dpsiM(WorkingIndex, *it)[idim];
      ratioG += psiMinv(i, WorkingIndex) * dpsiM(WorkingIndex, *it)[idim];
      it++;
    }
    grads(ReferenceDeterminant, WorkingIndex)[idim] = ratioG;
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

void MultiDiracDeterminant::evaluateGradsWithSpin(ParticleSet& P, int iat)
{
  assert(P.isSpinor() == is_spinor_);
  const int WorkingIndex = iat - FirstIndex;
  assert(WorkingIndex >= 0 && WorkingIndex < LastIndex - FirstIndex);

  const auto& confgList = *ciConfigList;
  for (size_t idim = 0; idim < OHMMS_DIM; idim++)
  {
    //dpsiMinv = psiMinv_temp;
    dpsiMinv         = psiMinv;
    auto it          = confgList[ReferenceDeterminant].occup.begin();
    ValueType ratioG = 0.0;
    for (size_t i = 0; i < NumPtcls; i++)
    {
      psiV_temp[i] = dpsiM(WorkingIndex, *it)[idim];
      ratioG += psiMinv(i, WorkingIndex) * dpsiM(WorkingIndex, *it)[idim];
      it++;
    }
    grads(ReferenceDeterminant, WorkingIndex)[idim] = ratioG;
    InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, WorkingIndex, ratioG);
    for (size_t i = 0; i < NumOrbitals; i++)
      TpsiM(i, WorkingIndex) = dpsiM(WorkingIndex, i)[idim];
    BuildDotProductsAndCalculateRatios(ReferenceDeterminant, WorkingIndex, grads, dpsiMinv, TpsiM, dotProducts,
                                       *detData, *uniquePairs, *DetSigns, idim);
  }

  //Now compute the spin gradient, same procedure as normal gradient components above
  dpsiMinv          = psiMinv;
  auto it           = confgList[ReferenceDeterminant].occup.begin();
  ValueType ratioSG = 0.0;
  for (size_t i = 0; i < NumPtcls; i++)
  {
    psiV_temp[i] = dspin_psiM(WorkingIndex, *it);
    ratioSG += psiMinv(i, WorkingIndex) * dspin_psiM(WorkingIndex, *it);
    it++;
  }
  spingrads(ReferenceDeterminant, WorkingIndex) = ratioSG;
  InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, WorkingIndex, ratioSG);
  for (size_t i = 0; i < NumOrbitals; i++)
    TpsiM(i, WorkingIndex) = dspin_psiM(WorkingIndex, i);
  BuildDotProductsAndCalculateRatios(ReferenceDeterminant, WorkingIndex, spingrads, dpsiMinv, TpsiM, dotProducts,
                                     *detData, *uniquePairs, *DetSigns);

  // check comment above
  for (size_t i = 0; i < NumOrbitals; i++)
    TpsiM(i, WorkingIndex) = psiM(WorkingIndex, i);
}

void MultiDiracDeterminant::mw_evaluateGrads(const RefVectorWithLeader<MultiDiracDeterminant>& det_list,
                                             const RefVectorWithLeader<ParticleSet>& P_list,
                                             int iat)
{
  const int nw                      = det_list.size();
  MultiDiracDeterminant& det_leader = det_list.getLeader();

  for (size_t iw = 0; iw < nw; iw++)
  {
    MultiDiracDeterminant& det = (det_list[iw]);
    const int WorkingIndex     = iat - det.FirstIndex;
    const auto& confgList      = *det.ciConfigList;

    for (size_t idim = 0; idim < OHMMS_DIM; idim++)
    {
      //dpsiMinv = psiMinv_temp;
      det.dpsiMinv     = det.psiMinv;
      auto it          = confgList[det.ReferenceDeterminant].occup.begin();
      ValueType ratioG = 0.0;
      for (size_t i = 0; i < det_leader.NumPtcls; i++)
      {
        det.psiV_temp[i] = det.dpsiM(WorkingIndex, *it)[idim];
        ratioG += det.psiMinv(i, WorkingIndex) * det.dpsiM(WorkingIndex, *it)[idim];
        it++;
      }
      det.grads(det.ReferenceDeterminant, WorkingIndex)[idim] = ratioG;
      InverseUpdateByColumn(det.dpsiMinv, det.psiV_temp, det.workV1, det.workV2, WorkingIndex, ratioG);
      for (size_t i = 0; i < det_leader.NumOrbitals; i++)
        det.TpsiM(i, WorkingIndex) = det.dpsiM(WorkingIndex, i)[idim];
      det.BuildDotProductsAndCalculateRatios(det.ReferenceDeterminant, WorkingIndex, det.grads, det.dpsiMinv, det.TpsiM,
                                             det.dotProducts, *det.detData, *det.uniquePairs, *det.DetSigns, idim);
    }

    for (size_t i = 0; i < det_leader.NumOrbitals; i++)
      det.TpsiM(i, WorkingIndex) = det.psiM(WorkingIndex, i);
  }
}

} // namespace qmcplusplus
