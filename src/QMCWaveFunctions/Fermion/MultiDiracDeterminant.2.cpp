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

void MultiDiracDeterminant::mw_BuildDotProductsAndCalculateRatios_impl(int nw,
                                                                       int ref,
                                                                       const std::vector<ValueType>& det0_list,
                                                                       const RefVector<ValueMatrix>& psiinv_list,
                                                                       const RefVector<ValueMatrix>& psi_list,
                                                                       const std::vector<int>& data,
                                                                       const std::vector<std::pair<int, int>>& pairs,
                                                                       const std::vector<RealType>& sign,
                                                                       const RefVector<ValueMatrix>& dotProducts_list,
                                                                       const RefVector<ValueVector>& ratios_list)
{
  const size_t npairs = pairs.size();
  //This is not sure but I think it is the case dur to the use of a const...
  const size_t num                      = psi_list[0].get().extent(1);
  const std::pair<int, int>* restrict p = pairs.data();
  const size_t nitems                   = sign.size();


  readMatTimer.start();
  ///To be flattned with Nb_unique_dets*NW. Needs reorg by excitation for mempry access by stride.
  for (size_t iw = 0; iw < nw; iw++)
  {
    for (size_t i = 0; i < npairs; ++i)
    {
      const int I                      = p[i].first;
      const int J                      = p[i].second;
      dotProducts_list[iw].get()(I, J) = simd::dot(psiinv_list[iw].get()[I], psi_list[iw].get()[J], num);
    }

    std::vector<int>::const_iterator it2 = data.begin();

    // explore Inclusive Scan for OpenMP
    for (size_t count = 0; count < nitems; ++count)
    {
      const size_t n = *it2;
      //ratios[count]=(count!=ref)?sign[count]*det0*CalculateRatioFromMatrixElements(n,dotProducts,it2+1):det0;
      if (count != ref)
      {
        ratios_list[iw].get()[count] =
            sign[count] * det0_list[iw] * CalculateRatioFromMatrixElements(n, dotProducts_list[iw].get(), it2 + 1);
      }
      it2 += 3 * n + 1;
    }
    ratios_list[iw].get()[ref] = det0_list[iw];
  }
  readMatTimer.stop();
}


void MultiDiracDeterminant::BuildDotProductsAndCalculateRatios(int ref,
                                                               const ValueMatrix& psiinv,
                                                               const ValueMatrix& psi,
                                                               const std::vector<int>& data,
                                                               const std::vector<std::pair<int, int>>& pairs,
                                                               const std::vector<RealType>& sign,
                                                               ValueMatrix& dotProducts,
                                                               ValueVector& ratios)
{
  BuildDotProductsAndCalculateRatios_impl(ref, ValueType(1), ratios.data(), psiinv, psi, dotProducts, data, pairs,
                                          sign);
}

void MultiDiracDeterminant::mw_BuildDotProductsAndCalculateRatios(int nw,
                                                                  int ref,
                                                                  const std::vector<ValueType>& det0_list,
                                                                  const RefVector<ValueMatrix>& psiinv_list,
                                                                  const RefVector<ValueMatrix>& psi_list,
                                                                  const std::vector<int>& data,
                                                                  const std::vector<std::pair<int, int>>& pairs,
                                                                  const std::vector<RealType>& sign,
                                                                  const RefVector<ValueMatrix>& dotProducts_list,
                                                                  const RefVector<ValueVector>& ratios_list)
{
  mw_BuildDotProductsAndCalculateRatios_impl(nw, ref, det0_list, psiinv_list, psi_list, data, pairs, sign,
                                             dotProducts_list, ratios_list);
}

void MultiDiracDeterminant::BuildDotProductsAndCalculateRatiosGrads(int ref,
                                                                    const ValueMatrix& psiinv,
                                                                    const ValueMatrix& psi,
                                                                    const std::vector<int>& data,
                                                                    const std::vector<std::pair<int, int>>& pairs,
                                                                    const std::vector<RealType>& sign,
                                                                    const ValueType& det0_grad,
                                                                    ValueMatrix& dotProducts,
                                                                    int dx,
                                                                    int iat,
                                                                    GradMatrix& grads)
{
  BuildDotProductsAndCalculateRatios_impl(ref, det0_grad, WorkSpace.data(), psiinv, psi, dotProducts, data, pairs,
                                          sign);
  for (size_t count = 0; count < getNumDets(); ++count)
    grads(count, iat)[dx] = WorkSpace[count];
}

void MultiDiracDeterminant::mw_BuildDotProductsAndCalculateRatiosGrads(int nw,
                                                                       int ref,
                                                                       int iat,
                                                                       int dx,
                                                                       int getNumDets,
                                                                       const std::vector<ValueType>& det0_grad_list,
                                                                       const RefVector<ValueMatrix>& psiinv_list,
                                                                       const RefVector<ValueMatrix>& psi_list,
                                                                       const std::vector<int>& data,
                                                                       const std::vector<std::pair<int, int>>& pairs,
                                                                       const std::vector<RealType>& sign,
                                                                       const RefVector<ValueVector>& WorkSpace_list,
                                                                       const RefVector<ValueMatrix>& dotProducts_list,
                                                                       const RefVector<GradMatrix>& grads_list)

{
  mw_BuildDotProductsAndCalculateRatios_impl(nw, ref, det0_grad_list, psiinv_list, psi_list, data, pairs, sign,
                                             dotProducts_list, WorkSpace_list);
  for (size_t iw = 0; iw < nw; iw++)
    for (size_t count = 0; count < getNumDets; ++count)
      grads_list[iw].get()(count, iat)[dx] = WorkSpace_list[iw].get()[count];
}

void MultiDiracDeterminant::BuildDotProductsAndCalculateRatiosValueMatrixOneParticle(
    int ref,
    const ValueMatrix& psiinv,
    const ValueMatrix& psi,
    const std::vector<int>& data,
    const std::vector<std::pair<int, int>>& pairs,
    const std::vector<RealType>& sign,
    ValueMatrix& dotProducts,
    int iat,
    ValueMatrix& ratios)
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
  RefVectorWithLeader<SPOSet> phi_list(*det_leader.getPhi());
  det_leader.RatioTimer.start();

  std::vector<ValueType> curRatio_list;
  std::vector<ValueType> det0_list(nw, 1.0);


  RefVector<ValueVector> psiV_list, psiV_temp_list, new_ratios_to_ref_list, workV1_list, workV2_list;
  RefVector<ValueMatrix> psiMinv_temp_list, psiMinv_list, TpsiM_list, dotProducts_list, psiM_list;


  phi_list.reserve(nw);
  psiV_list.reserve(nw);
  psiV_temp_list.reserve(nw);
  psiMinv_list.reserve(nw);
  psiMinv_temp_list.reserve(nw);
  workV1_list.reserve(nw);
  dotProducts_list.reserve(nw);
  workV2_list.reserve(nw);
  TpsiM_list.reserve(nw);
  psiM_list.reserve(nw);
  new_ratios_to_ref_list.reserve(nw);

  curRatio_list.resize(nw);


  for (size_t iw = 0; iw < nw; iw++)
  {
    MultiDiracDeterminant& det = (det_list[iw]);
    det.UpdateMode             = ORB_PBYP_PARTIAL;
    phi_list.push_back(*det.Phi);
    psiV_list.push_back(det.psiV);
    psiV_temp_list.push_back(det.psiV_temp);
    workV1_list.push_back(det.workV1);
    workV2_list.push_back(det.workV2);
    psiMinv_list.push_back(det.psiMinv);
    psiM_list.push_back(det.psiM);
    psiMinv_temp_list.push_back(det.psiMinv_temp);
    new_ratios_to_ref_list.push_back(det.new_ratios_to_ref_);
    dotProducts_list.push_back(det.dotProducts);
    TpsiM_list.push_back(det.TpsiM);
  }

  det_leader.UpdateMode  = ORB_PBYP_RATIO;
  const int WorkingIndex = iat - det_leader.FirstIndex;

  det_leader.evalOrbTimer.start();
  det_leader.getPhi()->mw_evaluateValue(phi_list, P_list, iat, psiV_list);
  det_leader.evalOrbTimer.stop();

  det_leader.ExtraStuffTimer.start();
  for (size_t iw = 0; iw < nw; iw++)
  {
    MultiDiracDeterminant& det = (det_list[iw]);

    const auto& confgList = *det.ciConfigList;
    auto it(confgList[det.ReferenceDeterminant].occup.begin());
    psiMinv_temp_list[iw].get()                                       = psiMinv_list[iw].get();
    new_ratios_to_ref_list[iw].get()[det_leader.ReferenceDeterminant] = ValueType(1);
    for (size_t i = 0; i < det_leader.NumPtcls; i++)
      psiV_temp_list[iw].get()[i] = psiV_list[iw].get()[*(it++)];
  }

  mw_DetRatioByColumn(nw, WorkingIndex, psiMinv_temp_list, psiV_temp_list, curRatio_list);
  mw_InverseUpdateByColumn(nw, psiMinv_temp_list, psiV_temp_list, workV1_list, workV2_list, WorkingIndex,
                           curRatio_list);


  for (size_t iw = 0; iw < nw; iw++)
    for (size_t i = 0; i < det_leader.NumOrbitals; i++)
      TpsiM_list[iw].get()(i, WorkingIndex) = psiV_list[iw].get()[i];


  det_leader.mw_BuildDotProductsAndCalculateRatios(nw, det_leader.ReferenceDeterminant, det0_list, psiMinv_temp_list,
                                                   TpsiM_list, *det_leader.detData, *det_leader.uniquePairs,
                                                   *det_leader.DetSigns, dotProducts_list, new_ratios_to_ref_list);

  det_leader.ExtraStuffTimer.stop();
  for (size_t iw = 0; iw < nw; iw++)
  {
    MultiDiracDeterminant& det = (det_list[iw]);
    det.curRatio               = curRatio_list[iw];
    for (size_t i = 0; i < det_leader.NumOrbitals; i++)
      TpsiM_list[iw].get()(i, WorkingIndex) = psiM_list[iw].get()(WorkingIndex, i);
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
  auto ratio_old_ref_det = DetRatioByColumn(psiMinv_temp, psiV_temp, WorkingIndex);
  curRatio               = ratio_old_ref_det;
  InverseUpdateByColumn(psiMinv_temp, psiV_temp, workV1, workV2, WorkingIndex, ratio_old_ref_det);
  for (size_t i = 0; i < NumOrbitals; i++)
    TpsiM(i, WorkingIndex) = psiV[i];
  ExtraStuffTimer.stop();
  BuildDotProductsAndCalculateRatios(ReferenceDeterminant, psiMinv_temp, TpsiM, *detData, *uniquePairs, *DetSigns,
                                     dotProducts, new_ratios_to_ref_);
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
  curRatio = DetRatioByColumn(psiMinv_temp, psiV_temp, WorkingIndex);
  InverseUpdateByColumn(psiMinv_temp, psiV_temp, workV1, workV2, WorkingIndex, curRatio);
  for (size_t i = 0; i < NumOrbitals; i++)
    TpsiM(i, WorkingIndex) = psiV[i];
  ExtraStuffTimer.stop();
  BuildDotProductsAndCalculateRatios(ReferenceDeterminant, psiMinv_temp, TpsiM, *detData, *uniquePairs, *DetSigns,
                                     dotProducts, new_ratios_to_ref_);
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
    BuildDotProductsAndCalculateRatiosGrads(ReferenceDeterminant, dpsiMinv, TpsiM, *detData, *uniquePairs, *DetSigns,
                                            ratioGradRef[idim] / curRatio, dotProducts, idim, WorkingIndex, new_grads);
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
  new_spingrads(ReferenceDeterminant, WorkingIndex) = ratioSpinGradRef / curRatio;
  InverseUpdateByColumn(psiMinv_temp, psiV_temp, workV1, workV2, WorkingIndex, curRatio);
  for (size_t i = 0; i < NumOrbitals; i++)
    TpsiM(i, WorkingIndex) = psiV[i];
  ExtraStuffTimer.stop();
  BuildDotProductsAndCalculateRatios(ReferenceDeterminant, psiMinv_temp, TpsiM, *detData, *uniquePairs, *DetSigns,
                                     dotProducts, new_ratios_to_ref_);
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
    BuildDotProductsAndCalculateRatiosGrads(ReferenceDeterminant, dpsiMinv, TpsiM, *detData, *uniquePairs, *DetSigns,
                                            ratioGradRef[idim] / curRatio, dotProducts, idim, WorkingIndex, new_grads);
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
  BuildDotProductsAndCalculateRatiosValueMatrixOneParticle(ReferenceDeterminant, dpsiMinv, TpsiM, *detData,
                                                           *uniquePairs, *DetSigns, dotProducts, WorkingIndex,
                                                           new_spingrads);

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
  RefVectorWithLeader<SPOSet> phi_list(*det_leader.getPhi());

  RefVector<ValueVector> psiV_list, psiV_temp_list, d2psiV_list, workV1_list, workV2_list, new_ratios_to_ref_list;
  RefVector<GradVector> dpsiV_list;
  RefVector<GradMatrix> new_grads_list;
  RefVector<ValueMatrix> psiMinv_temp_list, psiMinv_list, dpsiMinv_list, TpsiM_list, dotProducts_list, psiM_list;
  RefVector<ValueVector> WorkSpace_list;

  std::vector<ValueType> curRatio_list, det0_grad_list;
  std::vector<ValueType> det0_list(nw, 1.0);
  std::vector<GradType> ratioGradRef_list;
  std::vector<ValueType> ratioGradReflistIdim;


  phi_list.reserve(nw);
  psiV_list.reserve(nw);
  dpsiV_list.reserve(nw);
  d2psiV_list.reserve(nw);
  workV1_list.reserve(nw);
  workV2_list.reserve(nw);
  psiV_temp_list.reserve(nw);
  psiMinv_temp_list.reserve(nw);
  psiMinv_list.reserve(nw);
  psiM_list.reserve(nw);

  TpsiM_list.reserve(nw);
  new_ratios_to_ref_list.reserve(nw);
  new_grads_list.reserve(nw);
  dotProducts_list.reserve(nw);
  dpsiMinv_list.reserve(nw);
  WorkSpace_list.reserve(nw);

  curRatio_list.resize(nw);
  ratioGradReflistIdim.resize(nw);
  ratioGradRef_list.resize(nw);
  det0_grad_list.resize(nw);

  det_leader.UpdateMode  = ORB_PBYP_PARTIAL;
  const int WorkingIndex = iat - det_leader.FirstIndex;


  for (size_t iw = 0; iw < nw; iw++)
  {
    MultiDiracDeterminant& det = (det_list[iw]);
    det.UpdateMode             = ORB_PBYP_PARTIAL;
    phi_list.push_back(*det.Phi);
    psiV_list.push_back(det.psiV);
    dpsiV_list.push_back(det.dpsiV);
    d2psiV_list.push_back(det.d2psiV);
    workV1_list.push_back(det.workV1);
    workV2_list.push_back(det.workV2);
    psiV_temp_list.push_back(det.psiV_temp);
    psiMinv_list.push_back(det.psiMinv);
    psiM_list.push_back(det.psiM);
    psiMinv_temp_list.push_back(det.psiMinv_temp);
    new_ratios_to_ref_list.push_back(det.new_ratios_to_ref_);
    new_grads_list.push_back(det.new_grads);
    TpsiM_list.push_back(det.TpsiM);
    dotProducts_list.push_back(det.dotProducts);
    dpsiMinv_list.push_back(det.dpsiMinv);
    WorkSpace_list.push_back(det.WorkSpace);
  }

  det_leader.evalOrb1Timer.start();
  ///Should be optimized for real Batched + Offload
  det_leader.getPhi()->mw_evaluateVGL(phi_list, P_list, iat, psiV_list, dpsiV_list, d2psiV_list);
  det_leader.evalOrb1Timer.stop();

  for (size_t iw = 0; iw < nw; iw++)
  {
    MultiDiracDeterminant& det = (det_list[iw]);
    det_leader.ExtraStuffTimer.start();
    psiMinv_temp_list[iw].get() = psiMinv_list[iw].get();
    const auto& confgList       = *det.ciConfigList;
    auto it(confgList[det_leader.ReferenceDeterminant].occup.begin());

    for (size_t i = 0; i < det_leader.NumPtcls; i++)
    {
      psiV_temp_list[iw].get()[i] = psiV_list[iw].get()[*it];
      ratioGradRef_list[iw] += psiMinv_temp_list[iw].get()(i, WorkingIndex) * dpsiV_list[iw].get()[*it];
      it++;
    }
    psiV_temp_list.push_back(psiV_temp_list[iw].get());
    det_leader.ExtraStuffTimer.stop();
  }

  mw_DetRatioByColumn(nw, WorkingIndex, psiMinv_temp_list, psiV_temp_list, curRatio_list);

  for (size_t iw = 0; iw < nw; iw++)
    for (size_t i = 0; i < det_leader.NumOrbitals; i++)
      TpsiM_list[iw].get()(i, WorkingIndex) = psiV_list[iw].get()[i];


  mw_InverseUpdateByColumn(nw, psiMinv_temp_list, psiV_temp_list, workV1_list, workV2_list, WorkingIndex,
                           curRatio_list);
  det_leader.mw_BuildDotProductsAndCalculateRatios(nw, det_leader.ReferenceDeterminant, det0_list, psiMinv_temp_list,
                                                   TpsiM_list, *det_leader.detData, *det_leader.uniquePairs,
                                                   *det_leader.DetSigns, dotProducts_list, new_ratios_to_ref_list);

  for (size_t idim = 0; idim < OHMMS_DIM; idim++)
  {
    for (size_t iw = 0; iw < nw; iw++)
    {
      MultiDiracDeterminant& det = (det_list[iw]);
      ratioGradReflistIdim[iw]   = ratioGradRef_list[iw][idim];
      const auto& confgList      = *det.ciConfigList;
      auto it(confgList[det_leader.ReferenceDeterminant].occup.begin());
      //ExtraStuffTimer.start();
      dpsiMinv_list[iw].get() = psiMinv_list[iw].get();
      it                      = confgList[det_leader.ReferenceDeterminant].occup.begin();
      for (size_t i = 0; i < det_leader.NumPtcls; i++)
        psiV_temp_list[iw].get()[i] = dpsiV_list[iw].get()[*(it++)][idim];
    }
    mw_InverseUpdateByColumn(nw, dpsiMinv_list, psiV_temp_list, workV1_list, workV2_list, WorkingIndex,
                             ratioGradReflistIdim);
    for (size_t iw = 0; iw < nw; iw++)
    {
      for (size_t i = 0; i < det_leader.NumOrbitals; i++)
        TpsiM_list[iw].get()(i, WorkingIndex) = dpsiV_list[iw].get()[i][idim];
      det0_grad_list[iw] = ratioGradReflistIdim[iw] / curRatio_list[iw];
    }
    det_leader.mw_BuildDotProductsAndCalculateRatiosGrads(nw, det_leader.ReferenceDeterminant, WorkingIndex, idim,
                                                          det_leader.getNumDets(), det0_grad_list, dpsiMinv_list,
                                                          TpsiM_list, *det_leader.detData, *det_leader.uniquePairs,
                                                          *det_leader.DetSigns, WorkSpace_list, dotProducts_list,
                                                          new_grads_list);

    for (size_t iw = 0; iw < nw; iw++)
      for (int i = 0; i < det_leader.NumOrbitals; i++)
        TpsiM_list[iw].get()(i, WorkingIndex) = psiM_list[iw].get()(WorkingIndex, i);
  }


  for (size_t iw = 0; iw < nw; iw++)
  {
    MultiDiracDeterminant& det = (det_list[iw]);
    det.curRatio               = curRatio_list[iw];
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
    InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, WorkingIndex, ratioG);
    for (size_t i = 0; i < NumOrbitals; i++)
      TpsiM(i, WorkingIndex) = dpsiM(WorkingIndex, i)[idim];
    BuildDotProductsAndCalculateRatiosGrads(ReferenceDeterminant, dpsiMinv, TpsiM, *detData, *uniquePairs, *DetSigns,
                                            ratioG, dotProducts, idim, WorkingIndex, grads);
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
    InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, WorkingIndex, ratioG);
    for (size_t i = 0; i < NumOrbitals; i++)
      TpsiM(i, WorkingIndex) = dpsiM(WorkingIndex, i)[idim];
    BuildDotProductsAndCalculateRatiosGrads(ReferenceDeterminant, dpsiMinv, TpsiM, *detData, *uniquePairs, *DetSigns,
                                            ratioG, dotProducts, idim, WorkingIndex, grads);
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
  BuildDotProductsAndCalculateRatiosValueMatrixOneParticle(ReferenceDeterminant, dpsiMinv, TpsiM, *detData,
                                                           *uniquePairs, *DetSigns, dotProducts, WorkingIndex,
                                                           spingrads);

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
  const int WorkingIndex            = iat - det_leader.FirstIndex;


  RefVector<ValueVector> psiV_temp_list, workV1_list, workV2_list;
  RefVector<ValueMatrix> dpsiMinv_list, psiMinv_list, psiM_list, TpsiM_list, dotProducts_list;
  RefVector<GradMatrix> dpsiM_list;
  RefVector<GradMatrix> grads_list;
  RefVector<ValueVector> WorkSpace_list;
  std::vector<ValueType> ratioG_list;


  psiMinv_list.reserve(nw);
  dpsiMinv_list.reserve(nw);
  workV1_list.reserve(nw);
  workV2_list.reserve(nw);
  dpsiM_list.reserve(nw);
  psiV_temp_list.reserve(nw);
  grads_list.reserve(nw);
  dotProducts_list.reserve(nw);
  TpsiM_list.reserve(nw);
  psiM_list.reserve(nw);
  WorkSpace_list.reserve(nw);

  ratioG_list.resize(nw);

  for (size_t iw = 0; iw < nw; iw++)
  {
    MultiDiracDeterminant& det = (det_list[iw]);
    psiMinv_list.push_back(det.psiMinv);
    dpsiMinv_list.push_back(det.dpsiMinv);
    psiV_temp_list.push_back(det.psiV_temp);
    workV1_list.push_back(det.workV1);
    workV2_list.push_back(det.workV2);
    dpsiM_list.push_back(det.dpsiM);
    grads_list.push_back(det.grads);
    TpsiM_list.push_back(det.TpsiM);
    psiM_list.push_back(det.psiM);
    dotProducts_list.push_back(det.dotProducts);
    WorkSpace_list.push_back(det.WorkSpace);
  }

  for (size_t idim = 0; idim < OHMMS_DIM; idim++)
  {
    for (size_t iw = 0; iw < nw; iw++)
    {
      MultiDiracDeterminant& det = (det_list[iw]);
      const auto& confgList      = *det.ciConfigList;
      auto it                    = confgList[det_leader.ReferenceDeterminant].occup.begin();

      dpsiMinv_list[iw].get() = psiMinv_list[iw].get();

      ratioG_list[iw] = 0.0;
      for (size_t i = 0; i < det_leader.NumPtcls; i++)
      {
        psiV_temp_list[iw].get()[i] = dpsiM_list[iw].get()(WorkingIndex, *it)[idim];
        ratioG_list[iw] += psiMinv_list[iw].get()(i, WorkingIndex) * dpsiM_list[iw].get()(WorkingIndex, *it)[idim];
        it++;
      }
    }

    mw_InverseUpdateByColumn(nw, dpsiMinv_list, psiV_temp_list, workV1_list, workV2_list, WorkingIndex, ratioG_list);

    for (size_t iw = 0; iw < nw; iw++)
      for (size_t i = 0; i < det_leader.NumOrbitals; i++)
        TpsiM_list[iw].get()(i, WorkingIndex) = dpsiM_list[iw].get()(WorkingIndex, i)[idim];
    det_leader.mw_BuildDotProductsAndCalculateRatiosGrads(nw, det_leader.ReferenceDeterminant, WorkingIndex, idim,
                                                          det_leader.getNumDets(), ratioG_list, dpsiMinv_list,
                                                          TpsiM_list, *det_leader.detData, *det_leader.uniquePairs,
                                                          *det_leader.DetSigns, WorkSpace_list, dotProducts_list,
                                                          grads_list);

    for (size_t iw = 0; iw < nw; iw++)
      for (size_t i = 0; i < det_leader.NumOrbitals; i++)
        TpsiM_list[iw].get()(i, WorkingIndex) = psiM_list[iw].get()(WorkingIndex, i);
  }
}

} // namespace qmcplusplus
