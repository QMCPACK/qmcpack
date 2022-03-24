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
#include "OMPTarget/ompBLAS.hpp"

namespace qmcplusplus
{
/** shared function used by BuildDotProductsAndCalculateRatios */
void MultiDiracDeterminant::BuildDotProductsAndCalculateRatios_impl(
    int ref,
    ValueType det0,
    ValueType* restrict ratios,
    const OffloadMatrix<ValueType>& psiinv,
    const OffloadMatrix<ValueType>& psi,
    OffloadMatrix<ValueType>& dotProducts,
    const OffloadVector<int>& data,
    const VectorSoaContainer<int, 2, OffloadPinnedAllocator<int>>& pairs,
    const OffloadVector<RealType>& sign)
{
  buildTableTimer.start();
  const size_t num     = psi.extent(1);
  const size_t npairs  = pairs.size();
  const size_t nb_cols = dotProducts.cols();
  //MatrixOperators::product_ABt(psiinv,psi,dotProducts);
  const int* first  = pairs.data(0);
  const int* second = pairs.data(1);
  for (size_t i = 0; i < npairs; ++i)
  {
    const int I       = first[i];
    const int J       = second[i];
    dotProducts(I, J) = simd::dot(psiinv[I], psi[J], num);
  }
  dotProducts.updateTo();
  buildTableTimer.stop();
  readMatTimer.start();
  const int* it2      = data.data();
  const size_t nitems = sign.size();
  // explore Inclusive Scan for OpenMP
  for (size_t count = 0; count < nitems; ++count)
  {
    const size_t n = *it2;
    if (count != ref)
      ratios[count] = sign[count] * det0 *
          (n > MaxSmallDet ? det_calculator_.evaluate(dotProducts, it2 + 1, n)
                           : calcSmallDeterminant(n, dotProducts.data(), it2 + 1, nb_cols));
    it2 += 3 * n + 1;
  }

  ratios[ref] = det0;
  readMatTimer.stop();
}

void MultiDiracDeterminant::mw_BuildDotProductsAndCalculateRatios_impl(
    int nw,
    int ref,
    const OffloadVector<ValueType>& det0_list,
    const RefVector<OffloadMatrix<ValueType>>& psiinv_list,
    const RefVector<OffloadMatrix<ValueType>>& psi_list,
    const OffloadVector<int>& data,
    const VectorSoaContainer<int, 2, OffloadPinnedAllocator<int>>& pairs,
    const OffloadVector<RealType>& sign,
    const RefVector<OffloadMatrix<ValueType>>& dotProducts_list,
    const RefVector<OffloadVector<ValueType>>& ratios_list)
{
  const size_t npairs = pairs.size();
  const size_t num    = psi_list[0].get().extent(1);
  const size_t nitems = sign.size();

  const int* first  = pairs.data(0);
  const int* second = pairs.data(1);


  OffloadVector<ValueType*> psiinv_deviceptr_list(nw);
  OffloadVector<ValueType*> psi_deviceptr_list(nw);
  OffloadVector<ValueType*> dotProducts_deviceptr_list(nw);

  for (size_t iw = 0; iw < nw; iw++)
  {
	  psiinv_list[iw].get().updateTo();
	  psi_list[iw].get().updateTo();
  }
  for (size_t iw = 0; iw < nw; iw++)
  {
    psiinv_deviceptr_list[iw]      = psiinv_list[iw].get().device_data();
    psi_deviceptr_list[iw]         = psi_list[iw].get().device_data();
    dotProducts_deviceptr_list[iw] = dotProducts_list[iw].get().device_data();
  }

  readMatTimer.start();

  const size_t nb_cols_psi(psi_list[0].get().cols());
  const size_t nb_cols_psiinv(psiinv_list[0].get().cols());
  const size_t nb_cols_dotProd(dotProducts_list[0].get().cols());

  auto* dotProducts_list_ptr  = dotProducts_deviceptr_list.data();
  const auto* psiinv_list_ptr = psiinv_deviceptr_list.device_data();
  const auto* psi_list_ptr    = psi_deviceptr_list.device_data();
  psiinv_deviceptr_list.updateTo();
  psi_deviceptr_list.updateTo();


  PRAGMA_OFFLOAD("omp target teams distribute  map(always,to: dotProducts_list_ptr[:nw]) \
              is_device_ptr(psi_list_ptr,psiinv_list_ptr)")
  for (size_t iw = 0; iw < nw; iw++)
    for (size_t i = 0; i < npairs; ++i)
    {
      const int I = first[i];
      const int J = second[i];

      ValueType dotProducts_local = 0.0;
      PRAGMA_OFFLOAD("omp parallel for reduction(+ : dotProducts_local)")
      for (size_t ind = 0; ind < num; ind++)
        dotProducts_local += psiinv_list_ptr[iw][I * nb_cols_psiinv + ind] * psi_list_ptr[iw][J * nb_cols_psi + ind];
      dotProducts_list_ptr[iw][I * nb_cols_dotProd + J] = dotProducts_local;
    }

  const int max_ext_level = ndets_per_excitation_level_->size() - 1;

  // Compute workload changes drastically as the excitation level increases.
  // this may need different parallelization strategy.
  size_t det_offset  = 1;
  size_t data_offset = 1;

  auto update_offsets = [&](size_t ext_level) {
    det_offset += (*ndets_per_excitation_level_)[ext_level];
    data_offset += (*ndets_per_excitation_level_)[ext_level] * (3 * ext_level + 1);
  };

  if (max_ext_level >= 1)
  {
    mw_updateRatios<1>(det_offset, data_offset, ratios_list, data, sign, det0_list, dotProducts_list);
    update_offsets(1);
  }

  if (max_ext_level >= 2)
  {
    mw_updateRatios<2>(det_offset, data_offset, ratios_list, data, sign, det0_list, dotProducts_list);
    update_offsets(2);
  }

  if (max_ext_level >= 3)
  {
    mw_updateRatios<3>(det_offset, data_offset, ratios_list, data, sign, det0_list, dotProducts_list);
    update_offsets(3);
  }

  if (max_ext_level >= 4)
  {
    mw_updateRatios<4>(det_offset, data_offset, ratios_list, data, sign, det0_list, dotProducts_list);
    update_offsets(4);
  }

  if (max_ext_level >= 5)
  {
    mw_updateRatios<5>(det_offset, data_offset, ratios_list, data, sign, det0_list, dotProducts_list);
    update_offsets(5);
  }

  if (max_ext_level >= 6)
  {
    for (size_t iw = 0; iw < nw; iw++)
      dotProducts_list[iw].get().updateFrom();
  }
  for (size_t ext_level = 6; ext_level <= max_ext_level; ext_level++)
  {
    mw_updateRatios_generic(ext_level, det_offset, data_offset, ratios_list, det_calculator_, data, sign, det0_list,
                            dotProducts_list);
    update_offsets(ext_level);
  }

  for (size_t iw = 0; iw < nw; iw++)
    ratios_list[iw].get().updateFrom();

  readMatTimer.stop();
}


void MultiDiracDeterminant::BuildDotProductsAndCalculateRatios(
    int ref,
    const OffloadMatrix<ValueType>& psiinv,
    const OffloadMatrix<ValueType>& psi,
    const OffloadVector<int>& data,
    const VectorSoaContainer<int, 2, OffloadPinnedAllocator<int>>& pairs,
    const OffloadVector<RealType>& sign,
    OffloadMatrix<ValueType>& dotProducts,
    OffloadVector<ValueType>& ratios)
{
  BuildDotProductsAndCalculateRatios_impl(ref, ValueType(1), ratios.data(), psiinv, psi, dotProducts, data, pairs,
                                          sign);
}

void MultiDiracDeterminant::mw_BuildDotProductsAndCalculateRatios(
    int nw,
    int ref,
    const OffloadVector<ValueType>& det0_list,
    const RefVector<OffloadMatrix<ValueType>>& psiinv_list,
    const RefVector<OffloadMatrix<ValueType>>& psi_list,
    const OffloadVector<int>& data,
    const VectorSoaContainer<int, 2, OffloadPinnedAllocator<int>>& pairs,
    const OffloadVector<RealType>& sign,
    const RefVector<OffloadMatrix<ValueType>>& dotProducts_list,
    const RefVector<OffloadVector<ValueType>>& ratios_list)
{
  mw_BuildDotProductsAndCalculateRatios_impl(nw, ref, det0_list, psiinv_list, psi_list, data, pairs, sign,
                                             dotProducts_list, ratios_list);
}

void MultiDiracDeterminant::BuildDotProductsAndCalculateRatiosGrads(
    int ref,
    const OffloadMatrix<ValueType>& psiinv,
    const OffloadMatrix<ValueType>& psi,
    const OffloadVector<int>& data,
    const VectorSoaContainer<int, 2, OffloadPinnedAllocator<int>>& pairs,
    const OffloadVector<RealType>& sign,
    const ValueType& det0_grad,
    OffloadMatrix<ValueType>& dotProducts,
    int dx,
    int iat,
    OffloadMatrix<GradType>& grads)
{
  BuildDotProductsAndCalculateRatios_impl(ref, det0_grad, WorkSpace.data(), psiinv, psi, dotProducts, data, pairs,
                                          sign);
  for (size_t count = 0; count < getNumDets(); ++count)
    grads(count, iat)[dx] = WorkSpace[count];
}

void MultiDiracDeterminant::mw_BuildDotProductsAndCalculateRatiosGrads(
    int nw,
    int ref,
    int iat,
    int dx,
    int getNumDets,
    const OffloadVector<ValueType>& det0_grad_list,
    const RefVector<OffloadMatrix<ValueType>>& psiinv_list,
    const RefVector<OffloadMatrix<ValueType>>& psi_list,
    const OffloadVector<int>& data,
    const VectorSoaContainer<int, 2, OffloadPinnedAllocator<int>>& pairs,
    const OffloadVector<RealType>& sign,
    const RefVector<OffloadVector<ValueType>>& WorkSpace_list,
    const RefVector<OffloadMatrix<ValueType>>& dotProducts_list,
    const RefVector<OffloadMatrix<GradType>>& grads_list)

{
  mw_BuildDotProductsAndCalculateRatios_impl(nw, ref, det0_grad_list, psiinv_list, psi_list, data, pairs, sign,
                                             dotProducts_list, WorkSpace_list);

  ///NEEDS TO BE MOVED HERE TO OMP OFFLOAD
  for (size_t iw = 0; iw < nw; iw++)
    for (size_t count = 0; count < getNumDets; ++count)
      grads_list[iw].get()(count, iat)[dx] = WorkSpace_list[iw].get()[count];
}

void MultiDiracDeterminant::BuildDotProductsAndCalculateRatiosValueMatrixOneParticle(
    int ref,
    const OffloadMatrix<ValueType>& psiinv,
    const OffloadMatrix<ValueType>& psi,
    const OffloadVector<int>& data,
    const VectorSoaContainer<int, 2, OffloadPinnedAllocator<int>>& pairs,
    const OffloadVector<RealType>& sign,
    OffloadMatrix<ValueType>& dotProducts,
    int iat,
    OffloadMatrix<ValueType>& ratios)
{
  const ValueType det0 = ratios(ref, iat);
  BuildDotProductsAndCalculateRatios_impl(ref, det0, WorkSpace.data(), psiinv, psi, dotProducts, data, pairs, sign);
  //splatt
  for (size_t count = 0; count < getNumDets(); ++count)
    ratios(count, iat) = WorkSpace[count];
  ratios.updateTo();
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
      ratios(count,iat) = sign[count]*det0*CustomizedMatrixDet(n,dotProducts,it2+1);
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

  OffloadVector<ValueType> det0_list(nw, 1.0);
  OffloadVector<ValueType> curRatio_list(nw, 0.0);
  OffloadVector<size_t> confgListOccup(det_leader.NumPtcls,0.0);

  RefVector<OffloadVector<ValueType>> psiV_list, psiV_temp_list, new_ratios_to_ref_list;
  RefVector<OffloadMatrix<ValueType>> TpsiM_list, psiM_list, dotProducts_list;
  RefVector<OffloadMatrix<ValueType>> psiMinv_temp_list, psiMinv_list;
  RefVector<OffloadVector<ValueType>> workV1_list, workV2_list;


  phi_list.reserve(nw);
  psiV_list.reserve(nw);
  psiV_temp_list.reserve(nw);
  psiMinv_list.reserve(nw);
  psiMinv_temp_list.reserve(nw);
  dotProducts_list.reserve(nw);
  workV1_list.reserve(nw);
  workV2_list.reserve(nw);

  TpsiM_list.reserve(nw);
  psiM_list.reserve(nw);
  new_ratios_to_ref_list.reserve(nw);

  int success = 0;
  int dummy_handle=0;


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
  for (size_t iw = 0; iw < nw; iw++)
  {
    MultiDiracDeterminant& det = (det_list[iw]);
    Vector<ValueType> psiV_list_host_view(psiV_list[iw].get().data(), psiV_list[iw].get().size());
    det.getPhi()->evaluateValue(P_list[iw], iat, psiV_list_host_view);
    ///Big transfer
    psiV_list[iw].get().updateTo();
  }
  det_leader.evalOrbTimer.stop();

  const auto psiMinv_rows   = psiMinv_list[0].get().rows();
  const auto psiMinv_cols   = psiMinv_list[0].get().cols();
  const auto TpsiM_num_cols = TpsiM_list[0].get().cols(); 
  const auto TpsiM_num_rows = TpsiM_list[0].get().rows(); 
  const auto NumPtcls     = det_leader.NumPtcls;
  const auto NumOrbitals     = det_leader.NumOrbitals;
  
  const auto& confgList      = *det_leader.ciConfigList;
  for (size_t i = 0; i < det_leader.NumPtcls; i++)
      confgListOccup[i]=confgList[det_leader.ReferenceDeterminant].occup[i];



  OffloadVector<ValueType*> psiMinv_temp_deviceptr_list(nw);
  OffloadVector<ValueType*> psiMinv_deviceptr_list(nw);
  OffloadVector<ValueType*> psiV_deviceptr_list(nw);
  OffloadVector<ValueType*> TpsiM_deviceptr_list(nw);
  OffloadVector<ValueType*> psiV_temp_deviceptr_list(nw);

  for (size_t iw = 0; iw < nw; iw++)
  {
    psiV_deviceptr_list[iw]         = psiV_list[iw].get().device_data();
    TpsiM_deviceptr_list[iw]         = TpsiM_list[iw].get().device_data();
    psiMinv_deviceptr_list[iw]      = psiMinv_list[iw].get().device_data();
    psiMinv_temp_deviceptr_list[iw]      = psiMinv_temp_list[iw].get().device_data();
    psiV_temp_deviceptr_list[iw]    = psiV_temp_list[iw].get().device_data();
  }

  auto* psiV_list_ptr = psiV_deviceptr_list.device_data();
  auto* TpsiM_list_ptr = TpsiM_deviceptr_list.device_data();
  auto* psiMinv_list_ptr = psiMinv_deviceptr_list.device_data();
  auto* psiMinv_temp_list_ptr = psiMinv_temp_deviceptr_list.device_data();
  auto* psiV_temp_list_ptr    = psiV_temp_deviceptr_list.device_data();
  auto* curRatio_list_ptr = curRatio_list.data();
  auto* confgListOccup_ptr = confgListOccup.device_data();

  psiMinv_deviceptr_list.updateTo();
  psiMinv_temp_deviceptr_list.updateTo();
  TpsiM_deviceptr_list.updateTo();
  psiV_deviceptr_list.updateTo();
  psiV_temp_deviceptr_list.updateTo();
  confgListOccup.updateTo();
  det0_list.updateTo();




  det_leader.ExtraStuffTimer.start();
/*


  success=ompBLAS::copy_batched_offset(dummy_handle, det_leader.NumOrbitals, psiV_list_ptr, 0, 1, TpsiM_list_ptr, WorkingIndex, TpsiM_num_cols, nw);
  if (success != 0)
        throw std::runtime_error("In MultiDiracDeterminant ompBLAS::copy_batched_offset failed.");

 */ 




  success=ompBLAS::copy_batched(dummy_handle, psiMinv_rows*psiMinv_cols, psiMinv_list_ptr,1,psiMinv_temp_list_ptr, 1, nw) ;
  if (success != 0)
        throw std::runtime_error("In MultiDiracDeterminant ompBLAS::copy_batched_offset failed.");

  PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(2) is_device_ptr(psiV_temp_list_ptr, psiV_list_ptr, confgListOccup_ptr)")
  for (size_t iw = 0; iw < nw; iw++)
    for (size_t i = 0; i < NumPtcls; i++)
    {
      size_t J=confgListOccup_ptr[i];
      psiV_temp_list_ptr[iw][i] = psiV_list_ptr[iw][J];
    }

  PRAGMA_OFFLOAD("omp target teams distribute map(always,from:curRatio_list_ptr[:nw]) \
                  is_device_ptr(psiMinv_temp_list_ptr,psiV_temp_list_ptr)")
  for (size_t iw = 0; iw < nw; iw++)
  {
    ValueType c_ratio = 0.0;
    PRAGMA_OFFLOAD("omp parallel for reduction(+ : c_ratio)")
    for (size_t jc = 0; jc < psiMinv_cols; jc +=1)
    {
      size_t ic=jc*psiMinv_cols;
      c_ratio += (psiMinv_temp_list_ptr[iw] + WorkingIndex)[ic] * psiV_temp_list_ptr[iw][jc];
    }
    curRatio_list_ptr[iw] = c_ratio;
  }


  det_leader.omp_mw_InverseUpdateByColumn(nw, WorkingIndex,curRatio_list,psiV_temp_list,workV1_list,workV2_list,psiMinv_temp_list);


  PRAGMA_OFFLOAD("omp target teams distribute parallel for  is_device_ptr(TpsiM_list_ptr, psiV_list_ptr)")
  for (size_t iw = 0; iw < nw; iw++)
    for (size_t i = 0; i < NumOrbitals; i++)
    {
      TpsiM_list_ptr[iw][i*TpsiM_num_cols+WorkingIndex] = psiV_list_ptr[iw][i];
      printf("-1  iw=%zu   i=%zu   Tpsi=%f   psiV=%f \n", iw, i, TpsiM_list_ptr[iw][i*TpsiM_num_cols+WorkingIndex] , psiV_list_ptr[iw][i]);
    }

/*  success=ompBLAS::copy_batched_offset(dummy_handle, det_leader.NumOrbitals, psiV_list_ptr, 0, 1, TpsiM_list_ptr, WorkingIndex, TpsiM_num_cols, nw);
  if (success != 0)
        throw std::runtime_error("In MultiDiracDeterminant ompBLAS::copy_batched_offset failed.");
*/

  for (size_t iw = 0; iw < nw; iw++)
  {
    for (size_t i = 0; i < det_leader.NumOrbitals; i++)
    {
      TpsiM_list[iw].get()(i, WorkingIndex) = psiV_list[iw].get()[i];
      printf("-2  iw=%zu   i=%zu   Tpsi=%f   psiV=%f \n", iw, i, TpsiM_list[iw].get()(i, WorkingIndex) , psiV_list[iw].get()[i]);
    }
    TpsiM_list[iw].get().updateTo();
  }


  for (size_t iw = 0; iw < nw; iw++)
  {
      psiMinv_temp_list[iw].get().updateFrom();
      psiV_temp_list[iw].get().updateFrom();
      TpsiM_list[iw].get().updateFrom();
  }



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

  Vector<ValueType> psiV_host_view(psiV.data(), psiV.size());
  Phi->evaluateValue(P, iat, psiV_host_view);
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
  ///Using Host Views for Phi-evaluateVGL since not ported to GPU
  Vector<ValueType> psiV_host_view(psiV.data(), psiV.size());
  Vector<GradType> dpsiV_host_view(dpsiV.data(), dpsiV.size());
  Vector<ValueType> d2psiV_host_view(d2psiV.data(), d2psiV.size());
  Phi->evaluateVGL(P, iat, psiV_host_view, dpsiV_host_view, d2psiV_host_view);
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
  ///Creating Host view to call Phi->evaluateVGL
  Vector<ValueType> psiV_host_view(psiV.data(), psiV.size());
  Vector<GradType> dpsiV_host_view(dpsiV.data(), dpsiV.size());
  Vector<ValueType> d2psiV_host_view(d2psiV.data(), d2psiV.size());
  Phi->evaluateVGL_spin(P, iat, psiV_host_view, dpsiV_host_view, d2psiV_host_view, dspin_psiV);
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

  const int NumOrbitals=det_leader.NumOrbitals;
  RefVector<OffloadVector<ValueType>> psiV_list, psiV_temp_list, new_ratios_to_ref_list, WorkSpace_list;
  RefVector<OffloadVector<ValueType>> d2psiV_list, workV1_list, workV2_list;

  int success = 0;
  int dummy_handle=0;
  RefVector<OffloadMatrix<ValueType>> psiMinv_temp_list, psiMinv_list, dpsiMinv_list;
  RefVector<OffloadMatrix<ValueType>> dotProducts_list, psiM_list, TpsiM_list;

  RefVector<OffloadVector<GradType>> dpsiV_list;
  RefVector<OffloadMatrix<GradType>> new_grads_list;

  const size_t NumPtcls(det_leader.NumPtcls);
  OffloadVector<size_t> confgListOccup(NumPtcls,0.0);
  OffloadVector<ValueType> curRatio_list(nw,0.0), det0_grad_list, det0_list(nw, 1.0), ratioGradReflistIdim;
  std::vector<GradType> ratioGradRef_list;

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
  for (size_t iw = 0; iw < nw; iw++)
  {
    MultiDiracDeterminant& det = (det_list[iw]);
    Vector<ValueType> psiV_list_host_view(psiV_list[iw].get().data(), psiV_list[iw].get().size());
    Vector<GradType> dpsiV_list_host_view(dpsiV_list[iw].get().data(), dpsiV_list[iw].get().size());
    Vector<ValueType> d2psiV_list_host_view(d2psiV_list[iw].get().data(), d2psiV_list[iw].get().size());
    det.Phi->evaluateVGL(P_list[iw], iat, psiV_list_host_view, dpsiV_list_host_view, d2psiV_list_host_view);
    psiV_list[iw].get().updateTo();
    dpsiV_list[iw].get().updateTo();

  }
  det_leader.evalOrb1Timer.stop();

  const auto psiMinv_rows   = psiMinv_list[0].get().rows();
  const auto psiMinv_cols   = psiMinv_list[0].get().cols();
  const auto TpsiM_num_cols = TpsiM_list[0].get().cols(); 
  const auto& confgList      = *det_leader.ciConfigList;
  for (size_t i = 0; i < det_leader.NumPtcls; i++)
      confgListOccup[i]=confgList[det_leader.ReferenceDeterminant].occup[i];



  OffloadVector<ValueType*> psiV_deviceptr_list(nw);
  OffloadVector<ValueType*> psiV_temp_deviceptr_list(nw);
  OffloadVector<ValueType*> TpsiM_deviceptr_list(nw);
  OffloadVector<ValueType*> psiMinv_temp_deviceptr_list(nw);
  OffloadVector<ValueType*> psiMinv_deviceptr_list(nw);
  OffloadVector<ValueType*> dpsiMinv_deviceptr_list(nw);

  for (size_t iw = 0; iw < nw; iw++)
  {
    psiV_deviceptr_list[iw]         = psiV_list[iw].get().device_data();
    psiV_temp_deviceptr_list[iw]         = psiV_temp_list[iw].get().device_data();
    TpsiM_deviceptr_list[iw]         = TpsiM_list[iw].get().device_data();
    psiMinv_deviceptr_list[iw]      = psiMinv_list[iw].get().device_data();
    dpsiMinv_deviceptr_list[iw] = dpsiMinv_list[iw].get().device_data();
    psiMinv_temp_deviceptr_list[iw]      = psiMinv_temp_list[iw].get().device_data();


  }

  auto* psiV_list_ptr = psiV_deviceptr_list.device_data();
  auto* psiV_temp_list_ptr = psiV_temp_deviceptr_list.device_data();
  auto* TpsiM_list_ptr = TpsiM_deviceptr_list.device_data();
  auto* psiMinv_list_ptr = psiMinv_deviceptr_list.device_data();
  auto* dpsiMinv_list_ptr= dpsiMinv_deviceptr_list.device_data();
  auto* psiMinv_temp_list_ptr = psiMinv_temp_deviceptr_list.device_data();

  auto* curRatio_list_ptr = curRatio_list.data();
  auto* det0_grad_list_ptr = det0_grad_list.data();
  auto* confgListOccup_ptr = confgListOccup.device_data();
  auto* ratioGradReflistIdim_ptr=ratioGradReflistIdim.device_data();

  psiMinv_deviceptr_list.updateTo();
  dpsiMinv_deviceptr_list.updateTo();
  psiMinv_temp_deviceptr_list.updateTo();
  TpsiM_deviceptr_list.updateTo();
  psiV_deviceptr_list.updateTo();
  psiV_temp_deviceptr_list.updateTo();
  confgListOccup.updateTo();
  det0_list.updateTo();


  success=ompBLAS::copy_batched(dummy_handle, psiMinv_rows*psiMinv_cols, psiMinv_list_ptr,1,psiMinv_temp_list_ptr, 1, nw) ;
  if (success != 0)
        throw std::runtime_error("In MultiDiracDeterminant ompBLAS::copy_batched_offset failed.");

  PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(2) is_device_ptr(psiV_temp_list_ptr, psiV_list_ptr, confgListOccup_ptr)")
  for (size_t iw = 0; iw < nw; iw++)
    for (size_t i = 0; i < NumPtcls; i++)
    {
      size_t J=confgListOccup_ptr[i];
      psiV_temp_list_ptr[iw][i] = psiV_list_ptr[iw][J];
    }

  for (size_t iw = 0; iw < nw; iw++)
  {
   psiMinv_temp_list[iw].get().updateFrom();
   psiV_temp_list[iw].get().updateFrom();
  }

  for (size_t iw = 0; iw < nw; iw++)
    for (size_t i = 0; i < NumPtcls; i++)
    {
      size_t J=confgListOccup[i];
      ratioGradRef_list[iw] += psiMinv_temp_list[iw].get()(i, WorkingIndex) * dpsiV_list[iw].get()[J];
    }


  PRAGMA_OFFLOAD("omp target teams distribute map(always,from:curRatio_list_ptr[:nw]) \
                  is_device_ptr(psiV_temp_list_ptr,psiMinv_temp_list_ptr)")
  for (size_t iw = 0; iw < nw; iw++)
  {
    ValueType c_ratio = 0.0;
    PRAGMA_OFFLOAD("omp parallel for reduction(+ : c_ratio)")
    for (size_t jc = 0; jc < psiMinv_cols; jc +=1)
    {
      size_t ic=jc*psiMinv_cols;
      c_ratio += (psiMinv_temp_list_ptr[iw] + WorkingIndex)[ic] * psiV_temp_list_ptr[iw][jc];
    }
    curRatio_list_ptr[iw] = c_ratio;
  }

  for (size_t iw = 0; iw < nw; iw++)
    for (size_t i = 0; i < det_leader.NumOrbitals; i++)
      TpsiM_list[iw].get()(i, WorkingIndex) = psiV_list[iw].get()[i];


  det_leader.omp_mw_InverseUpdateByColumn(nw, WorkingIndex,curRatio_list,psiV_temp_list,workV1_list,workV2_list,psiMinv_temp_list);
  for (size_t iw = 0; iw < nw; iw++)
   psiMinv_temp_list[iw].get().updateFrom();

  det_leader.mw_BuildDotProductsAndCalculateRatios(nw, det_leader.ReferenceDeterminant, det0_list, psiMinv_temp_list,
                                                   TpsiM_list, *det_leader.detData, *det_leader.uniquePairs,
                                                   *det_leader.DetSigns, dotProducts_list, new_ratios_to_ref_list);

  for (size_t idim = 0; idim < OHMMS_DIM; idim++)
  {

    success=ompBLAS::copy_batched(dummy_handle, psiMinv_rows*psiMinv_cols, psiMinv_list_ptr,1,dpsiMinv_list_ptr, 1, nw) ;
    if (success != 0)
          throw std::runtime_error("In MultiDiracDeterminant ompBLAS::copy_batched_offset failed.");

    for (size_t iw = 0; iw < nw; iw++)
    {
      ratioGradReflistIdim[iw]   = ratioGradRef_list[iw][idim];

      for (size_t i = 0; i < NumPtcls; i++)
      {
        size_t J=confgListOccup[i];
        psiV_temp_list[iw].get()[i] = dpsiV_list[iw].get()[J][idim];
      }
    }

    for (size_t iw = 0; iw < nw; iw++)
      psiV_temp_list[iw].get().updateTo();
    ratioGradReflistIdim.updateTo();

    det_leader.omp_mw_InverseUpdateByColumn(nw, WorkingIndex,ratioGradReflistIdim,psiV_temp_list,workV1_list,workV2_list,dpsiMinv_list);



    for (size_t iw = 0; iw < nw; iw++)
    {
      for (size_t i = 0; i < det_leader.NumOrbitals; i++)
        TpsiM_list[iw].get()(i, WorkingIndex) = dpsiV_list[iw].get()[i][idim];
    }

  PRAGMA_OFFLOAD("omp target teams distribute parallel for  map(always,from:det0_grad_list_ptr[:nw]) \
                  map(always, to: curRatio_list_ptr[:nw]) is_device_ptr(ratioGradReflistIdim_ptr)")
  for (size_t iw = 0; iw < nw; iw++)
      det0_grad_list_ptr[iw] = ratioGradReflistIdim_ptr[iw] / curRatio_list_ptr[iw];

    for (size_t iw = 0; iw < nw; iw++){
	    dpsiMinv_list[iw].get().updateFrom();
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

/*
  det_leader.ExtraStuffTimer.start();
  success=ompBLAS::copy_batched(dummy_handle, psiMinv_rows*psiMinv_cols, psiMinv_list_ptr,1,psiMinv_temp_list_ptr, 1, nw) ;
  if (success != 0)
        throw std::runtime_error("In MultiDiracDeterminant ompBLAS::copy_batched_offset failed.");

  PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(2) is_device_ptr(psiV_temp_list_ptr, psiV_list_ptr, confgListOccup_ptr)")
  for (size_t iw = 0; iw < nw; iw++)
    for (size_t i = 0; i < NumPtcls; i++)
    {
      size_t J=confgListOccup_ptr[i];
      psiV_temp_list_ptr[iw][i] = psiV_list_ptr[iw][J];
    }


  psiMinv_temp_deviceptr_list.updateFrom();
  for (size_t iw = 0; iw < nw; iw++)
  {
    ///GradType ratioGradref=(0,0);
    for (size_t i = 0; i < det_leader.NumPtcls; i++)
    {
      size_t J=confgListOccup[i];
      ratioGradRef_list[iw] += psiMinv_temp_list[iw].get()(i, WorkingIndex) * dpsiV_list[iw].get()[J];
    }
    //ratioGradRef_list[iw]=ratioGradref;
  }



  success=ompBLAS::copy_batched_offset(dummy_handle, det_leader.NumOrbitals, psiV_list_ptr, 0, 1, TpsiM_list_ptr, WorkingIndex, TpsiM_num_cols, nw);
  if (success != 0)
        throw std::runtime_error("In MultiDiracDeterminant ompBLAS::copy_batched_offset failed.");


  auto* ratioGradReflistIdim_ptr = ratioGradReflistIdim.data();
  auto* det0_grad_list_ptr  = det0_grad_list.data();

  PRAGMA_OFFLOAD("omp target teams distribute map(always,from:curRatio_list_ptr[:nw]) \
                  map(always, to: psiV_temp_list_ptr[:nw]) is_device_ptr(psiMinv_temp_list_ptr)")
  for (size_t iw = 0; iw < nw; iw++)
  {
    ValueType c_ratio = 0.0;
    PRAGMA_OFFLOAD("omp parallel for reduction(+ : c_ratio)")
    for (size_t jc = 0; jc < psiMinv_cols; jc +=1)
    {
      size_t ic=jc*psiMinv_cols;
      c_ratio += (psiMinv_temp_list_ptr[iw] + WorkingIndex)[ic] * psiV_temp_list_ptr[iw][jc];
    }
    curRatio_list_ptr[iw] = c_ratio;
  }

  det_leader.omp_mw_InverseUpdateByColumn(nw, WorkingIndex,curRatio_list,psiV_temp_list,workV1_list,workV2_list,psiMinv_temp_list);

  det_leader.mw_BuildDotProductsAndCalculateRatios(nw, det_leader.ReferenceDeterminant, det0_list, psiMinv_temp_list,
                                                   TpsiM_list, *det_leader.detData, *det_leader.uniquePairs,
                                                   *det_leader.DetSigns, dotProducts_list, new_ratios_to_ref_list);



  for (size_t idim = 0; idim < OHMMS_DIM; idim++)
  {

    
    success=ompBLAS::copy_batched(dummy_handle, psiMinv_rows*psiMinv_cols, psiMinv_list_ptr,1,dpsiMinv_list_ptr, 1, nw) ;
    if (success != 0)
          throw std::runtime_error("In MultiDiracDeterminant ompBLAS::copy_batched_offset failed.");

    for (size_t iw = 0; iw < nw; iw++)
    {
      ratioGradReflistIdim[iw]   = ratioGradRef_list[iw][idim];
      for (size_t i = 0; i < NumPtcls; i++)
      {
        size_t J=confgListOccup[i];
        psiV_temp_list[iw].get()[i] = dpsiV_list[iw].get()[J][idim];
      }
      psiV_temp_list[iw].get().updateTo();
    }
    ratioGradReflistIdim.updateTo();

    det_leader.omp_mw_InverseUpdateByColumn(nw, WorkingIndex,ratioGradReflistIdim,psiV_temp_list,workV1_list,workV2_list,dpsiMinv_list);


    for (size_t iw = 0; iw < nw; iw++)
    {
      for (size_t i = 0; i < NumOrbitals; i++)
        TpsiM_list[iw](i, WorkingIndex) = dpsiV_list[iw].get()[i][idim];
      det0_grad_list[iw] = ratioGradReflistIdim[iw] / curRatio_list[iw];

      TpsiM_list[iw].get().updateTo();
    }
    det0_grad_list.updateTo();

    det_leader.mw_BuildDotProductsAndCalculateRatiosGrads(nw, det_leader.ReferenceDeterminant, WorkingIndex, idim,
                                                          det_leader.getNumDets(), det0_grad_list, dpsiMinv_list,
                                                          TpsiM_list, *det_leader.detData, *det_leader.uniquePairs,
                                                          *det_leader.DetSigns, WorkSpace_list, dotProducts_list,
                                                          new_grads_list);

    for (size_t iw = 0; iw < nw; iw++)
      for (int i = 0; i < det_leader.NumOrbitals; i++)
        TpsiM_list[iw].get()(i, WorkingIndex) = psiM_list[iw].get()(WorkingIndex, i);
  }
*/

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
  int success=0;
  int dummy_handle=0;

  RefVector<OffloadMatrix<ValueType>> dpsiMinv_list, psiMinv_list;
  RefVector<OffloadMatrix<GradType>> dpsiM_list;
  RefVector<OffloadVector<ValueType>> psiV_temp_list, WorkSpace_list, workV1_list, workV2_list;
  RefVector<OffloadMatrix<ValueType>> dotProducts_list, TpsiM_list, psiM_list;
  RefVector<OffloadMatrix<GradType>> grads_list;

  OffloadVector<ValueType> ratioG_list;

  const size_t NumPtcls(det_leader.NumPtcls);
  OffloadVector<size_t> confgListOccup(NumPtcls,0.0);

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

  const auto psiMinv_rows   = psiMinv_list[0].get().rows();
  const auto psiMinv_cols   = psiMinv_list[0].get().cols();

  const auto& confgList      = *det_leader.ciConfigList;
  for (size_t i = 0; i < det_leader.NumPtcls; i++)
      confgListOccup[i]=confgList[det_leader.ReferenceDeterminant].occup[i];


  OffloadVector<ValueType*> psiMinv_deviceptr_list(nw);
  OffloadVector<ValueType*> dpsiMinv_deviceptr_list(nw);
  for (size_t iw = 0; iw < nw; iw++)
  {
    psiMinv_deviceptr_list[iw]      = psiMinv_list[iw].get().device_data();
    dpsiMinv_deviceptr_list[iw] = dpsiMinv_list[iw].get().device_data();
  }

  auto* psiMinv_list_ptr = psiMinv_deviceptr_list.device_data();
  auto* dpsiMinv_list_ptr= dpsiMinv_deviceptr_list.device_data();

  psiMinv_deviceptr_list.updateTo();
  dpsiMinv_deviceptr_list.updateTo();
  confgListOccup.updateTo();

  for (size_t idim = 0; idim < OHMMS_DIM; idim++)
  {

    success=ompBLAS::copy_batched(dummy_handle, psiMinv_rows*psiMinv_cols, psiMinv_list_ptr,1,dpsiMinv_list_ptr, 1, nw) ;
    if (success != 0)
          throw std::runtime_error("In MultiDiracDeterminant ompBLAS::copy_batched_offset failed.");


    for (size_t iw = 0; iw < nw; iw++)
      dpsiMinv_list[iw].get().updateFrom();
      
    for (size_t iw = 0; iw < nw; iw++)
    {
      ratioG_list[iw] = 0.0;
      for (size_t i = 0; i < det_leader.NumPtcls; i++)
      {
        size_t J=confgListOccup[i];
        psiV_temp_list[iw].get()[i] = dpsiM_list[iw].get()(WorkingIndex, J)[idim];
        ratioG_list[iw] += psiMinv_list[iw].get()(i, WorkingIndex) * dpsiM_list[iw].get()(WorkingIndex, J)[idim];
      }
      psiV_temp_list[iw].get().updateTo();
    }
    ratioG_list.updateTo();

    det_leader.omp_mw_InverseUpdateByColumn(nw, WorkingIndex,ratioG_list,psiV_temp_list,workV1_list,workV2_list,dpsiMinv_list);

    for (size_t iw = 0; iw < nw; iw++)
        dpsiMinv_list[iw].get().updateFrom();


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

/*
  for (size_t idim = 0; idim < OHMMS_DIM; idim++)
  {
    success=ompBLAS::copy_batched(dummy_handle, psiMinv_rows*psiMinv_cols, psiMinv_list_ptr,1,dpsiMinv_list_ptr, 1, nw) ;
    if (success != 0)
          throw std::runtime_error("In MultiDiracDeterminant ompBLAS::copy_batched_offset failed.");



    for (size_t iw = 0; iw < nw; iw++)
    {
      dpsiMinv_list[iw].get().updateFrom();
      ratioG_list[iw] = 0.0;
      for (size_t i = 0; i < NumPtcls; i++)
      {
        size_t J=confgListOccup[i];
        psiV_temp_list[iw].get()[i] = dpsiM_list[iw].get()(WorkingIndex, J)[idim];
        ratioG_list[iw] += psiMinv_list[iw].get()(i, WorkingIndex) * dpsiM_list[iw].get()(WorkingIndex, J)[idim];
      }
      psiV_temp_list[iw].get().updateTo();
    }
    ratioG_list.updateTo();
    det_leader.omp_mw_InverseUpdateByColumn(nw, WorkingIndex,ratioG_list,psiV_temp_list,workV1_list,workV2_list,dpsiMinv_list);


    for (size_t iw = 0; iw < nw; iw++)
      for (size_t i = 0; i < det_leader.NumOrbitals; i++)
        TpsiM_list[iw].get()(i, WorkingIndex) = dpsiM_list[iw].get()(WorkingIndex, i)[idim];

    for (size_t iw = 0; iw < nw; iw++)
                 TpsiM_list[iw].get().updateTo();

    det_leader.mw_BuildDotProductsAndCalculateRatiosGrads(nw, det_leader.ReferenceDeterminant, WorkingIndex, idim,
                                                          det_leader.getNumDets(), ratioG_list, dpsiMinv_list,
                                                          TpsiM_list, *det_leader.detData, *det_leader.uniquePairs,
                                                          *det_leader.DetSigns, WorkSpace_list, dotProducts_list,
                                                          grads_list);

    for (size_t iw = 0; iw < nw; iw++)
      for (size_t i = 0; i < det_leader.NumOrbitals; i++)
        TpsiM_list[iw].get()(i, WorkingIndex) = psiM_list[iw].get()(WorkingIndex, i);
  }
  */
}

void MultiDiracDeterminant::mw_updateRatios_generic(int ext_level,
                                                    const size_t det_offset,
                                                    const size_t data_offset,
                                                    const RefVector<OffloadVector<ValueType>>& ratios_list,
                                                    SmallMatrixDetCalculator<ValueType>& det_calculator,
                                                    const OffloadVector<int>& data,
                                                    const OffloadVector<RealType>& sign,
                                                    const OffloadVector<ValueType>& det0_list,
                                                    const RefVector<OffloadMatrix<ValueType>>& dotProducts_list) const
{
  const size_t nw = ratios_list.size();
  const int* it2  = data.data() + data_offset;
  for (size_t iw = 0; iw < nw; iw++)
    for (size_t count = 0; count < (*ndets_per_excitation_level_)[ext_level]; ++count)
    {
      size_t det_id                 = det_offset + count;
      ratios_list[iw].get()[det_id] = sign[det_id] * det0_list[iw] *
          det_calculator.evaluate(dotProducts_list[iw].get(), it2 + 1 + count * (3 * ext_level + 1), ext_level);
    }
}

template<unsigned EXT_LEVEL>
void MultiDiracDeterminant::mw_updateRatios(const size_t det_offset,
                                            const size_t data_offset,
                                            const RefVector<OffloadVector<ValueType>>& ratios_list,
                                            const OffloadVector<int>& data,
                                            const OffloadVector<RealType>& sign,
                                            const OffloadVector<ValueType>& det0_list,
                                            const RefVector<OffloadMatrix<ValueType>>& dotProducts_list) const
{
  const size_t nw        = ratios_list.size();
  const size_t size_sign = sign.size();
  const size_t nb_cols_dotProd(dotProducts_list[0].get().cols());
  const size_t ndet_ext = (*ndets_per_excitation_level_)[EXT_LEVEL];


  OffloadVector<ValueType*> ratios_deviceptr_list(nw);
  OffloadVector<ValueType*> dotProducts_deviceptr_list(nw);

  for (size_t iw = 0; iw < nw; iw++)
  {
    ratios_deviceptr_list[iw]      = ratios_list[iw].get().device_data();
    dotProducts_deviceptr_list[iw] = dotProducts_list[iw].get().device_data();
  }


  auto* ratios_list_ptr            = ratios_deviceptr_list.data();
  const auto* sign_ptr             = sign.data();
  const int* data_ptr              = data.data();
  const auto* det0_list_ptr        = det0_list.data();
  const auto* dotProducts_list_ptr = dotProducts_deviceptr_list.data();


  PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(2) map(always, to:ratios_list_ptr[:nw]) \
		                                                       map(always, to:dotProducts_list_ptr[:nw])")
  for (size_t iw = 0; iw < nw; iw++)
    for (size_t count = 0; count < ndet_ext; ++count)
    {
      size_t det_id = det_offset + count;
      ValueType ratios_local;
      ///Initialization here to avoid one additional transfer and allow the use of collapse(2)
      ratios_list_ptr[iw][0] = det0_list_ptr[iw];
      ratios_local           = sign_ptr[det_id] * det0_list_ptr[iw] *
          CustomizedMatrixDet<EXT_LEVEL>::evaluate(dotProducts_list_ptr[iw],
                                                   (data_ptr + data_offset) + 1 + count * (3 * EXT_LEVEL + 1),
                                                   nb_cols_dotProd);
      ratios_list_ptr[iw][det_id] = ratios_local;
    }
}


 void MultiDiracDeterminant::omp_mw_InverseUpdateByColumn(int nw,
		                                         const int idx,
                                                         const OffloadVector<ValueType>& curRatio_list,
                                                         const RefVector<OffloadVector<ValueType>>& psiV_list,
                                                         RefVector<OffloadVector<ValueType>>& workV1_list,
                                                         RefVector<OffloadVector<ValueType>>& workV2_list,
                                                         RefVector<OffloadMatrix<ValueType>>& psiMinv_list) const
{



  const ValueType cone(1);
  constexpr ValueType czero(0);
  OffloadVector<ValueType> czero_vec(nw,czero);
  ValueType* czero_ptr = czero_vec.device_data();
  czero_vec.updateTo();

  int success = 0;

  constexpr ValueType cminus_one(-1.0);
  OffloadVector<ValueType> cminus_one_vec(nw,cminus_one);
  ValueType* cminus_one_ptr = cminus_one_vec.device_data();
  cminus_one_vec.updateTo();

  int dummy_handle=0;
  const auto psiMinv_rows = psiMinv_list[0].get().rows();
  OffloadVector<ValueType> invCurRatio_list(nw, 1.0);

  OffloadVector<ValueType*> workV1_deviceptr_list(nw);
  OffloadVector<ValueType*> workV2_deviceptr_list(nw);
  OffloadVector<ValueType*> psiV_deviceptr_list(nw);
  OffloadVector<ValueType*> psiMinv_deviceptr_list(nw);

  for (size_t iw = 0; iw < nw; iw++)
  {
    psiV_deviceptr_list[iw]    = psiV_list[iw].get().device_data();
    psiMinv_deviceptr_list[iw] = psiMinv_list[iw].get().device_data();
    workV1_deviceptr_list[iw] = workV1_list[iw].get().device_data();
    workV2_deviceptr_list[iw] = workV2_list[iw].get().device_data();
  }

  auto* psiV_list_ptr    = psiV_deviceptr_list.device_data();
  auto* psiMinv_list_ptr = psiMinv_deviceptr_list.device_data();
  auto* workV1_list_ptr   = workV1_deviceptr_list.device_data();
  auto* workV2_list_ptr   = workV2_deviceptr_list.device_data();
  auto* curRatio_list_ptr = curRatio_list.data();
  auto* invCurRatio_list_ptr = invCurRatio_list.device_data();


  psiV_deviceptr_list.updateTo();
  psiMinv_deviceptr_list.updateTo();
  workV1_deviceptr_list.updateTo();
  workV2_deviceptr_list.updateTo();

  PRAGMA_OFFLOAD("omp target teams distribute parallel for is_device_ptr(invCurRatio_list_ptr)")
  for (size_t iw = 0; iw < nw; iw++)
    invCurRatio_list_ptr[iw] = cone / curRatio_list_ptr[iw];
     

  success=ompBLAS::gemv_batched(dummy_handle, 'N', psiMinv_rows, psiMinv_rows, invCurRatio_list_ptr, psiMinv_list_ptr, psiMinv_rows, psiV_list_ptr, 1, czero_ptr, workV1_list_ptr, 1,nw);
  if (success != 0)
        throw std::runtime_error("In MultiDiracDeterminant ompBLAS::gemv_batched failed.");

  PRAGMA_OFFLOAD("omp target teams distribute parallel for is_device_ptr(workV1_list_ptr,invCurRatio_list_ptr)")
  for (size_t iw = 0; iw < nw; iw++)
    workV1_list_ptr[iw][idx] = cone - invCurRatio_list_ptr[iw];

  success=ompBLAS::copy_batched_offset(dummy_handle, psiMinv_rows, psiMinv_list_ptr, idx, psiMinv_rows, workV2_list_ptr, 0, 1, nw);
  if (success != 0)
        throw std::runtime_error("In MultiDiracDeterminant ompBLAS::copy_batched_offset failed.");

  success=ompBLAS::ger_batched(dummy_handle,psiMinv_rows, psiMinv_rows, cminus_one_ptr, workV1_list_ptr, 1, workV2_list_ptr, 1, 
             psiMinv_list_ptr, psiMinv_rows,nw);
  if (success != 0)
        throw std::runtime_error("In MultiDiracDeterminant ompBLAS::ger_batched failed.");


}

} // namespace qmcplusplus
