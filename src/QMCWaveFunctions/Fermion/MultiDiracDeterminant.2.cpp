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
#include "OMPTarget/ompReductionComplex.hpp"
#include "OhmmsPETE/ompReductionTinyVector.hpp"

namespace qmcplusplus
{
/** shared function used by buildTableMatrix_calculateRatios */
void MultiDiracDeterminant::buildTableMatrix_calculateRatios_impl(
    int ref,
    ValueType det0,
    ValueType* restrict ratios,
    const OffloadMatrix<ValueType>& psiinv,
    const OffloadMatrix<ValueType>& psi,
    OffloadMatrix<ValueType>& table_matrix,
    const OffloadVector<int>& data,
    const VectorSoaContainer<int, 2, OffloadPinnedAllocator<int>>& pairs,
    const OffloadVector<RealType>& sign)
{
  {
    ScopedTimer local(buildTable_timer);
    const size_t num    = psi.extent(1);
    const size_t npairs = pairs.size();
    //MatrixOperators::product_ABt(psiinv,psi,table_matrix);
    const int* first  = pairs.data(0);
    const int* second = pairs.data(1);
    for (size_t i = 0; i < npairs; ++i)
    {
      const int I        = first[i];
      const int J        = second[i];
      table_matrix(I, J) = simd::dot(psiinv[I], psi[J], num);
    }
  }

  {
    ScopedTimer local(table2ratios_timer);
    const int* it2       = data.data();
    const size_t nitems  = sign.size();
    const size_t nb_cols = table_matrix.cols();
    // explore Inclusive Scan for OpenMP
    for (size_t count = 0; count < nitems; ++count)
    {
      const size_t n = *it2;
      if (count != ref)
        ratios[count] = sign[count] * det0 *
            (n > MaxSmallDet ? det_calculator_.evaluate(table_matrix, it2 + 1, n)
                             : calcSmallDeterminant(n, table_matrix.data(), it2 + 1, nb_cols));
      it2 += 3 * n + 1;
    }

    ratios[ref] = det0;
  }
}

void MultiDiracDeterminant::mw_buildTableMatrix_calculateRatios_impl(
    MultiDiracDetMultiWalkerResource& mw_res,
    int ref,
    const OffloadVector<ValueType>& det0_list,
    const RefVector<OffloadMatrix<ValueType>>& psiinv_list,
    const RefVector<OffloadMatrix<ValueType>>& psi_list,
    const OffloadVector<int>& data,
    const VectorSoaContainer<int, 2, OffloadPinnedAllocator<int>>& pairs,
    const OffloadVector<RealType>& sign,
    const RefVector<OffloadMatrix<ValueType>>& table_matrix_list,
    const RefVector<OffloadVector<ValueType>>& ratios_list)
{
  const size_t nw                   = ratios_list.size();
  auto& psiinv_deviceptr_list       = mw_res.psiinv_deviceptr_list;
  auto& psi_deviceptr_list          = mw_res.psi_deviceptr_list;
  auto& table_matrix_deviceptr_list = mw_res.table_matrix_deviceptr_list;
  auto& ratios_deviceptr_list       = mw_res.ratios_deviceptr_list;
  const size_t nb_cols_table_matrix(table_matrix_list[0].get().cols());

  {
    ScopedTimer local(buildTable_timer);
    const size_t npairs = pairs.size();
    const size_t num    = psi_list[0].get().extent(1);
    const size_t nitems = sign.size();

    const int* first  = pairs.data(0);
    const int* second = pairs.data(1);

    psiinv_deviceptr_list.resize(nw);
    psi_deviceptr_list.resize(nw);
    table_matrix_deviceptr_list.resize(nw);
    ratios_deviceptr_list.resize(nw);

    for (size_t iw = 0; iw < nw; iw++)
    {
      psiinv_deviceptr_list[iw]       = psiinv_list[iw].get().device_data();
      psi_deviceptr_list[iw]          = psi_list[iw].get().device_data();
      table_matrix_deviceptr_list[iw] = table_matrix_list[iw].get().device_data();
      ratios_deviceptr_list[iw]       = ratios_list[iw].get().device_data();
    }

    const size_t nb_cols_psi(psi_list[0].get().cols());
    const size_t nb_cols_psiinv(psiinv_list[0].get().cols());

    auto* ratios_list_ptr       = ratios_deviceptr_list.data();
    auto* table_matrix_list_ptr = table_matrix_deviceptr_list.data();
    const auto* psiinv_list_ptr = psiinv_deviceptr_list.data();
    const auto* psi_list_ptr    = psi_deviceptr_list.data();

    {
      ScopedTimer local_timer(offload_timer);
      PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(2) \
          map(always, to: psiinv_list_ptr[:nw], psi_list_ptr[:nw]) \
          map(always, to: ratios_list_ptr[:nw], table_matrix_list_ptr[:nw]) \
	  map(to:first[:npairs], second[:npairs])")
      for (size_t iw = 0; iw < nw; iw++)
        for (size_t i = 0; i < npairs; ++i)
        {
          const int I = first[i];
          const int J = second[i];

          ValueType table_matrix_local = 0.0;
          for (size_t ind = 0; ind < num; ind++)
            table_matrix_local +=
                psiinv_list_ptr[iw][I * nb_cols_psiinv + ind] * psi_list_ptr[iw][J * nb_cols_psi + ind];
          table_matrix_list_ptr[iw][I * nb_cols_table_matrix + J] = table_matrix_local;
        }
    }
  }
  {
    ScopedTimer local(table2ratios_timer);
    const int max_ext_level = ndets_per_excitation_level_->size() - 1;

    // Compute workload changes drastically as the excitation level increases.
    // this may need different parallelization strategy.
    size_t det_offset  = 1;
    size_t data_offset = 1;

    auto update_offsets = [&](size_t ext_level) {
      det_offset += (*ndets_per_excitation_level_)[ext_level];
      data_offset += (*ndets_per_excitation_level_)[ext_level] * (3 * ext_level + 1);
    };

    mw_updateRatios_det0(det0_list, ratios_deviceptr_list);

    if (max_ext_level >= 1)
    {
      mw_updateRatios<1>(det_offset, data_offset, data, sign, table_matrix_deviceptr_list, nb_cols_table_matrix,
                         ratios_deviceptr_list);
      update_offsets(1);
    }

    if (max_ext_level >= 2)
    {
      mw_updateRatios<2>(det_offset, data_offset, data, sign, table_matrix_deviceptr_list, nb_cols_table_matrix,
                         ratios_deviceptr_list);
      update_offsets(2);
    }

    if (max_ext_level >= 3)
    {
      mw_updateRatios<3>(det_offset, data_offset, data, sign, table_matrix_deviceptr_list, nb_cols_table_matrix,
                         ratios_deviceptr_list);
      update_offsets(3);
    }

    if (max_ext_level >= 4)
    {
      mw_updateRatios<4>(det_offset, data_offset, data, sign, table_matrix_deviceptr_list, nb_cols_table_matrix,
                         ratios_deviceptr_list);
      update_offsets(4);
    }

    if (max_ext_level >= 5)
    {
      mw_updateRatios<5>(det_offset, data_offset, data, sign, table_matrix_deviceptr_list, nb_cols_table_matrix,
                         ratios_deviceptr_list);
      update_offsets(5);
    }

    if (max_ext_level >= 6)
    {
      for (size_t iw = 0; iw < nw; iw++)
        table_matrix_list[iw].get().updateFrom();
      for (size_t ext_level = 6; ext_level <= max_ext_level; ext_level++)
      {
        mw_updateRatios_generic(ext_level, det_offset, data_offset, det_calculator_, data, sign, table_matrix_list,
                                ratios_list);
        update_offsets(ext_level);
      }
      // FIXME need to transfer the part of det ratios ext_level >= 6 to the device.
    }
  }
}


void MultiDiracDeterminant::buildTableMatrix_calculateRatios(
    int ref,
    const OffloadMatrix<ValueType>& psiinv,
    const OffloadMatrix<ValueType>& psi,
    const OffloadVector<int>& data,
    const VectorSoaContainer<int, 2, OffloadPinnedAllocator<int>>& pairs,
    const OffloadVector<RealType>& sign,
    OffloadMatrix<ValueType>& table_matrix,
    OffloadVector<ValueType>& ratios)
{
  ScopedTimer local_timer(calculateRatios_timer);
  buildTableMatrix_calculateRatios_impl(ref, ValueType(1), ratios.data(), psiinv, psi, table_matrix, data, pairs, sign);
}

void MultiDiracDeterminant::mw_buildTableMatrix_calculateRatios(
    MultiDiracDetMultiWalkerResource& mw_res,
    int ref,
    const OffloadVector<ValueType>& det0_list,
    const RefVector<OffloadMatrix<ValueType>>& psiinv_list,
    const RefVector<OffloadMatrix<ValueType>>& psi_list,
    const OffloadVector<int>& data,
    const VectorSoaContainer<int, 2, OffloadPinnedAllocator<int>>& pairs,
    const OffloadVector<RealType>& sign,
    const RefVector<OffloadMatrix<ValueType>>& table_matrix_list,
    const RefVector<OffloadVector<ValueType>>& ratios_list)
{
  ScopedTimer local_timer(calculateRatios_timer);
  mw_buildTableMatrix_calculateRatios_impl(mw_res, ref, det0_list, psiinv_list, psi_list, data, pairs, sign,
                                           table_matrix_list, ratios_list);

  {
    ScopedTimer local_timer(transferD2H_timer);
    for (size_t iw = 0; iw < ratios_list.size(); iw++)
      ratios_list[iw].get().updateFrom();
  }
}

void MultiDiracDeterminant::buildTableMatrix_calculateGradRatios(
    int ref,
    const OffloadMatrix<ValueType>& psiinv,
    const OffloadMatrix<ValueType>& psi,
    const OffloadVector<int>& data,
    const VectorSoaContainer<int, 2, OffloadPinnedAllocator<int>>& pairs,
    const OffloadVector<RealType>& sign,
    const ValueType& det0_grad,
    OffloadMatrix<ValueType>& table_matrix,
    int dx,
    int iat,
    Matrix<GradType>& grads)
{
  ScopedTimer local_timer(calculateGradRatios_timer);
  buildTableMatrix_calculateRatios_impl(ref, det0_grad, WorkSpace.data(), psiinv, psi, table_matrix, data, pairs, sign);
  for (size_t count = 0; count < getNumDets(); ++count)
    grads(count, iat)[dx] = WorkSpace[count];
}

void MultiDiracDeterminant::mw_buildTableMatrix_calculateGradRatios(
    MultiDiracDetMultiWalkerResource& mw_res,
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
    const RefVector<OffloadMatrix<ValueType>>& table_matrix_list,
    UnpinnedOffloadMatrix<ValueType>& mw_grads)
{
  ScopedTimer local_timer(calculateGradRatios_timer);
  mw_buildTableMatrix_calculateRatios_impl(mw_res, ref, det0_grad_list, psiinv_list, psi_list, data, pairs, sign,
                                           table_matrix_list, WorkSpace_list);

  const size_t nw = WorkSpace_list.size();
  OffloadVector<ValueType*> WorkSpace_deviceptr_list(nw);
  for (size_t iw = 0; iw < nw; iw++)
    WorkSpace_deviceptr_list[iw] = WorkSpace_list[iw].get().device_data();

  auto* WorkSpace_list_ptr = WorkSpace_deviceptr_list.data();
  auto* mw_grads_ptr       = mw_grads.data();
  const size_t Grads_cols  = mw_grads.cols();

  PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(2)  map(from:mw_grads_ptr[:mw_grads.size()]) \
		                                                        map(always, to:WorkSpace_list_ptr[:nw])")
  for (size_t iw = 0; iw < nw; iw++)
    for (size_t count = 0; count < getNumDets; ++count)
      mw_grads_ptr[(3 * iw + dx) * Grads_cols + count] = WorkSpace_list_ptr[iw][count];
}

void MultiDiracDeterminant::buildTableMatrix_calculateRatiosValueMatrixOneParticle(
    int ref,
    const OffloadMatrix<ValueType>& psiinv,
    const OffloadMatrix<ValueType>& psi,
    const OffloadVector<int>& data,
    const VectorSoaContainer<int, 2, OffloadPinnedAllocator<int>>& pairs,
    const OffloadVector<RealType>& sign,
    OffloadMatrix<ValueType>& table_matrix,
    int iat,
    Matrix<ValueType>& ratios)
{
  const ValueType det0 = ratios(ref, iat);
  buildTableMatrix_calculateRatios_impl(ref, det0, WorkSpace.data(), psiinv, psi, table_matrix, data, pairs, sign);
  //splatt
  for (size_t count = 0; count < getNumDets(); ++count)
    ratios(count, iat) = WorkSpace[count];
#if 0
    ValueType det0 = ratios(ref,iat);
    int num=psi.extent(1);
    std::vector<std::pair<int,int> >::iterator it(pairs.begin()), last(pairs.end());
    while(it != last)
    {
      table_matrix((*it).first,(*it).second) = simd::dot(psiinv[(*it).first],psi[(*it).second],num);
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
      ratios(count,iat) = sign[count]*det0*CustomizedMatrixDet(n,table_matrix,it2+1);
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

  ScopedTimer local_timer(det_leader.evaluateDetsForPtclMove_timer);

  RefVector<OffloadVector<ValueType>> psiV_list, new_ratios_to_ref_list;
  RefVector<OffloadMatrix<ValueType>> TpsiM_list, psiM_list, table_matrix_list;
  RefVector<OffloadMatrix<ValueType>> psiMinv_temp_list, psiMinv_list;

  phi_list.reserve(nw);
  psiV_list.reserve(nw);
  psiMinv_list.reserve(nw);
  psiMinv_temp_list.reserve(nw);
  table_matrix_list.reserve(nw);
  TpsiM_list.reserve(nw);
  psiM_list.reserve(nw);
  new_ratios_to_ref_list.reserve(nw);

  for (size_t iw = 0; iw < nw; iw++)
  {
    MultiDiracDeterminant& det = (det_list[iw]);
    det.UpdateMode             = ORB_PBYP_PARTIAL;
    phi_list.push_back(*det.Phi);
    psiV_list.push_back(det.psiV);
    psiMinv_list.push_back(det.psiMinv);
    psiM_list.push_back(det.psiM);
    psiMinv_temp_list.push_back(det.psiMinv_temp);
    new_ratios_to_ref_list.push_back(det.new_ratios_to_ref_);
    table_matrix_list.push_back(det.table_matrix);
    TpsiM_list.push_back(det.TpsiM);
  }


  det_leader.UpdateMode  = ORB_PBYP_RATIO;
  const int WorkingIndex = iat - det_leader.FirstIndex;

  {
    ScopedTimer local_timer(det_leader.evalOrbValue_timer);
    for (size_t iw = 0; iw < nw; iw++)
    {
      MultiDiracDeterminant& det = (det_list[iw]);
      Vector<ValueType> psiV_list_host_view(psiV_list[iw].get().data(), psiV_list[iw].get().size());
      det.getPhi()->evaluateValue(P_list[iw], iat, psiV_list_host_view);
      ///Transfer of data from host to Device
      {
        ScopedTimer local_timer(det.transferH2D_timer);
        psiV_list[iw].get().updateTo();
      }
    }
  }

  size_t success          = 0;
  int dummy_handle        = 0;
  const auto psiMinv_rows = psiMinv_list[0].get().rows();
  const auto psiMinv_cols = psiMinv_list[0].get().cols();
  const auto TpsiM_cols   = TpsiM_list[0].get().cols();
  const auto psiM_cols    = psiM_list[0].get().cols();
  const auto TpsiM_rows   = TpsiM_list[0].get().rows();
  const auto NumPtcls     = det_leader.NumPtcls;
  const auto NumOrbitals  = det_leader.NumOrbitals;
  const auto& confgList   = *det_leader.ciConfigList;

  auto& mw_res             = det_leader.mw_res_handle_.getResource();
  auto* psiV_list_devptr   = mw_res.psiV_deviceptr_list.device_data();
  auto* psiV_temp_list_ptr = mw_res.psiV_temp_deviceptr_list.data();

  auto* psiMinv_list_devptr      = mw_res.psiMinv_deviceptr_list.device_data();
  auto* psiMinv_temp_list_devptr = mw_res.psiMinv_temp_deviceptr_list.device_data();

  auto* TpsiM_list_devptr = mw_res.TpsiM_deviceptr_list.device_data();
  auto* psiM_list_ptr     = mw_res.psiM_deviceptr_list.data();

  auto& curRatio_list = mw_res.curRatio_list;
  curRatio_list.resize(nw);
  auto* curRatio_list_ptr = curRatio_list.data();

  auto& inv_curRatio_list = mw_res.inv_curRatio_list;
  inv_curRatio_list.resize(nw);
  auto* inv_curRatio_list_ptr = inv_curRatio_list.data();

  auto* confgListOccup_ptr = det_leader.refdet_occup->data();

  {
    success = ompBLAS::copy_batched(dummy_handle, psiMinv_rows * psiMinv_cols, psiMinv_list_devptr, 1,
                                    psiMinv_temp_list_devptr, 1, nw);
    if (success != 0)
      throw std::runtime_error("In MultiDiracDeterminant ompBLAS::copy_batched_offset failed.");

    success = ompBLAS::copy_batched_offset(dummy_handle, det_leader.NumOrbitals, psiV_list_devptr, 0, 1,
                                           TpsiM_list_devptr, WorkingIndex, TpsiM_cols, nw);
    if (success != 0)
      throw std::runtime_error("In MultiDiracDeterminant ompBLAS::copy_batched_offset failed.");

    PRAGMA_OFFLOAD("omp target teams distribute map(always, from:curRatio_list_ptr[:nw]) \
                    is_device_ptr(psiV_list_devptr, psiMinv_temp_list_devptr)")
    for (size_t iw = 0; iw < nw; iw++)
    {
      ValueType c_ratio = 0.0;
      PRAGMA_OFFLOAD("omp parallel for reduction(+ : c_ratio)")
      for (size_t jc = 0; jc < NumPtcls; jc++)
      {
        const size_t J             = confgListOccup_ptr[jc];
        psiV_temp_list_ptr[iw][jc] = psiV_list_devptr[iw][J];
        size_t ic                  = jc * psiMinv_cols;
        c_ratio += (psiMinv_temp_list_devptr[iw] + WorkingIndex)[ic] * psiV_temp_list_ptr[iw][jc];
      }
      curRatio_list_ptr[iw]     = c_ratio;
      inv_curRatio_list_ptr[iw] = ValueType(1) / c_ratio;
    }

    det_leader.mw_InverseUpdateByColumn(det_leader.mw_res_handle_, WorkingIndex, inv_curRatio_list,
                                        mw_res.psiV_temp_deviceptr_list, mw_res.psiMinv_temp_deviceptr_list,
                                        psiMinv_rows);
    ///This is needed by acceptMove. Eventually acceptMove will need to become mw_acceptMove.
    {
      ScopedTimer local_timer(det_leader.transferD2H_timer);
      for (size_t iw = 0; iw < nw; iw++)
        psiMinv_temp_list[iw].get().updateFrom();
    }

    auto& det0_list = mw_res.cone_vec;
    det_leader.mw_buildTableMatrix_calculateRatios(det_leader.mw_res_handle_, det_leader.ReferenceDeterminant,
                                                   det0_list, psiMinv_temp_list, TpsiM_list, *det_leader.detData,
                                                   *det_leader.uniquePairs, *det_leader.DetSigns, table_matrix_list,
                                                   new_ratios_to_ref_list);

    // restore the modified column of TpsiM.
    PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(2) is_device_ptr(TpsiM_list_devptr) \
		                                    map(always, to:psiM_list_ptr[:nw])")
    for (size_t iw = 0; iw < nw; iw++)
      for (size_t i = 0; i < NumOrbitals; i++)
        TpsiM_list_devptr[iw][i * TpsiM_cols + WorkingIndex] = psiM_list_ptr[iw][i + psiM_cols * WorkingIndex];
  }

  for (size_t iw = 0; iw < nw; iw++)
  {
    MultiDiracDeterminant& det = (det_list[iw]);
    det.curRatio               = curRatio_list_ptr[iw];
  }
}

void MultiDiracDeterminant::evaluateDetsForPtclMove(const ParticleSet& P, int iat, int refPtcl)
{
  ScopedTimer local_timer(evaluateDetsForPtclMove_timer);

  UpdateMode = ORB_PBYP_RATIO;
  {
    ScopedTimer orb_timer(evalOrbValue_timer);
    Vector<ValueType> psiV_host_view(psiV.data(), psiV.size());
    Phi->evaluateValue(P, iat, psiV_host_view);
  }
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
  {
    ScopedTimer inverse(updateInverse_timer);
    psiMinv_temp = psiMinv;
    for (size_t i = 0; i < NumPtcls; i++)
      psiV_temp[i] = psiV[*(it++)];
    auto ratio_old_ref_det = DetRatioByColumn(psiMinv_temp, psiV_temp, WorkingIndex);
    curRatio               = ratio_old_ref_det;
    InverseUpdateByColumn(psiMinv_temp, psiV_temp, workV1, workV2, WorkingIndex, ratio_old_ref_det);
    for (size_t i = 0; i < NumOrbitals; i++)
      TpsiM(i, WorkingIndex) = psiV[i];
  }
  buildTableMatrix_calculateRatios(ReferenceDeterminant, psiMinv_temp, TpsiM, *detData, *uniquePairs, *DetSigns,
                                   table_matrix, new_ratios_to_ref_);
  // check comment above
  for (size_t i = 0; i < NumOrbitals; i++)
    TpsiM(i, WorkingIndex) = psiM(WorkingIndex, i);
}

void MultiDiracDeterminant::evaluateDetsAndGradsForPtclMove(const ParticleSet& P, int iat)
{
  ScopedTimer local_timer(evaluateDetsAndGradsForPtclMove_timer);
  UpdateMode = ORB_PBYP_PARTIAL;
  {
    ScopedTimer local_timer(evalOrbVGL_timer);
    // Using Host Views for Phi-evaluateVGL since not ported to GPU
    Vector<ValueType> psiV_host_view(psiV.data(), psiV.size());
    Vector<GradType> dpsiV_host_view(dpsiV.data(), dpsiV.size());
    Vector<ValueType> d2psiV_host_view(d2psiV.data(), d2psiV.size());
    Phi->evaluateVGL(P, iat, psiV_host_view, dpsiV_host_view, d2psiV_host_view);
  }
  const int WorkingIndex = iat - FirstIndex;
  const auto& confgList  = *ciConfigList;
  assert(WorkingIndex >= 0 && WorkingIndex < LastIndex - FirstIndex);

  GradType ratioGradRef;
  {
    ScopedTimer inverse(updateInverse_timer);
    //mmorales: check comment above
    psiMinv_temp = psiMinv;
    //std::vector<int>::iterator it(confgList[ReferenceDeterminant].occup.begin());
    auto it(confgList[ReferenceDeterminant].occup.begin());
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
  }
  buildTableMatrix_calculateRatios(ReferenceDeterminant, psiMinv_temp, TpsiM, *detData, *uniquePairs, *DetSigns,
                                   table_matrix, new_ratios_to_ref_);
  for (size_t idim = 0; idim < OHMMS_DIM; idim++)
  {
    {
      ScopedTimer inverse(updateInverse_timer);
      //dpsiMinv = psiMinv_temp;
      dpsiMinv = psiMinv;
      auto it(confgList[ReferenceDeterminant].occup.begin());
      for (size_t i = 0; i < NumPtcls; i++)
        psiV_temp[i] = dpsiV[*(it++)][idim];
      InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, WorkingIndex, ratioGradRef[idim]);
      for (size_t i = 0; i < NumOrbitals; i++)
        TpsiM(i, WorkingIndex) = dpsiV[i][idim];
    }
    buildTableMatrix_calculateGradRatios(ReferenceDeterminant, dpsiMinv, TpsiM, *detData, *uniquePairs, *DetSigns,
                                         ratioGradRef[idim] / curRatio, table_matrix, idim, WorkingIndex, new_grads);
  }
  // check comment above
  for (int i = 0; i < NumOrbitals; i++)
    TpsiM(i, WorkingIndex) = psiM(WorkingIndex, i);
}

void MultiDiracDeterminant::evaluateDetsAndGradsForPtclMoveWithSpin(const ParticleSet& P, int iat)
{
  assert(P.isSpinor() == is_spinor_);
  ScopedTimer local_timer(evaluateDetsAndGradsForPtclMove_timer);
  UpdateMode = ORB_PBYP_PARTIAL;
  {
    ScopedTimer orb_timer(evalOrbVGL_timer);
    // Creating Host view to call Phi->evaluateVGL
    Vector<ValueType> psiV_host_view(psiV.data(), psiV.size());
    Vector<GradType> dpsiV_host_view(dpsiV.data(), dpsiV.size());
    Vector<ValueType> d2psiV_host_view(d2psiV.data(), d2psiV.size());
    Phi->evaluateVGL_spin(P, iat, psiV_host_view, dpsiV_host_view, d2psiV_host_view, dspin_psiV);
  }
  const int WorkingIndex = iat - FirstIndex;
  assert(WorkingIndex >= 0 && WorkingIndex < LastIndex - FirstIndex);
  const auto& confgList = *ciConfigList;
  GradType ratioGradRef;
  ValueType ratioSpinGradRef = 0.0;
  {
    ScopedTimer inverse(updateInverse_timer);
    //mmorales: check comment above
    psiMinv_temp = psiMinv;
    //std::vector<int>::iterator it(confgList[ReferenceDeterminant].occup.begin());
    auto it(confgList[ReferenceDeterminant].occup.begin());
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
  }
  buildTableMatrix_calculateRatios(ReferenceDeterminant, psiMinv_temp, TpsiM, *detData, *uniquePairs, *DetSigns,
                                   table_matrix, new_ratios_to_ref_);
  for (size_t idim = 0; idim < OHMMS_DIM; idim++)
  {
    {
      ScopedTimer inverse(updateInverse_timer);
      //dpsiMinv = psiMinv_temp;
      dpsiMinv = psiMinv;
      auto it(confgList[ReferenceDeterminant].occup.begin());
      for (size_t i = 0; i < NumPtcls; i++)
        psiV_temp[i] = dpsiV[*(it++)][idim];
      InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, WorkingIndex, ratioGradRef[idim]);
      for (size_t i = 0; i < NumOrbitals; i++)
        TpsiM(i, WorkingIndex) = dpsiV[i][idim];
    }
    buildTableMatrix_calculateGradRatios(ReferenceDeterminant, dpsiMinv, TpsiM, *detData, *uniquePairs, *DetSigns,
                                         ratioGradRef[idim] / curRatio, table_matrix, idim, WorkingIndex, new_grads);
  }

  //Now compute the spin gradient, same procedure as normal gradient components above
  {
    ScopedTimer inverse(updateInverse_timer);
    dpsiMinv = psiMinv;
    auto it(confgList[ReferenceDeterminant].occup.begin());
    for (size_t i = 0; i < NumPtcls; i++)
      psiV_temp[i] = dspin_psiV[*(it++)];
    InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, WorkingIndex, ratioSpinGradRef);
    for (size_t i = 0; i < NumOrbitals; i++)
      TpsiM(i, WorkingIndex) = dspin_psiV[i];
  }
  buildTableMatrix_calculateRatiosValueMatrixOneParticle(ReferenceDeterminant, dpsiMinv, TpsiM, *detData, *uniquePairs,
                                                         *DetSigns, table_matrix, WorkingIndex, new_spingrads);

  // check comment above
  for (int i = 0; i < NumOrbitals; i++)
    TpsiM(i, WorkingIndex) = psiM(WorkingIndex, i);
}

void MultiDiracDeterminant::mw_evaluateDetsAndGradsForPtclMove(
    const RefVectorWithLeader<MultiDiracDeterminant>& det_list,
    const RefVectorWithLeader<ParticleSet>& P_list,
    int iat,
    UnpinnedOffloadMatrix<ValueType>& mw_grads)
{
  const int nw                      = det_list.size();
  MultiDiracDeterminant& det_leader = det_list.getLeader();
  RefVectorWithLeader<SPOSet> phi_list(*det_leader.getPhi());

  ScopedTimer local_timer(det_leader.evaluateDetsAndGradsForPtclMove_timer);
  int success      = 0;
  int dummy_handle = 0;
  const size_t NumOrbitals(det_leader.NumOrbitals);
  const size_t NumPtcls(det_leader.NumPtcls);

  RefVector<OffloadVector<ValueType>> psiV_list, psiV_temp_list, new_ratios_to_ref_list, WorkSpace_list;
  RefVector<OffloadVector<ValueType>> d2psiV_list;

  RefVector<OffloadMatrix<ValueType>> psiMinv_temp_list, psiMinv_list, dpsiMinv_list;
  RefVector<OffloadMatrix<ValueType>> table_matrix_list, psiM_list, TpsiM_list;

  RefVector<OffloadVector<GradType>> dpsiV_list;

  phi_list.reserve(nw);
  psiV_list.reserve(nw);
  dpsiV_list.reserve(nw);
  d2psiV_list.reserve(nw);
  psiV_temp_list.reserve(nw);
  psiMinv_temp_list.reserve(nw);
  psiMinv_list.reserve(nw);
  psiM_list.reserve(nw);

  TpsiM_list.reserve(nw);
  new_ratios_to_ref_list.reserve(nw);
  table_matrix_list.reserve(nw);
  dpsiMinv_list.reserve(nw);
  WorkSpace_list.reserve(nw);

  auto& mw_res            = det_leader.mw_res_handle_.getResource();
  auto& ratioGradRef_list = mw_res.ratioGradRef_list;
  auto& det0_grad_list    = mw_res.det0_grad_list;
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
    psiV_temp_list.push_back(det.psiV_temp);
    psiMinv_list.push_back(det.psiMinv);
    psiM_list.push_back(det.psiM);
    psiMinv_temp_list.push_back(det.psiMinv_temp);
    new_ratios_to_ref_list.push_back(det.new_ratios_to_ref_);
    TpsiM_list.push_back(det.TpsiM);
    table_matrix_list.push_back(det.table_matrix);
    dpsiMinv_list.push_back(det.dpsiMinv);
    WorkSpace_list.push_back(det.WorkSpace);
  }

  {
    ScopedTimer local_timer(det_leader.evalOrbVGL_timer);
    for (size_t iw = 0; iw < nw; iw++)
    {
      MultiDiracDeterminant& det = (det_list[iw]);
      Vector<ValueType> psiV_list_host_view(psiV_list[iw].get().data(), psiV_list[iw].get().size());
      Vector<GradType> dpsiV_list_host_view(dpsiV_list[iw].get().data(), dpsiV_list[iw].get().size());
      Vector<ValueType> d2psiV_list_host_view(d2psiV_list[iw].get().data(), d2psiV_list[iw].get().size());
      det.Phi->evaluateVGL(P_list[iw], iat, psiV_list_host_view, dpsiV_list_host_view, d2psiV_list_host_view);

      {
        ScopedTimer local_timer(det_leader.transferH2D_timer);
        psiV_list[iw].get().updateTo();
        dpsiV_list[iw].get().updateTo();
      }
    }
  }

  const auto psiMinv_rows   = psiMinv_list[0].get().rows();
  const auto psiMinv_cols   = psiMinv_list[0].get().cols();
  const auto TpsiM_num_cols = TpsiM_list[0].get().cols();
  const auto psiM_num_cols  = psiM_list[0].get().cols();
  const auto& confgList     = *det_leader.ciConfigList;

  auto* psiV_list_devptr         = mw_res.psiV_deviceptr_list.device_data();
  auto* psiV_temp_list_ptr       = mw_res.psiV_temp_deviceptr_list.data();
  auto* TpsiM_list_devptr        = mw_res.TpsiM_deviceptr_list.device_data();
  auto* psiM_list_ptr            = mw_res.psiM_deviceptr_list.data();
  auto* psiMinv_list_devptr      = mw_res.psiMinv_deviceptr_list.device_data();
  auto* dpsiMinv_list_devptr     = mw_res.dpsiMinv_deviceptr_list.device_data();
  auto* psiMinv_temp_list_devptr = mw_res.psiMinv_temp_deviceptr_list.device_data();

  auto& curRatio_list = mw_res.curRatio_list;
  curRatio_list.resize(nw);
  auto* curRatio_list_ptr = curRatio_list.data();

  auto& inv_curRatio_list = mw_res.inv_curRatio_list;
  inv_curRatio_list.resize(nw);
  auto* inv_curRatio_list_ptr = inv_curRatio_list.data();

  auto* det0_grad_list_ptr = det0_grad_list.data();
  auto* confgListOccup_ptr = det_leader.refdet_occup->data();

  auto* dpsiV_list_ptr        = mw_res.dpsiV_deviceptr_list.data();
  auto* ratioGradRef_list_ptr = ratioGradRef_list.data();

  {
    success = ompBLAS::copy_batched(dummy_handle, psiMinv_rows * psiMinv_cols, psiMinv_list_devptr, 1,
                                    psiMinv_temp_list_devptr, 1, nw);
    if (success != 0)
      throw std::runtime_error("In MultiDiracDeterminant ompBLAS::copy_batched_offset failed.");


    PRAGMA_OFFLOAD("omp target teams distribute is_device_ptr(psiV_list_devptr, psiMinv_temp_list_devptr) \
		                    map(always, from:curRatio_list_ptr[:nw])")
    for (size_t iw = 0; iw < nw; iw++)
    {
      GradType ratioGradRef_local(0);
      PRAGMA_OFFLOAD("omp parallel for reduction(+ : ratioGradRef_local)")
      for (size_t i = 0; i < NumPtcls; i++)
      {
        const size_t J            = confgListOccup_ptr[i];
        psiV_temp_list_ptr[iw][i] = psiV_list_devptr[iw][J];
        ratioGradRef_local += psiMinv_temp_list_devptr[iw][i * psiMinv_cols + WorkingIndex] * dpsiV_list_ptr[iw][J];
      }
      ratioGradRef_list_ptr[iw] = ratioGradRef_local;

      ValueType c_ratio = 0.0;
      PRAGMA_OFFLOAD("omp parallel for reduction(+ : c_ratio)")
      for (size_t jc = 0; jc < psiMinv_cols; jc += 1)
      {
        const size_t ic = jc * psiMinv_cols;
        c_ratio += (psiMinv_temp_list_devptr[iw] + WorkingIndex)[ic] * psiV_temp_list_ptr[iw][jc];
      }
      curRatio_list_ptr[iw]     = c_ratio;
      inv_curRatio_list_ptr[iw] = ValueType(1) / c_ratio;
    }

    success = ompBLAS::copy_batched_offset(dummy_handle, det_leader.NumOrbitals, psiV_list_devptr, 0, 1,
                                           TpsiM_list_devptr, WorkingIndex, TpsiM_num_cols, nw);
    if (success != 0)
      throw std::runtime_error("In MultiDiracDeterminant ompBLAS::copy_batched_offset failed.");


    det_leader.mw_InverseUpdateByColumn(det_leader.mw_res_handle_, WorkingIndex, inv_curRatio_list,
                                        mw_res.psiV_temp_deviceptr_list, mw_res.psiMinv_temp_deviceptr_list,
                                        psiMinv_rows);

    ///This is needed by Host in acceptMove. Eventually acceptMove will need to become mw_acceptMove.
    {
      ScopedTimer local_timer(det_leader.transferD2H_timer);
      for (size_t iw = 0; iw < nw; iw++)
        psiMinv_temp_list[iw].get().updateFrom();
    }

    auto& det0_list = mw_res.cone_vec;
    det_leader.mw_buildTableMatrix_calculateRatios(det_leader.mw_res_handle_, det_leader.ReferenceDeterminant,
                                                   det0_list, psiMinv_temp_list, TpsiM_list, *det_leader.detData,
                                                   *det_leader.uniquePairs, *det_leader.DetSigns, table_matrix_list,
                                                   new_ratios_to_ref_list);

    for (size_t idim = 0; idim < OHMMS_DIM; idim++)
    {
      success = ompBLAS::copy_batched(dummy_handle, psiMinv_rows * psiMinv_cols, psiMinv_list_devptr, 1,
                                      dpsiMinv_list_devptr, 1, nw);
      if (success != 0)
        throw std::runtime_error("In MultiDiracDeterminant ompBLAS::copy_batched_offset failed.");

      PRAGMA_OFFLOAD("omp target teams distribute map(to: ratioGradRef_list_ptr[:nw])")
      for (size_t iw = 0; iw < nw; iw++)
      {
        inv_curRatio_list_ptr[iw] = ValueType(1) / ratioGradRef_list_ptr[iw][idim];

        for (size_t i = 0; i < NumPtcls; i++)
        {
          const size_t J            = confgListOccup_ptr[i];
          psiV_temp_list_ptr[iw][i] = dpsiV_list_ptr[iw][J][idim];
        }
      }

      det_leader.mw_InverseUpdateByColumn(det_leader.mw_res_handle_, WorkingIndex, inv_curRatio_list,
                                          mw_res.psiV_temp_deviceptr_list, mw_res.dpsiMinv_deviceptr_list,
                                          psiMinv_rows);

      PRAGMA_OFFLOAD("omp target teams distribute map(to:dpsiV_list_ptr[:nw], curRatio_list_ptr[:nw]) \
		                       map(always,from:det0_grad_list_ptr[:nw]) \
		                       is_device_ptr(TpsiM_list_devptr)")
      for (size_t iw = 0; iw < nw; iw++)
      {
        det0_grad_list_ptr[iw] = ratioGradRef_list_ptr[iw][idim] / curRatio_list_ptr[iw];
        for (size_t i = 0; i < NumOrbitals; i++)
          TpsiM_list_devptr[iw][i * TpsiM_num_cols + WorkingIndex] = dpsiV_list_ptr[iw][i][idim];
      }

      det_leader.mw_buildTableMatrix_calculateGradRatios(det_leader.mw_res_handle_, det_leader.ReferenceDeterminant,
                                                         WorkingIndex, idim, det_leader.getNumDets(), det0_grad_list,
                                                         dpsiMinv_list, TpsiM_list, *det_leader.detData,
                                                         *det_leader.uniquePairs, *det_leader.DetSigns, WorkSpace_list,
                                                         table_matrix_list, mw_grads);
    }

    // restore the modified column of TpsiM.
    PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(2) is_device_ptr(TpsiM_list_devptr) \
		                                    map(always, to:psiM_list_ptr[:nw])")
    for (size_t iw = 0; iw < nw; iw++)
      for (size_t i = 0; i < NumOrbitals; i++)
        TpsiM_list_devptr[iw][i * TpsiM_num_cols + WorkingIndex] = psiM_list_ptr[iw][i + psiM_num_cols * WorkingIndex];
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
    buildTableMatrix_calculateGradRatios(ReferenceDeterminant, dpsiMinv, TpsiM, *detData, *uniquePairs, *DetSigns,
                                         ratioG, table_matrix, idim, WorkingIndex, grads);
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
    buildTableMatrix_calculateGradRatios(ReferenceDeterminant, dpsiMinv, TpsiM, *detData, *uniquePairs, *DetSigns,
                                         ratioG, table_matrix, idim, WorkingIndex, grads);
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
  buildTableMatrix_calculateRatiosValueMatrixOneParticle(ReferenceDeterminant, dpsiMinv, TpsiM, *detData, *uniquePairs,
                                                         *DetSigns, table_matrix, WorkingIndex, spingrads);

  // check comment above
  for (size_t i = 0; i < NumOrbitals; i++)
    TpsiM(i, WorkingIndex) = psiM(WorkingIndex, i);
}

void MultiDiracDeterminant::mw_evaluateGrads(const RefVectorWithLeader<MultiDiracDeterminant>& det_list,
                                             const RefVectorWithLeader<ParticleSet>& P_list,
                                             int iat,
                                             UnpinnedOffloadMatrix<ValueType>& mw_grads)
{
  const int nw                      = det_list.size();
  MultiDiracDeterminant& det_leader = det_list.getLeader();
  const int WorkingIndex            = iat - det_leader.FirstIndex;
  int success                       = 0;
  int dummy_handle                  = 0;
  const size_t NumOrbitals(det_leader.NumOrbitals);
  const size_t NumPtcls(det_leader.NumPtcls);

  ScopedTimer local_timer(det_leader.evaluateGrads_timer);

  RefVector<OffloadMatrix<ValueType>> dpsiMinv_list, psiMinv_list;
  RefVector<OffloadMatrix<GradType>> dpsiM_list;
  RefVector<OffloadVector<ValueType>> psiV_temp_list, WorkSpace_list;
  RefVector<OffloadMatrix<ValueType>> table_matrix_list, TpsiM_list, psiM_list;
  //RefVector<OffloadMatrix<GradType>> grads_list;

  psiMinv_list.reserve(nw);
  dpsiMinv_list.reserve(nw);
  dpsiM_list.reserve(nw);
  psiV_temp_list.reserve(nw);
  //grads_list.reserve(nw);
  table_matrix_list.reserve(nw);
  TpsiM_list.reserve(nw);
  psiM_list.reserve(nw);
  WorkSpace_list.reserve(nw);

  for (size_t iw = 0; iw < nw; iw++)
  {
    MultiDiracDeterminant& det = (det_list[iw]);
    psiMinv_list.push_back(det.psiMinv);
    dpsiMinv_list.push_back(det.dpsiMinv);
    psiV_temp_list.push_back(det.psiV_temp);
    dpsiM_list.push_back(det.dpsiM);
    //grads_list.push_back(det.grads);
    TpsiM_list.push_back(det.TpsiM);
    psiM_list.push_back(det.psiM);
    table_matrix_list.push_back(det.table_matrix);
    WorkSpace_list.push_back(det.WorkSpace);
  }

  const auto psiMinv_rows = psiMinv_list[0].get().rows();
  const auto psiMinv_cols = psiMinv_list[0].get().cols();
  const auto TpsiM_cols   = TpsiM_list[0].get().cols();
  const auto psiM_cols    = psiM_list[0].get().cols();
  const auto dpsiM_cols   = dpsiM_list[0].get().cols();
  const auto dpsiM_rows   = dpsiM_list[0].get().rows();
  const auto& confgList   = *det_leader.ciConfigList;

  auto& mw_res      = det_leader.mw_res_handle_.getResource();
  auto& ratioG_list = mw_res.curRatio_list;
  ratioG_list.resize(nw);
  auto* ratioG_list_ptr = ratioG_list.data();

  auto& inv_ratioG_list = mw_res.inv_curRatio_list;
  inv_ratioG_list.resize(nw);
  auto* inv_ratioG_list_ptr = inv_ratioG_list.data();

  auto* psiMinv_list_devptr  = mw_res.psiMinv_deviceptr_list.device_data();
  auto* dpsiMinv_list_devptr = mw_res.dpsiMinv_deviceptr_list.device_data();
  auto* TpsiM_list_devptr    = mw_res.TpsiM_deviceptr_list.data();
  auto* psiM_list_ptr        = mw_res.psiM_deviceptr_list.data();
  auto* psiV_temp_list_ptr   = mw_res.psiV_temp_deviceptr_list.data();
  auto* dpsiM_list_ptr       = mw_res.dpsiM_deviceptr_list.data();
  auto* confgListOccup_ptr   = det_leader.refdet_occup->data();

  {
    for (size_t idim = 0; idim < OHMMS_DIM; idim++)
    {
      success = ompBLAS::copy_batched(dummy_handle, psiMinv_rows * psiMinv_cols, psiMinv_list_devptr, 1,
                                      dpsiMinv_list_devptr, 1, nw);
      if (success != 0)
        throw std::runtime_error("In MultiDiracDeterminant ompBLAS::copy_batched_offset failed.");

      PRAGMA_OFFLOAD("omp target teams distribute  map(always, to: psiV_temp_list_ptr[:nw]) \
		                                  map(always, to: dpsiM_list_ptr[:nw])")
      for (size_t iw = 0; iw < nw; iw++)
        for (size_t i = 0; i < NumPtcls; i++)
        {
          size_t J                  = confgListOccup_ptr[i];
          psiV_temp_list_ptr[iw][i] = dpsiM_list_ptr[iw][WorkingIndex * dpsiM_cols + J][idim];
        }

      PRAGMA_OFFLOAD("omp target teams distribute is_device_ptr(psiMinv_list_devptr) \
		                                  map(always, from: ratioG_list_ptr[:nw])")
      for (size_t iw = 0; iw < nw; iw++)
      {
        ValueType ratioG_local(0);
        PRAGMA_OFFLOAD("omp parallel for reduction(+ : ratioG_local)")
        for (size_t i = 0; i < NumPtcls; i++)
        {
          size_t J = confgListOccup_ptr[i];
          ratioG_local += psiMinv_list_devptr[iw][i * psiMinv_cols + WorkingIndex] *
              dpsiM_list_ptr[iw][WorkingIndex * dpsiM_cols + J][idim];
        }
        ratioG_list_ptr[iw]     = ratioG_local;
        inv_ratioG_list_ptr[iw] = ValueType(1) / ratioG_local;
      }

      det_leader.mw_InverseUpdateByColumn(det_leader.mw_res_handle_, WorkingIndex, inv_ratioG_list,
                                          mw_res.psiV_temp_deviceptr_list, mw_res.dpsiMinv_deviceptr_list,
                                          psiMinv_rows);

      PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(2) map(to:dpsiM_list_ptr[:nw]) \
		                                              map(always, to: TpsiM_list_devptr[:nw])")
      for (size_t iw = 0; iw < nw; iw++)
        for (size_t i = 0; i < NumOrbitals; i++)
          TpsiM_list_devptr[iw][i * TpsiM_cols + WorkingIndex] =
              dpsiM_list_ptr[iw][dpsiM_cols * WorkingIndex + i][idim];


      det_leader.mw_buildTableMatrix_calculateGradRatios(det_leader.mw_res_handle_, det_leader.ReferenceDeterminant,
                                                         WorkingIndex, idim, det_leader.getNumDets(), ratioG_list,
                                                         dpsiMinv_list, TpsiM_list, *det_leader.detData,
                                                         *det_leader.uniquePairs, *det_leader.DetSigns, WorkSpace_list,
                                                         table_matrix_list, mw_grads);
    }

    // restore the modified column of TpsiM.
    PRAGMA_OFFLOAD("omp target teams distribute parallel for map(from:TpsiM_list_devptr[:nw]) \
                                                       map(always,to:psiM_list_ptr[:nw])")
    for (size_t iw = 0; iw < nw; iw++)
      for (size_t i = 0; i < NumOrbitals; i++)
        TpsiM_list_devptr[iw][i * TpsiM_cols + WorkingIndex] = psiM_list_ptr[iw][i + psiM_cols * WorkingIndex];
  }
}

void MultiDiracDeterminant::mw_updateRatios_generic(int ext_level,
                                                    const size_t det_offset,
                                                    const size_t data_offset,
                                                    SmallMatrixDetCalculator<ValueType>& det_calculator,
                                                    const OffloadVector<int>& data,
                                                    const OffloadVector<RealType>& sign,
                                                    const RefVector<OffloadMatrix<ValueType>>& table_matrix_list,
                                                    const RefVector<OffloadVector<ValueType>>& ratios_list) const
{
  const size_t nw = ratios_list.size();
  const int* it2  = data.data() + data_offset;
  for (size_t iw = 0; iw < nw; iw++)
    for (size_t count = 0; count < (*ndets_per_excitation_level_)[ext_level]; ++count)
    {
      size_t det_id                 = det_offset + count;
      ratios_list[iw].get()[det_id] = sign[det_id] * ratios_list[iw].get()[0] *
          det_calculator.evaluate(table_matrix_list[iw].get(), it2 + 1 + count * (3 * ext_level + 1), ext_level);
    }
}

void MultiDiracDeterminant::mw_updateRatios_det0(const OffloadVector<ValueType>& det0_list,
                                                 const OffloadVector<ValueType*>& ratios_deviceptr_list) const
{
  ScopedTimer local_timer(updateRatios_timer);

  const size_t nw           = ratios_deviceptr_list.size();
  auto* ratios_list_ptr     = ratios_deviceptr_list.data();
  const auto* det0_list_ptr = det0_list.data();

  PRAGMA_OFFLOAD("omp target teams distribute parallel for map(always, to: det0_list_ptr[:nw])")
  for (size_t iw = 0; iw < nw; iw++)
    ratios_list_ptr[iw][0] = det0_list_ptr[iw];
}

template<unsigned EXT_LEVEL>
void MultiDiracDeterminant::mw_updateRatios(const size_t det_offset,
                                            const size_t data_offset,
                                            const OffloadVector<int>& data,
                                            const OffloadVector<RealType>& sign,
                                            const OffloadVector<ValueType*>& table_matrix_deviceptr_list,
                                            const size_t num_table_matrix_cols,
                                            const OffloadVector<ValueType*>& ratios_deviceptr_list) const
{
  const size_t nw        = ratios_deviceptr_list.size();
  const size_t size_sign = sign.size();
  const size_t ndet_ext  = (*ndets_per_excitation_level_)[EXT_LEVEL];

  ScopedTimer local_timer(updateRatios_timer);

  auto* ratios_list_ptr             = ratios_deviceptr_list.data();
  const auto* sign_ptr              = sign.data();
  const int* data_ptr               = data.data();
  const auto* table_matrix_list_ptr = table_matrix_deviceptr_list.data();

  PRAGMA_OFFLOAD("omp target teams distribute parallel for collapse(2)")
  for (size_t iw = 0; iw < nw; iw++)
    for (size_t count = 0; count < ndet_ext; ++count)
    {
      size_t det_id               = det_offset + count;
      ratios_list_ptr[iw][det_id] = sign_ptr[det_id] * ratios_list_ptr[iw][0] *
          CustomizedMatrixDet<EXT_LEVEL>::evaluate(table_matrix_list_ptr[iw],
                                                   (data_ptr + data_offset) + 1 + count * (3 * EXT_LEVEL + 1),
                                                   num_table_matrix_cols);
    }
}

void MultiDiracDeterminant::mw_InverseUpdateByColumn(MultiDiracDetMultiWalkerResource& mw_res,
                                                     const int working_index,
                                                     const OffloadVector<ValueType>& inv_curRatio_list,
                                                     const OffloadVector<ValueType*>& psiV_deviceptr_list,
                                                     const OffloadVector<ValueType*>& psiMinv_deviceptr_list,
                                                     const size_t psiMinv_rows) const
{
  ScopedTimer local_timer(updateInverse_timer);

  constexpr ValueType cone(1);

  const size_t nw = inv_curRatio_list.size();

  ValueType* czero_ptr      = mw_res.czero_vec.device_data();
  ValueType* cminus_one_ptr = mw_res.cminus_one_vec.device_data();

  int success      = 0;
  int dummy_handle = 0;

  auto* psiV_list_devptr         = psiV_deviceptr_list.device_data();
  auto* psiMinv_list_devptr      = psiMinv_deviceptr_list.device_data();
  auto* workV1_list_ptr          = mw_res.workV1_deviceptr_list.device_data();
  auto* workV2_list_ptr          = mw_res.workV2_deviceptr_list.device_data();
  auto* inv_curRatio_list_devptr = inv_curRatio_list.device_data();

  success =
      ompBLAS::gemv_batched(dummy_handle, 'N', psiMinv_rows, psiMinv_rows, inv_curRatio_list_devptr,
                            psiMinv_list_devptr, psiMinv_rows, psiV_list_devptr, 1, czero_ptr, workV1_list_ptr, 1, nw);
  if (success != 0)
    throw std::runtime_error("In MultiDiracDeterminant ompBLAS::gemv_batched failed.");

  PRAGMA_OFFLOAD("omp target teams distribute parallel for is_device_ptr(workV1_list_ptr, inv_curRatio_list_devptr)")
  for (size_t iw = 0; iw < nw; iw++)
    workV1_list_ptr[iw][working_index] = cone - inv_curRatio_list_devptr[iw];
  if (success != 0)
    throw std::runtime_error("In MultiDiracDeterminant ompBLAS::copy_batched_offset failed.");

  success = ompBLAS::copy_batched_offset(dummy_handle, psiMinv_rows, psiMinv_list_devptr, working_index, psiMinv_rows,
                                         workV2_list_ptr, 0, 1, nw);
  if (success != 0)
    throw std::runtime_error("In MultiDiracDeterminant ompBLAS::copy_batched_offset failed.");

  success = ompBLAS::ger_batched(dummy_handle, psiMinv_rows, psiMinv_rows, cminus_one_ptr, workV1_list_ptr, 1,
                                 workV2_list_ptr, 1, psiMinv_list_devptr, psiMinv_rows, nw);
  if (success != 0)
    throw std::runtime_error("In MultiDiracDeterminant ompBLAS::ger_batched failed.");
}

} // namespace qmcplusplus
