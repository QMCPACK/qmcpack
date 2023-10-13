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


#include "MultiDiracDeterminant.h"
#include "QMCWaveFunctions/Fermion/ci_configuration2.h"
#ifndef QMC_COMPLEX
#include "QMCWaveFunctions/RotatedSPOs.h"
#endif
#include "Message/Communicate.h"
#include "Numerics/DeterminantOperators.h"
#include "CPU/BLAS.hpp"
#include "Numerics/MatrixOperators.h"
#include <algorithm>
#include <vector>

// mmorales:
// NOTE NOTE NOTE:
// right now the code assumes that all the orbitals in the active space are used,
// this means that there can be problems if some of the orbitals are not used

namespace qmcplusplus
{
void MultiDiracDeterminant::createDetData(const int ref_det_id,
                                          const std::vector<ci_configuration2>& configlist_unsorted,
                                          const std::vector<size_t>& C2nodes_unsorted,
                                          std::vector<size_t>& C2nodes_sorted)
{
  auto& ref                        = configlist_unsorted[ref_det_id];
  auto& configlist_sorted          = *ciConfigList;
  auto& data                       = *detData;
  auto& pairs                      = *uniquePairs;
  auto& sign                       = *DetSigns;
  auto& ndets_per_excitation_level = *ndets_per_excitation_level_;

  const size_t nci = configlist_unsorted.size();
  std::vector<std::pair<int, int>> pairs_local;

  size_t nex_max = 0;
  std::vector<size_t> pos(NumPtcls);
  std::vector<size_t> ocp(NumPtcls);
  std::vector<size_t> uno(NumPtcls);
  // map key is exc. lvl
  std::map<int, std::vector<int>> dataMap;
  std::map<int, std::vector<int>> sortMap;
  std::vector<RealType> tmp_sign(nci, 0);
  for (size_t i = 0; i < nci; i++)
  {
    size_t nex;
    tmp_sign[i] = ref.calculateExcitations(configlist_unsorted[i], nex, pos, ocp, uno);
    nex_max     = std::max(nex, nex_max);
    dataMap[nex].push_back(nex);
    sortMap[nex].push_back(i);
    for (int k = 0; k < nex; k++)
      dataMap[nex].push_back(pos[k]);
    for (int k = 0; k < nex; k++)
      dataMap[nex].push_back(uno[k]);
    for (int k = 0; k < nex; k++)
      dataMap[nex].push_back(ocp[k]);
    // determine unique pairs, to avoid redundant calculation of matrix elements
    // if storing the entire MOxMO matrix is too much, then make an array and a mapping to it.
    // is there an easier way??
    for (int k1 = 0; k1 < nex; k1++)
      for (int k2 = 0; k2 < nex; k2++)
      {
        //           std::pair<int,int> temp(ocp[k1],uno[k2]);
        std::pair<int, int> temp(pos[k1], uno[k2]);
        if (find(pairs_local.begin(), pairs_local.end(), temp) == pairs_local.end()) //pair is new
          pairs_local.push_back(temp);
      }
  }
  pairs.resize(pairs_local.size());
  int* first  = pairs.data(0);
  int* second = pairs.data(1);

  for (size_t i = 0; i < pairs_local.size(); i++)
  {
    first[i]  = pairs_local[i].first;
    second[i] = pairs_local[i].second;
  }


  app_log() << "Number of terms in pairs array: " << pairs.size() << std::endl;
  ndets_per_excitation_level.resize(nex_max + 1, 0);
  //reorder configs and det data
  std::vector<size_t> det_idx_order;           // old indices in new order
  std::vector<size_t> det_idx_reverse(nci, 0); // new indices in old order

  // populate data, ordered by exc. lvl.
  // make mapping from new to old det idx
  std::vector<int> data_local;
  data_local.clear();
  data.clear();

  for (const auto& [nex, det_idx_old] : sortMap)
  {
    data_local.insert(data_local.end(), dataMap[nex].begin(), dataMap[nex].end());
    det_idx_order.insert(det_idx_order.end(), det_idx_old.begin(), det_idx_old.end());
    ndets_per_excitation_level[nex] = det_idx_old.size();
  }
  data.resize(data_local.size());
  for (size_t i = 0; i < data_local.size(); i++)
    data[i] = data_local[i];

  assert(det_idx_order.size() == nci);

  // make reverse mapping (old to new) and reorder confgList by exc. lvl.
  configlist_sorted.resize(nci);
  sign.resize(nci);
  for (size_t i = 0; i < nci; i++)
  {
    det_idx_reverse[det_idx_order[i]] = i;
    configlist_sorted[i]              = configlist_unsorted[det_idx_order[i]];
    sign[i]                           = tmp_sign[det_idx_order[i]];
  }

  auto& refdet_occup_ref(*refdet_occup);
  refdet_occup_ref.resize(NumPtcls);
  for (size_t i = 0; i < NumPtcls; i++)
    refdet_occup_ref[i] = configlist_unsorted[ReferenceDeterminant].occup[i];

  {
    ScopedTimer local_timer(transferH2D_timer);
    sign.updateTo();
    pairs.updateTo();
    data.updateTo();
    refdet_occup_ref.updateTo();
  }
  // update C2nodes for new det ordering
  C2nodes_sorted.resize(C2nodes_unsorted.size());
  for (int i = 0; i < C2nodes_unsorted.size(); i++)
    C2nodes_sorted[i] = det_idx_reverse[C2nodes_unsorted[i]];

  /*
       std::cout <<"ref: " <<ref << std::endl;
       std::cout <<"list: " << std::endl;
       for(int i=0; i<confgList.size(); i++)
         std::cout <<confgList[i] << std::endl;

       std::cout <<"pairs: " << std::endl;
       for(int i=0; i<pairs.size(); i++)
         std::cout <<pairs[i].first <<"   " <<pairs[i].second << std::endl;
  */

  // make sure internal objects depending on the number of unique determinants are resized
  resize();
}

void MultiDiracDeterminant::evaluateForWalkerMove(const ParticleSet& P, bool fromScratch)
{
  ScopedTimer local_timer(evalWalker_timer);
  if (fromScratch)
  {
    ///Force host view as no implementation of evaluate_notranspose
    Matrix<ValueType> psiM_host_view(psiM.data(), psiM.rows(), psiM.cols());
    Matrix<GradType> dpsiM_host_view(dpsiM.data(), dpsiM.rows(), dpsiM.cols());
    Matrix<ValueType> d2psiM_host_view(d2psiM.data(), d2psiM.rows(), d2psiM.cols());
    Phi->evaluate_notranspose(P, FirstIndex, LastIndex, psiM_host_view, dpsiM_host_view, d2psiM_host_view);
    {
      ScopedTimer local_timer(transferH2D_timer);
      psiM.updateTo();
      dpsiM.updateTo();
      //d2psiM.updateTo();
    }
  }


  const auto& confgList = *ciConfigList;

  {
    ScopedTimer local_timer(inverse_timer);
    auto it(confgList[ReferenceDeterminant].occup.begin());
    for (size_t i = 0; i < NumPtcls; i++)
    {
      for (size_t j = 0; j < NumPtcls; j++)
        psiMinv(j, i) = psiM(j, *it);
      it++;
    }

    for (size_t i = 0; i < NumPtcls; i++)
      for (size_t j = 0; j < NumOrbitals; j++)
        TpsiM(j, i) = psiM(i, j);

    std::complex<RealType> logValueRef;
    InvertWithLog(psiMinv.data(), NumPtcls, NumPtcls, WorkSpace.data(), Pivot.data(), logValueRef);
    log_value_ref_det_ = logValueRef;
  }

  const RealType detsign = (*DetSigns)[ReferenceDeterminant];
  buildTableMatrix_calculateRatios(ReferenceDeterminant, psiMinv, TpsiM, *detData, *uniquePairs, *DetSigns,
                                   table_matrix, ratios_to_ref_);
  ///Pinning ratios_to_ref_ to the device.


  for (size_t iat = 0; iat < NumPtcls; iat++)
  {
    auto it(confgList[ReferenceDeterminant].occup.begin());
    GradType gradRatio;
    ValueType ratioLapl = 0.0;
    for (size_t i = 0; i < NumPtcls; i++)
    {
      gradRatio += psiMinv(i, iat) * dpsiM(iat, *it);
      ratioLapl += psiMinv(i, iat) * d2psiM(iat, *it);
      it++;
    }
    lapls(ReferenceDeterminant, iat) = ratioLapl;
    for (size_t idim = 0; idim < OHMMS_DIM; idim++)
    {
      dpsiMinv = psiMinv;
      it       = confgList[ReferenceDeterminant].occup.begin();
      for (size_t i = 0; i < NumPtcls; i++)
        psiV_temp[i] = dpsiM(iat, *(it++))[idim];
      InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, iat, gradRatio[idim]);
      for (size_t i = 0; i < NumOrbitals; i++)
        TpsiM(i, iat) = dpsiM(iat, i)[idim];
      buildTableMatrix_calculateGradRatios(ReferenceDeterminant, dpsiMinv, TpsiM, *detData, *uniquePairs, *DetSigns,
                                           gradRatio[idim], table_matrix, idim, iat, grads);
    }
    dpsiMinv = psiMinv;
    it       = confgList[ReferenceDeterminant].occup.begin();
    for (size_t i = 0; i < NumPtcls; i++)
      psiV_temp[i] = d2psiM(iat, *(it++));
    InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, iat, ratioLapl);
    for (size_t i = 0; i < NumOrbitals; i++)
      TpsiM(i, iat) = d2psiM(iat, i);
    buildTableMatrix_calculateRatiosValueMatrixOneParticle(ReferenceDeterminant, dpsiMinv, TpsiM, *detData,
                                                           *uniquePairs, *DetSigns, table_matrix, iat, lapls);
    // restore matrix
    for (size_t i = 0; i < NumOrbitals; i++)
      TpsiM(i, iat) = psiM(iat, i);
  }

  {
    ScopedTimer local_timer(transferH2D_timer);
    ratios_to_ref_.updateTo();
    psiMinv.updateTo();
    TpsiM.updateTo();
  }

  psiMinv_temp = psiMinv;
}

void MultiDiracDeterminant::evaluateForWalkerMoveWithSpin(const ParticleSet& P, bool fromScratch)
{
  ScopedTimer local_timer(evalWalker_timer);
  if (fromScratch)
  {
    ///Force host view as no implementation of evaluate_notranspose
    Matrix<ValueType> psiM_host_view(psiM.data(), psiM.rows(), psiM.cols());
    Matrix<GradType> dpsiM_host_view(dpsiM.data(), dpsiM.rows(), dpsiM.cols());
    Matrix<ValueType> d2psiM_host_view(d2psiM.data(), d2psiM.rows(), d2psiM.cols());
    Phi->evaluate_notranspose_spin(P, FirstIndex, LastIndex, psiM_host_view, dpsiM_host_view, d2psiM_host_view,
                                   dspin_psiM);
    {
      ScopedTimer local_timer(transferH2D_timer);
      psiM.updateTo();
      dpsiM.updateTo();
    }
  }


  const auto& confgList = *ciConfigList;
  std::complex<RealType> logValueRef;

  {
    ScopedTimer local_timer(inverse_timer);
    auto it(confgList[ReferenceDeterminant].occup.begin());
    for (size_t i = 0; i < NumPtcls; i++)
    {
      for (size_t j = 0; j < NumPtcls; j++)
        psiMinv(j, i) = psiM(j, *it);
      it++;
    }
    for (size_t i = 0; i < NumPtcls; i++)
    {
      for (size_t j = 0; j < NumOrbitals; j++)
        TpsiM(j, i) = psiM(i, j);
    }
    InvertWithLog(psiMinv.data(), NumPtcls, NumPtcls, WorkSpace.data(), Pivot.data(), logValueRef);
    log_value_ref_det_ = logValueRef;
  } ///Stop inverse_timerScop
  const RealType detsign = (*DetSigns)[ReferenceDeterminant];
  buildTableMatrix_calculateRatios(ReferenceDeterminant, psiMinv, TpsiM, *detData, *uniquePairs, *DetSigns,
                                   table_matrix, ratios_to_ref_);

  for (size_t iat = 0; iat < NumPtcls; iat++)
  {
    auto it(confgList[ReferenceDeterminant].occup.begin());
    GradType gradRatio;
    ValueType ratioLapl     = 0.0;
    ValueType spingradRatio = 0.0;
    for (size_t i = 0; i < NumPtcls; i++)
    {
      gradRatio += psiMinv(i, iat) * dpsiM(iat, *it);
      ratioLapl += psiMinv(i, iat) * d2psiM(iat, *it);
      spingradRatio += psiMinv(i, iat) * dspin_psiM(iat, *it);
      it++;
    }
    lapls(ReferenceDeterminant, iat)     = ratioLapl;
    spingrads(ReferenceDeterminant, iat) = spingradRatio;
    for (size_t idim = 0; idim < OHMMS_DIM; idim++)
    {
      dpsiMinv = psiMinv;
      it       = confgList[ReferenceDeterminant].occup.begin();
      for (size_t i = 0; i < NumPtcls; i++)
        psiV_temp[i] = dpsiM(iat, *(it++))[idim];
      InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, iat, gradRatio[idim]);
      for (size_t i = 0; i < NumOrbitals; i++)
        TpsiM(i, iat) = dpsiM(iat, i)[idim];
      buildTableMatrix_calculateGradRatios(ReferenceDeterminant, dpsiMinv, TpsiM, *detData, *uniquePairs, *DetSigns,
                                           gradRatio[idim], table_matrix, idim, iat, grads);
    }
    dpsiMinv = psiMinv;
    it       = confgList[ReferenceDeterminant].occup.begin();
    for (size_t i = 0; i < NumPtcls; i++)
      psiV_temp[i] = d2psiM(iat, *(it++));
    InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, iat, ratioLapl);
    for (size_t i = 0; i < NumOrbitals; i++)
      TpsiM(i, iat) = d2psiM(iat, i);
    buildTableMatrix_calculateRatiosValueMatrixOneParticle(ReferenceDeterminant, dpsiMinv, TpsiM, *detData,
                                                           *uniquePairs, *DetSigns, table_matrix, iat, lapls);

    //Adding the spin gradient
    dpsiMinv = psiMinv;
    it       = confgList[ReferenceDeterminant].occup.begin();
    for (size_t i = 0; i < NumPtcls; i++)
      psiV_temp[i] = dspin_psiM(iat, *(it++));
    InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, iat, spingradRatio);
    for (size_t i = 0; i < NumOrbitals; i++)
      TpsiM(i, iat) = dspin_psiM(iat, i);
    buildTableMatrix_calculateRatiosValueMatrixOneParticle(ReferenceDeterminant, dpsiMinv, TpsiM, *detData,
                                                           *uniquePairs, *DetSigns, table_matrix, iat, spingrads);

    // restore matrix
    for (size_t i = 0; i < NumOrbitals; i++)
      TpsiM(i, iat) = psiM(iat, i);
  }

  {
    ScopedTimer local_timer(transferH2D_timer);
    ratios_to_ref_.updateTo();
    psiMinv.updateTo();
    TpsiM.updateTo();
  }
  psiMinv_temp = psiMinv;
}


MultiDiracDeterminant::LogValueType MultiDiracDeterminant::updateBuffer(ParticleSet& P,
                                                                        WFBufferType& buf,
                                                                        bool fromscratch)
{
  assert(P.isSpinor() == is_spinor_);
  if (is_spinor_)
    evaluateForWalkerMoveWithSpin(P, (fromscratch || UpdateMode == ORB_PBYP_RATIO));
  else
    evaluateForWalkerMove(P, (fromscratch || UpdateMode == ORB_PBYP_RATIO));
  buf.put(psiM.first_address(), psiM.last_address());
  buf.put(FirstAddressOfdpsiM, LastAddressOfdpsiM);
  buf.put(d2psiM.first_address(), d2psiM.last_address());
  buf.put(psiMinv.first_address(), psiMinv.last_address());
  buf.put(log_value_ref_det_);
  buf.put(ratios_to_ref_.first_address(), ratios_to_ref_.last_address());
  buf.put(FirstAddressOfGrads, LastAddressOfGrads);
  buf.put(lapls.first_address(), lapls.last_address());
  if (is_spinor_)
  {
    buf.put(dspin_psiM.first_address(), dspin_psiM.last_address());
    buf.put(spingrads.first_address(), spingrads.last_address());
  }
  return 1.0;
}

void MultiDiracDeterminant::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  assert(P.isSpinor() == is_spinor_);
  buf.get(psiM.first_address(), psiM.last_address());
  buf.get(FirstAddressOfdpsiM, LastAddressOfdpsiM);
  buf.get(d2psiM.first_address(), d2psiM.last_address());
  buf.get(psiMinv.first_address(), psiMinv.last_address());
  buf.get(log_value_ref_det_);
  buf.get(ratios_to_ref_.first_address(), ratios_to_ref_.last_address());
  buf.get(FirstAddressOfGrads, LastAddressOfGrads);
  buf.get(lapls.first_address(), lapls.last_address());
  if (is_spinor_)
  {
    buf.get(dspin_psiM.first_address(), dspin_psiM.last_address());
    buf.get(spingrads.first_address(), spingrads.last_address());
  }
  // only used with ORB_PBYP_ALL,
  psiMinv_temp = psiMinv;
  int n1       = psiM.extent(0);
  int n2       = psiM.extent(1);
  for (int i = 0; i < n1; i++)
    for (int j = 0; j < n2; j++)
      TpsiM(j, i) = psiM(i, j);
}

/** move was accepted, update the real container
*/
void MultiDiracDeterminant::acceptMove(ParticleSet& P, int iat, bool safe_to_delay)
{
  const int WorkingIndex = iat - FirstIndex;
  assert(WorkingIndex >= 0 && WorkingIndex < LastIndex - FirstIndex);
  assert(P.isSpinor() == is_spinor_);
  log_value_ref_det_ += convertValueToLog(curRatio);
  curRatio = ValueType(1);
  switch (UpdateMode)
  {
  case ORB_PBYP_RATIO:
    psiMinv = psiMinv_temp;
    // Ye: During acceptMove and restore, TpsiM is updated on the host, thus need to update the device copy.
    // Ideally, this should be done directly on the device.
    // However, acceptMove/restore are shared by both single walker and batched APIs.
    // So the data motion must be kept consistently in all implementations.
    // Right now in batched APIs, ratio and ratioGad implementation also doesn't have the same data motion.
    // Thus also need a fix.
    for (int i = 0; i < NumOrbitals; i++)
      TpsiM(i, WorkingIndex) = psiV[i];
    std::copy(psiV.begin(), psiV.end(), psiM[iat - FirstIndex]);
    std::copy(new_ratios_to_ref_.begin(), new_ratios_to_ref_.end(), ratios_to_ref_.begin());
    {
      ScopedTimer local_timer(transferH2D_timer);
      ratios_to_ref_.updateTo();
      TpsiM.updateTo();
      psiMinv.updateTo();
      psiM.updateTo();
      dpsiM.updateTo();
    }
    break;
  case ORB_PBYP_PARTIAL:
    psiMinv = psiMinv_temp;
    for (int i = 0; i < NumOrbitals; i++)
      TpsiM(i, WorkingIndex) = psiV[i];
    std::copy(new_ratios_to_ref_.begin(), new_ratios_to_ref_.end(), ratios_to_ref_.begin());
    std::copy(psiV.begin(), psiV.end(), psiM[WorkingIndex]);
    std::copy(dpsiV.begin(), dpsiV.end(), dpsiM[WorkingIndex]);
    std::copy(d2psiV.begin(), d2psiV.end(), d2psiM[WorkingIndex]);
    {
      ScopedTimer local_timer(transferH2D_timer);
      ratios_to_ref_.updateTo();
      TpsiM.updateTo();
      psiMinv.updateTo();
      psiM.updateTo();
      dpsiM.updateTo();
    }
    if (is_spinor_)
      std::copy(dspin_psiV.begin(), dspin_psiV.end(), dspin_psiM[WorkingIndex]);
    break;
  default:
    psiMinv = psiMinv_temp;
    for (int i = 0; i < NumOrbitals; i++)
      TpsiM(i, WorkingIndex) = psiV[i];

    std::copy(new_ratios_to_ref_.begin(), new_ratios_to_ref_.end(), ratios_to_ref_.begin());
    std::copy(new_grads.begin(), new_grads.end(), grads.begin());
    std::copy(new_lapls.begin(), new_lapls.end(), lapls.begin());
    std::copy(psiV.begin(), psiV.end(), psiM[WorkingIndex]);
    std::copy(dpsiV.begin(), dpsiV.end(), dpsiM[WorkingIndex]);
    std::copy(d2psiV.begin(), d2psiV.end(), d2psiM[WorkingIndex]);
    {
      ScopedTimer local_timer(transferH2D_timer);
      ratios_to_ref_.updateTo();
      TpsiM.updateTo();
      psiMinv.updateTo();
      psiM.updateTo();
      dpsiM.updateTo();
    }
    if (is_spinor_)
    {
      std::copy(new_spingrads.begin(), new_spingrads.end(), spingrads.begin());
      std::copy(dspin_psiV.begin(), dspin_psiV.end(), dspin_psiM[WorkingIndex]);
    }
    break;
  }
}

/** move was rejected. copy the real container to the temporary to move on
*/
void MultiDiracDeterminant::restore(int iat)
{
  const int WorkingIndex = iat - FirstIndex;
  assert(WorkingIndex >= 0 && WorkingIndex < LastIndex - FirstIndex);
  psiMinv_temp = psiMinv;
  for (int i = 0; i < NumOrbitals; i++)
    TpsiM(i, WorkingIndex) = psiM(WorkingIndex, i);
  {
    ScopedTimer local_timer(transferH2D_timer);
    TpsiM.updateTo();
  }
  curRatio = ValueType(1);
  /*
      switch(UpdateMode)
      {
        case ORB_PBYP_RATIO:
          psiMinv_temp = psiMinv;
          for(int i=0; i<NumOrbitals; i++)
            TpsiM(i,WorkingIndex) = psiM(WorkingIndex,i);
          break;
        case ORB_PBYP_PARTIAL:
          psiMinv_temp = psiMinv;
          for(int i=0; i<NumOrbitals; i++)
            TpsiM(i,WorkingIndex) = psiM(WorkingIndex,i);
          break;
        default:
          break;
      }
  */
}

void MultiDiracDeterminant::mw_accept_rejectMove(const RefVectorWithLeader<MultiDiracDeterminant>& wfc_list,
                                                 const RefVectorWithLeader<ParticleSet>& p_list,
                                                 int iat,
                                                 const std::vector<bool>& isAccepted)
{
  // TODO to be expanded to serve offload needs without relying on calling acceptMove and restore
  for (int iw = 0; iw < wfc_list.size(); iw++)
    if (isAccepted[iw])
      wfc_list[iw].acceptMove(p_list[iw], iat, false);
    else
      wfc_list[iw].restore(iat);
}

// this has been fixed
MultiDiracDeterminant::MultiDiracDeterminant(const MultiDiracDeterminant& s)
    : WaveFunctionComponent(s),
      inverse_timer(s.inverse_timer),
      buildTable_timer(s.buildTable_timer),
      table2ratios_timer(s.table2ratios_timer),
      evalWalker_timer(s.evalWalker_timer),
      evalOrbValue_timer(s.evalOrbValue_timer),
      evalOrbVGL_timer(s.evalOrbVGL_timer),
      updateInverse_timer(s.updateInverse_timer),
      calculateRatios_timer(s.calculateRatios_timer),
      calculateGradRatios_timer(s.calculateGradRatios_timer),
      updateRatios_timer(s.updateRatios_timer),
      evaluateDetsForPtclMove_timer(s.evaluateDetsForPtclMove_timer),
      evaluateDetsAndGradsForPtclMove_timer(s.evaluateDetsAndGradsForPtclMove_timer),
      evaluateGrads_timer(s.evaluateGrads_timer),
      offload_timer(s.offload_timer),
      transferH2D_timer(s.transferH2D_timer),
      transferD2H_timer(s.transferD2H_timer),
      lookup_tbl(s.lookup_tbl),
      Phi(s.Phi->makeClone()),
      NumOrbitals(Phi->getOrbitalSetSize()),
      FirstIndex(s.FirstIndex),
      NumPtcls(s.NumPtcls),
      LastIndex(s.LastIndex),
      ciConfigList(s.ciConfigList),
      refdet_occup(s.refdet_occup),
      is_spinor_(s.is_spinor_),
      detData(s.detData),
      uniquePairs(s.uniquePairs),
      DetSigns(s.DetSigns),
      ndets_per_excitation_level_(s.ndets_per_excitation_level_)
{
  resize();
}

std::unique_ptr<SPOSet> MultiDiracDeterminant::clonePhi() const { return Phi->makeClone(); }

std::unique_ptr<WaveFunctionComponent> MultiDiracDeterminant::makeClone(ParticleSet& tqp) const
{
  APP_ABORT(" Illegal action. Cannot use MultiDiracDeterminant::makeClone");
  return std::unique_ptr<MultiDiracDeterminant>();
}

/** constructor
 *@param spos the single-particle orbital set
 *@param first index of the first particle
 *@param spinor flag to determinane if spin arrays need to be resized and used
 */
MultiDiracDeterminant::MultiDiracDeterminant(std::unique_ptr<SPOSet>&& spos, bool spinor, int first, int nel)
    : inverse_timer(createGlobalTimer(getClassName() + "::invertRefDet")),
      buildTable_timer(createGlobalTimer(getClassName() + "::buildTable")),
      table2ratios_timer(createGlobalTimer(getClassName() + "::table2ratios")),
      evalWalker_timer(createGlobalTimer(getClassName() + "::evalWalker")),
      evalOrbValue_timer(createGlobalTimer(getClassName() + "::evalOrbValue")),
      evalOrbVGL_timer(createGlobalTimer(getClassName() + "::evalOrbVGL")),
      updateInverse_timer(createGlobalTimer(getClassName() + "::updateRefDetInv")),
      calculateRatios_timer(createGlobalTimer(getClassName() + "::calcRatios")),
      calculateGradRatios_timer(createGlobalTimer(getClassName() + "::calcGradRatios")),
      updateRatios_timer(createGlobalTimer(getClassName() + "::updateRatios")),
      evaluateDetsForPtclMove_timer(createGlobalTimer(getClassName() + "::evaluateDet")),
      evaluateDetsAndGradsForPtclMove_timer(createGlobalTimer(getClassName() + "::evaluateDetAndGrad")),
      evaluateGrads_timer(createGlobalTimer(getClassName() + "::evaluateGrad")),
      offload_timer(createGlobalTimer(getClassName() + "::offload")),
      transferH2D_timer(createGlobalTimer(getClassName() + "::transferH2D")),
      transferD2H_timer(createGlobalTimer(getClassName() + "::transferD2H")),
      Phi(std::move(spos)),
      NumOrbitals(Phi->getOrbitalSetSize()),
      FirstIndex(first),
      NumPtcls(nel),
      LastIndex(first + nel),
      is_spinor_(spinor)
{
  ciConfigList                = std::make_shared<std::vector<ci_configuration2>>();
  refdet_occup                = std::make_shared<OffloadVector<size_t>>();
  detData                     = std::make_shared<OffloadVector<int>>();
  uniquePairs                 = std::make_shared<VectorSoaContainer<int, 2, OffloadPinnedAllocator<int>>>();
  DetSigns                    = std::make_shared<OffloadVector<RealType>>();
  ndets_per_excitation_level_ = std::make_shared<std::vector<int>>();
}

///default destructor
MultiDiracDeterminant::~MultiDiracDeterminant() = default;

void MultiDiracDeterminant::registerData(ParticleSet& P, WFBufferType& buf)
{
  assert(P.isSpinor() == is_spinor_);

  //extra pointers
  FirstAddressOfGrads = &(grads(0, 0)[0]);
  LastAddressOfGrads  = FirstAddressOfGrads + NumPtcls * DIM * getNumDets();
  FirstAddressOfdpsiM = &(dpsiM(0, 0)[0]);
  LastAddressOfdpsiM  = FirstAddressOfdpsiM + NumPtcls * NumOrbitals * DIM;

  //add the data:
  buf.add(psiM.first_address(), psiM.last_address());
  buf.add(FirstAddressOfdpsiM, LastAddressOfdpsiM);
  buf.add(d2psiM.first_address(), d2psiM.last_address());
  buf.add(psiMinv.first_address(), psiMinv.last_address());
  buf.add(log_value_ref_det_);
  buf.add(ratios_to_ref_.first_address(), ratios_to_ref_.last_address());
  buf.add(FirstAddressOfGrads, LastAddressOfGrads);
  buf.add(lapls.first_address(), lapls.last_address());
  if (is_spinor_)
  {
    buf.add(dspin_psiM.first_address(), dspin_psiM.last_address());
    buf.add(spingrads.first_address(), spingrads.last_address());
  }
}

void MultiDiracDeterminant::createResource(ResourceCollection& collection) const
{
  collection.addResource(std::make_unique<MultiDiracDetMultiWalkerResource>());
}

void MultiDiracDeterminant::acquireResource(ResourceCollection& collection,
                                            const RefVectorWithLeader<MultiDiracDeterminant>& wfc_list) const
{
  auto& wfc_leader          = wfc_list.getCastedLeader<MultiDiracDeterminant>();
  wfc_leader.mw_res_handle_ = collection.lendResource<MultiDiracDetMultiWalkerResource>();
  auto& mw_res              = wfc_leader.mw_res_handle_.getResource();

  const size_t nw = wfc_list.size();
  mw_res.resizeConstants(nw);

  auto& psiV_temp_deviceptr_list    = mw_res.psiV_temp_deviceptr_list;
  auto& psiMinv_temp_deviceptr_list = mw_res.psiMinv_temp_deviceptr_list;
  auto& dpsiMinv_deviceptr_list     = mw_res.dpsiMinv_deviceptr_list;
  auto& workV1_deviceptr_list       = mw_res.workV1_deviceptr_list;
  auto& workV2_deviceptr_list       = mw_res.workV2_deviceptr_list;

  auto& psiV_deviceptr_list    = mw_res.psiV_deviceptr_list;
  auto& dpsiV_deviceptr_list   = mw_res.dpsiV_deviceptr_list;
  auto& TpsiM_deviceptr_list   = mw_res.TpsiM_deviceptr_list;
  auto& psiM_deviceptr_list    = mw_res.psiM_deviceptr_list;
  auto& psiMinv_deviceptr_list = mw_res.psiMinv_deviceptr_list;
  auto& dpsiM_deviceptr_list   = mw_res.dpsiM_deviceptr_list;

  psiV_temp_deviceptr_list.resize(nw);
  psiMinv_temp_deviceptr_list.resize(nw);
  dpsiMinv_deviceptr_list.resize(nw);
  workV1_deviceptr_list.resize(nw);
  workV2_deviceptr_list.resize(nw);

  psiV_deviceptr_list.resize(nw);
  dpsiV_deviceptr_list.resize(nw);
  TpsiM_deviceptr_list.resize(nw);
  psiM_deviceptr_list.resize(nw);
  psiMinv_deviceptr_list.resize(nw);
  dpsiM_deviceptr_list.resize(nw);

  for (size_t iw = 0; iw < nw; iw++)
  {
    auto& det                       = wfc_list.getCastedElement<MultiDiracDeterminant>(iw);
    psiV_temp_deviceptr_list[iw]    = det.psiV_temp.device_data();
    psiMinv_temp_deviceptr_list[iw] = det.psiMinv_temp.device_data();
    dpsiMinv_deviceptr_list[iw]     = det.dpsiMinv.device_data();
    workV1_deviceptr_list[iw]       = det.workV1.device_data();
    workV2_deviceptr_list[iw]       = det.workV2.device_data();

    psiV_deviceptr_list[iw]    = det.psiV.device_data();
    dpsiV_deviceptr_list[iw]   = det.dpsiV.device_data();
    TpsiM_deviceptr_list[iw]   = det.TpsiM.device_data();
    psiM_deviceptr_list[iw]    = det.psiM.device_data();
    psiMinv_deviceptr_list[iw] = det.psiMinv.device_data();
    dpsiM_deviceptr_list[iw]   = det.dpsiM.device_data();
  }

  psiV_temp_deviceptr_list.updateTo();
  psiMinv_temp_deviceptr_list.updateTo();
  dpsiMinv_deviceptr_list.updateTo();
  workV1_deviceptr_list.updateTo();
  workV2_deviceptr_list.updateTo();

  psiV_deviceptr_list.updateTo();
  dpsiV_deviceptr_list.updateTo();
  TpsiM_deviceptr_list.updateTo();
  psiM_deviceptr_list.updateTo();
  psiMinv_deviceptr_list.updateTo();
  dpsiM_deviceptr_list.updateTo();
}

void MultiDiracDeterminant::releaseResource(ResourceCollection& collection,
                                            const RefVectorWithLeader<MultiDiracDeterminant>& wfc_list) const
{
  auto& wfc_leader = wfc_list.getCastedLeader<MultiDiracDeterminant>();
  collection.takebackResource(wfc_leader.mw_res_handle_);
}

///reset the size: with the number of particles and number of orbtials
void MultiDiracDeterminant::resize()
{
  const int nel = NumPtcls;
  assert(NumPtcls > 0);
  const int NumDets = getNumDets();
  assert(NumDets > 0);

  psiV_temp.resize(nel);
  psiV.resize(NumOrbitals);
  dpsiV.resize(NumOrbitals);
  d2psiV.resize(NumOrbitals);
  psiM.resize(nel, NumOrbitals);
  dpsiM.resize(nel, NumOrbitals);
  d2psiM.resize(nel, NumOrbitals);
  TpsiM.resize(NumOrbitals, nel);
  psiMinv.resize(nel, nel);
  dpsiMinv.resize(nel, nel);
  psiMinv_temp.resize(nel, nel);
  //scratch spaces: stateless
  WorkSpace.resize(std::max(nel, NumDets));
  Pivot.resize(nel);
  workV1.resize(nel);
  workV2.resize(nel);
  ratios_to_ref_.resize(NumDets);
  new_ratios_to_ref_.resize(NumDets);
  grads.resize(NumDets, nel);
  new_grads.resize(NumDets, nel);
  lapls.resize(NumDets, nel);
  new_lapls.resize(NumDets, nel);
  table_matrix.resize(NumOrbitals, NumOrbitals);
  det_calculator_.resize(nel);

  if (is_spinor_)
  {
    dspin_psiV.resize(NumOrbitals);
    dspin_psiM.resize(nel, NumOrbitals);
    spingrads.resize(NumDets, nel);
    new_spingrads.resize(NumDets, nel);
  }
}

void MultiDiracDeterminant::buildOptVariables(std::vector<size_t>& C2node)
{
  if (!isOptimizable())
    return;

  const size_t nel = NumPtcls;
  const size_t nmo = NumOrbitals;
  //a vector in which the element's index value correspond to Molecular Orbitals.
  //The element value at an index indicates how many times an electron is excited from or to that orbital in the Multi-Slater expansion i.e the indices with non-zero elements are active space orbitals
  std::vector<int> occupancy_vector(nmo, 0);

  // Function to fill occupancy_vectors and also return number of unique determinants
  const size_t unique_dets = build_occ_vec(*detData, nel, nmo, occupancy_vector);

  // When calculating the parameter derivative of the Multi-Slater component of the wavefunction, each unique deterimant can contribute multiple times.
  // The lookup_tbls are used so that a parameter derivative of a unique determinant is only done once and then scaled according to how many times it appears in the Multi-Slater expansion
  lookup_tbl.resize(unique_dets);
  //construct lookup table
  for (int i(0); i < C2node.size(); i++)
  {
    lookup_tbl[C2node[i]].push_back(i);
  }

  // create active rotation parameter indices
  std::vector<std::pair<int, int>> m_act_rot_inds;
  std::vector<std::pair<int, int>> other_rot_inds;

  for (int i = 0; i < nmo; i++)
    for (int j = i + 1; j < nmo; j++)
    {
      bool core_i(!occupancy_vector[i] and i <= nel - 1); // true if orbital i is a 'core' orbital
      bool core_j(!occupancy_vector[j] and j <= nel - 1); // true if orbital j is a 'core' orbital
      bool virt_i(!occupancy_vector[i] and i > nel - 1);  // true if orbital i is a 'virtual' orbital
      bool virt_j(!occupancy_vector[j] and j > nel - 1);  // true if orbital j is a 'virtual' orbital
      if (!((core_i and core_j) or (virt_i and virt_j)))
      {
        m_act_rot_inds.push_back(
            std::pair<
                int,
                int>(i,
                     j)); // orbital rotation parameter accepted as long as rotation isn't core-core or virtual-virtual
      }
      else
      {
        other_rot_inds.push_back(std::pair<int, int>(i, j));
      }
    }

  std::vector<std::pair<int, int>> full_rot_inds;

  // Copy the adjustable rotations first
  full_rot_inds = m_act_rot_inds;
  // Add the other rotations at the end
  full_rot_inds.insert(full_rot_inds.end(), other_rot_inds.begin(), other_rot_inds.end());


#ifndef QMC_COMPLEX
  RotatedSPOs* rot_spo = dynamic_cast<RotatedSPOs*>(Phi.get());
  if (rot_spo)
    rot_spo->buildOptVariables(m_act_rot_inds, full_rot_inds);
#endif
}

int MultiDiracDeterminant::build_occ_vec(const OffloadVector<int>& data,
                                         const size_t nel,
                                         const size_t nmo,
                                         std::vector<int>& occ_vec)
{
  size_t it = 0;
  int count = 0; //number of determinants
  while (it < data.size())
  {
    int k = data[it]; // number of excitations with respect to the reference matrix
    if (count == 0)
    {
      it += 3 * k + 1;
      count++;
    }
    else
    {
      for (int i = 0; i < k; i++)
      {
        //for determining active orbitals
        occ_vec[data[it + 1 + i]]++;
        occ_vec[data[it + 1 + k + i]]++;
      }
      it += 3 * k + 1;
      count++;
    }
  }
  return count;
}


void MultiDiracDeterminant::evaluateDerivatives(ParticleSet& P,
                                                const opt_variables_type& optvars,
                                                Vector<ValueType>& dlogpsi,
                                                Vector<ValueType>& dhpsioverpsi,
                                                const MultiDiracDeterminant& pseudo_dn,
                                                const ValueType& psiCurrent,
                                                const std::vector<ValueType>& Coeff,
                                                const std::vector<size_t>& C2node_up,
                                                const std::vector<size_t>& C2node_dn)
{
  if (!isOptimizable())
    return;

  const OffloadVector<ValueType>& detValues_up = getRatiosToRefDet();
  const OffloadVector<ValueType>& detValues_dn = pseudo_dn.getRatiosToRefDet();
  const Matrix<GradType>& grads_up             = grads;
  const Matrix<GradType>& grads_dn             = pseudo_dn.grads;
  const Matrix<ValueType>& lapls_up            = lapls;
  const Matrix<ValueType>& lapls_dn            = pseudo_dn.lapls;
  const OffloadMatrix<ValueType>& M_up         = psiM;
  const OffloadMatrix<ValueType>& M_dn         = pseudo_dn.psiM;
  const OffloadMatrix<ValueType>& Minv_up      = psiMinv;
  const OffloadMatrix<ValueType>& Minv_dn      = pseudo_dn.psiMinv;
  const OffloadMatrix<GradType>& B_grad        = dpsiM;
  const OffloadMatrix<ValueType>& B_lapl       = d2psiM;
  std::vector<int> detData_local(detData->size());
  for (size_t i = 0; i < detData->size(); i++)
    detData_local[i] = (*detData)[i];


  const size_t N1  = FirstIndex;
  const size_t N2  = pseudo_dn.FirstIndex;
  const size_t NP1 = NumPtcls;
  const size_t NP2 = pseudo_dn.NumPtcls;
  Vector<ValueType> detValues_up_host_view(const_cast<ValueType*>(detValues_up.data()), detValues_up.size());
  Vector<ValueType> detValues_dn_host_view(const_cast<ValueType*>(detValues_dn.data()), detValues_dn.size());
  Matrix<ValueType> M_up_host_view(const_cast<ValueType*>(M_up.data()), M_up.rows(), M_up.cols());
  Matrix<ValueType> M_dn_host_view(const_cast<ValueType*>(M_dn.data()), M_dn.rows(), M_dn.cols());
  Matrix<ValueType> Minv_up_host_view(const_cast<ValueType*>(Minv_up.data()), Minv_up.rows(), Minv_up.cols());
  Matrix<ValueType> Minv_dn_host_view(const_cast<ValueType*>(Minv_dn.data()), Minv_dn.rows(), Minv_dn.cols());
  Matrix<GradType> B_grad_host_view(const_cast<GradType*>(B_grad.data()), B_grad.rows(), B_grad.cols());
  Matrix<ValueType> B_lapl_host_view(const_cast<ValueType*>(B_lapl.data()), B_lapl.rows(), B_lapl.cols());
  Phi->evaluateDerivatives(P, optvars, dlogpsi, dhpsioverpsi, psiCurrent, Coeff, C2node_up, C2node_dn,
                           detValues_up_host_view, detValues_dn_host_view, grads_up, grads_dn, lapls_up, lapls_dn,
                           M_up_host_view, M_dn_host_view, Minv_up_host_view, Minv_dn_host_view, B_grad_host_view,
                           B_lapl_host_view, detData_local, N1, N2, NP1, NP2, lookup_tbl);
}


void MultiDiracDeterminant::evaluateDerivativesWF(ParticleSet& P,
                                                  const opt_variables_type& optvars,
                                                  Vector<ValueType>& dlogpsi,
                                                  const MultiDiracDeterminant& pseudo_dn,
                                                  const PsiValueType& psiCurrent,
                                                  const std::vector<ValueType>& Coeff,
                                                  const std::vector<size_t>& C2node_up,
                                                  const std::vector<size_t>& C2node_dn)
{
  if (!isOptimizable())
    return;

  const OffloadVector<ValueType>& detValues_up = getRatiosToRefDet();
  const OffloadVector<ValueType>& detValues_dn = pseudo_dn.getRatiosToRefDet();
  const OffloadMatrix<ValueType>& M_up         = psiM;
  const OffloadMatrix<ValueType>& M_dn         = pseudo_dn.psiM;
  const OffloadMatrix<ValueType>& Minv_up      = psiMinv;
  const OffloadMatrix<ValueType>& Minv_dn      = pseudo_dn.psiMinv;

  std::vector<int> detData_local(detData->size());
  for (size_t i = 0; i < detData->size(); i++)
    detData_local[i] = (*detData)[i];
  Vector<ValueType> detValues_up_host_view(const_cast<ValueType*>(detValues_up.data()), detValues_up.size());
  Vector<ValueType> detValues_dn_host_view(const_cast<ValueType*>(detValues_dn.data()), detValues_dn.size());
  Matrix<ValueType> M_up_host_view(const_cast<ValueType*>(M_up.data()), M_up.rows(), M_up.cols());
  Matrix<ValueType> M_dn_host_view(const_cast<ValueType*>(M_dn.data()), M_dn.rows(), M_dn.cols());
  Matrix<ValueType> Minv_up_host_view(const_cast<ValueType*>(Minv_up.data()), Minv_up.rows(), Minv_up.cols());
  Matrix<ValueType> Minv_dn_host_view(const_cast<ValueType*>(Minv_dn.data()), Minv_dn.rows(), Minv_dn.cols());
  Phi->evaluateDerivativesWF(P, optvars, dlogpsi, psiCurrent, Coeff, C2node_up, C2node_dn, detValues_up_host_view,
                             detValues_dn_host_view, M_up_host_view, M_dn_host_view, Minv_up_host_view,
                             Minv_dn_host_view, detData_local, lookup_tbl);
}

} // namespace qmcplusplus
