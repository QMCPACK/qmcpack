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
  ReferenceDeterminant = ref_det_id;

  auto& ref                        = configlist_unsorted[ref_det_id];
  auto& configlist_sorted          = *ciConfigList;
  auto& data                       = *detData;
  auto& pairs                      = *uniquePairs;
  auto& sign                       = *DetSigns;
  auto& ndets_per_excitation_level = *ndets_per_excitation_level_;

  const size_t nci = configlist_unsorted.size();

  size_t nex_max = 0;
  std::vector<size_t> pos(NumPtcls);
  std::vector<size_t> ocp(NumPtcls);
  std::vector<size_t> uno(NumPtcls);
  // map key is exc. lvl
  std::map<int, std::vector<int>> dataMap;
  std::map<int, std::vector<int>> sortMap;
  std::vector<RealType> tmp_sign(nci, 0);
  pairs.clear();
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
        if (find(pairs.begin(), pairs.end(), temp) == pairs.end()) //pair is new
          pairs.push_back(temp);
      }
  }
  app_log() << "Number of terms in pairs array: " << pairs.size() << std::endl;
  ndets_per_excitation_level.resize(nex_max + 1, 0);
  //reorder configs and det data
  std::vector<size_t> det_idx_order;           // old indices in new order
  std::vector<size_t> det_idx_reverse(nci, 0); // new indices in old order

  // populate data, ordered by exc. lvl.
  // make mapping from new to old det idx
  data.clear();
  for (const auto& [nex, det_idx_old] : sortMap)
  {
    data.insert(data.end(), dataMap[nex].begin(), dataMap[nex].end());
    det_idx_order.insert(det_idx_order.end(), det_idx_old.begin(), det_idx_old.end());
    ndets_per_excitation_level[nex] = det_idx_old.size();
  }
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
  evalWTimer.start();
  if (fromScratch)
    Phi->evaluate_notranspose(P, FirstIndex, LastIndex, psiM, dpsiM, d2psiM);

  InverseTimer.start();

  const auto& confgList = *ciConfigList;

  //std::vector<int>::iterator it(confgList[ReferenceDeterminant].occup.begin());
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
  InverseTimer.stop();
  const RealType detsign = (*DetSigns)[ReferenceDeterminant];
  BuildDotProductsAndCalculateRatios(ReferenceDeterminant, psiMinv, TpsiM, *detData, *uniquePairs, *DetSigns,
                                     dotProducts, ratios_to_ref_);
  for (size_t iat = 0; iat < NumPtcls; iat++)
  {
    it = confgList[ReferenceDeterminant].occup.begin();
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
      //MultiDiracDeterminant::InverseUpdateByColumn_GRAD(dpsiMinv,dpsiV,workV1,workV2,iat,gradRatio[idim],idim);
      for (size_t i = 0; i < NumOrbitals; i++)
        TpsiM(i, iat) = dpsiM(iat, i)[idim];
      BuildDotProductsAndCalculateRatiosGrads(ReferenceDeterminant, dpsiMinv, TpsiM, *detData, *uniquePairs, *DetSigns,
                                              gradRatio[idim], dotProducts, idim, iat, grads);
    }
    dpsiMinv = psiMinv;
    it       = confgList[ReferenceDeterminant].occup.begin();
    for (size_t i = 0; i < NumPtcls; i++)
      psiV_temp[i] = d2psiM(iat, *(it++));
    InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, iat, ratioLapl);
    //MultiDiracDeterminant::InverseUpdateByColumn(dpsiMinv,d2psiM,workV1,workV2,iat,ratioLapl,confgList[ReferenceDeterminant].occup.begin());
    for (size_t i = 0; i < NumOrbitals; i++)
      TpsiM(i, iat) = d2psiM(iat, i);
    BuildDotProductsAndCalculateRatiosValueMatrixOneParticle(ReferenceDeterminant, dpsiMinv, TpsiM, *detData,
                                                             *uniquePairs, *DetSigns, dotProducts, iat, lapls);
    // restore matrix
    for (size_t i = 0; i < NumOrbitals; i++)
      TpsiM(i, iat) = psiM(iat, i);
  }

  psiMinv_temp = psiMinv;
  evalWTimer.stop();
}

void MultiDiracDeterminant::evaluateForWalkerMoveWithSpin(const ParticleSet& P, bool fromScratch)
{
  evalWTimer.start();
  if (fromScratch)
    Phi->evaluate_notranspose_spin(P, FirstIndex, LastIndex, psiM, dpsiM, d2psiM, dspin_psiM);

  InverseTimer.start();

  const auto& confgList = *ciConfigList;

  //std::vector<int>::iterator it(confgList[ReferenceDeterminant].occup.begin());
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

  std::complex<RealType> logValueRef;
  InvertWithLog(psiMinv.data(), NumPtcls, NumPtcls, WorkSpace.data(), Pivot.data(), logValueRef);
  log_value_ref_det_ = logValueRef;
  InverseTimer.stop();
  const RealType detsign = (*DetSigns)[ReferenceDeterminant];
  BuildDotProductsAndCalculateRatios(ReferenceDeterminant, psiMinv, TpsiM, *detData, *uniquePairs, *DetSigns,
                                     dotProducts, ratios_to_ref_);
  for (size_t iat = 0; iat < NumPtcls; iat++)
  {
    it = confgList[ReferenceDeterminant].occup.begin();
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
      //MultiDiracDeterminant::InverseUpdateByColumn_GRAD(dpsiMinv,dpsiV,workV1,workV2,iat,gradRatio[idim],idim);
      for (size_t i = 0; i < NumOrbitals; i++)
        TpsiM(i, iat) = dpsiM(iat, i)[idim];
      BuildDotProductsAndCalculateRatiosGrads(ReferenceDeterminant, dpsiMinv, TpsiM, *detData, *uniquePairs, *DetSigns,
                                              gradRatio[idim], dotProducts, idim, iat, grads);
    }
    dpsiMinv = psiMinv;
    it       = confgList[ReferenceDeterminant].occup.begin();
    for (size_t i = 0; i < NumPtcls; i++)
      psiV_temp[i] = d2psiM(iat, *(it++));
    InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, iat, ratioLapl);
    //MultiDiracDeterminant::InverseUpdateByColumn(dpsiMinv,d2psiM,workV1,workV2,iat,ratioLapl,confgList[ReferenceDeterminant].occup.begin());
    for (size_t i = 0; i < NumOrbitals; i++)
      TpsiM(i, iat) = d2psiM(iat, i);
    BuildDotProductsAndCalculateRatiosValueMatrixOneParticle(ReferenceDeterminant, dpsiMinv, TpsiM, *detData,
                                                             *uniquePairs, *DetSigns, dotProducts, iat, lapls);

    //Adding the spin gradient
    dpsiMinv = psiMinv;
    it       = confgList[ReferenceDeterminant].occup.begin();
    for (size_t i = 0; i < NumPtcls; i++)
      psiV_temp[i] = dspin_psiM(iat, *(it++));
    InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, iat, spingradRatio);
    for (size_t i = 0; i < NumOrbitals; i++)
      TpsiM(i, iat) = dspin_psiM(iat, i);
    BuildDotProductsAndCalculateRatiosValueMatrixOneParticle(ReferenceDeterminant, dpsiMinv, TpsiM, *detData,
                                                             *uniquePairs, *DetSigns, dotProducts, iat, spingrads);

    // restore matrix
    for (size_t i = 0; i < NumOrbitals; i++)
      TpsiM(i, iat) = psiM(iat, i);
  }
  psiMinv_temp = psiMinv;
  evalWTimer.stop();
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
    for (int i = 0; i < NumOrbitals; i++)
      TpsiM(i, WorkingIndex) = psiV[i];
    std::copy(psiV.begin(), psiV.end(), psiM[iat - FirstIndex]);
    std::copy(new_ratios_to_ref_.begin(), new_ratios_to_ref_.end(), ratios_to_ref_.begin());
    break;
  case ORB_PBYP_PARTIAL:
    psiMinv = psiMinv_temp;
    for (int i = 0; i < NumOrbitals; i++)
      TpsiM(i, WorkingIndex) = psiV[i];
    std::copy(new_ratios_to_ref_.begin(), new_ratios_to_ref_.end(), ratios_to_ref_.begin());
    std::copy(psiV.begin(), psiV.end(), psiM[WorkingIndex]);
    std::copy(dpsiV.begin(), dpsiV.end(), dpsiM[WorkingIndex]);
    std::copy(d2psiV.begin(), d2psiV.end(), d2psiM[WorkingIndex]);
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

// this has been fixed
MultiDiracDeterminant::MultiDiracDeterminant(const MultiDiracDeterminant& s)
    : WaveFunctionComponent(s),
      UpdateTimer(*timer_manager.createTimer(ClassName + "::update")),
      RatioTimer(*timer_manager.createTimer(ClassName + "::ratio")),
      MWRatioTimer(*timer_manager.createTimer(ClassName + "::mwratio")),
      InverseTimer(*timer_manager.createTimer(ClassName + "::inverse")),
      buildTableTimer(*timer_manager.createTimer(ClassName + "::buildTable")),
      readMatTimer(*timer_manager.createTimer(ClassName + "::readMat")),
      evalWTimer(*timer_manager.createTimer(ClassName + "::evalW")),
      evalOrbTimer(*timer_manager.createTimer(ClassName + "::evalOrb")),
      evalOrb1Timer(*timer_manager.createTimer(ClassName + "::evalOrbGrad")),
      readMatGradTimer(*timer_manager.createTimer(ClassName + "::readMatGrad")),
      buildTableGradTimer(*timer_manager.createTimer(ClassName + "::buildTableGrad")),
      ExtraStuffTimer(*timer_manager.createTimer(ClassName + "::RefDetInvUpdate")),
      Phi(s.Phi->makeClone()),
      NumOrbitals(Phi->getOrbitalSetSize()),
      FirstIndex(s.FirstIndex),
      NumPtcls(s.NumPtcls),
      LastIndex(s.LastIndex),
      ciConfigList(s.ciConfigList),
      ReferenceDeterminant(s.ReferenceDeterminant),
      is_spinor_(s.is_spinor_),
      detData(s.detData),
      uniquePairs(s.uniquePairs),
      DetSigns(s.DetSigns),
      ndets_per_excitation_level_(s.ndets_per_excitation_level_)
{
  Optimizable = s.Optimizable;

  resize();
  registerTimers();
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
    : WaveFunctionComponent("MultiDiracDeterminant"),
      UpdateTimer(*timer_manager.createTimer(ClassName + "::update")),
      RatioTimer(*timer_manager.createTimer(ClassName + "::ratio")),
      MWRatioTimer(*timer_manager.createTimer(ClassName + "::mwratio")),
      InverseTimer(*timer_manager.createTimer(ClassName + "::inverse")),
      buildTableTimer(*timer_manager.createTimer(ClassName + "::buildTable")),
      readMatTimer(*timer_manager.createTimer(ClassName + "::readMat")),
      evalWTimer(*timer_manager.createTimer(ClassName + "::evalW")),
      evalOrbTimer(*timer_manager.createTimer(ClassName + "::evalOrb")),
      evalOrb1Timer(*timer_manager.createTimer(ClassName + "::evalOrbGrad")),
      readMatGradTimer(*timer_manager.createTimer(ClassName + "::readMatGrad")),
      buildTableGradTimer(*timer_manager.createTimer(ClassName + "::buildTableGrad")),
      ExtraStuffTimer(*timer_manager.createTimer(ClassName + "::RefDetInvUpdate")),
      Phi(std::move(spos)),
      NumOrbitals(Phi->getOrbitalSetSize()),
      FirstIndex(first),
      NumPtcls(nel),
      LastIndex(first + nel),
      ReferenceDeterminant(0),
      is_spinor_(spinor)
{
  (Phi->isOptimizable() == true) ? Optimizable = true : Optimizable = false;

  ciConfigList                = std::make_shared<std::vector<ci_configuration2>>();
  detData                     = std::make_shared<std::vector<int>>();
  uniquePairs                 = std::make_shared<std::vector<std::pair<int, int>>>();
  DetSigns                    = std::make_shared<std::vector<RealType>>();
  ndets_per_excitation_level_ = std::make_shared<std::vector<int>>();

  registerTimers();
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
  dotProducts.resize(NumOrbitals, NumOrbitals);
  DetCalculator.resize(nel);

  if (is_spinor_)
  {
    dspin_psiV.resize(NumOrbitals);
    dspin_psiM.resize(nel, NumOrbitals);
    spingrads.resize(NumDets, nel);
    new_spingrads.resize(NumDets, nel);
  }
}

void MultiDiracDeterminant::registerTimers()
{
  UpdateTimer.reset();
  RatioTimer.reset();
  MWRatioTimer.reset();
  InverseTimer.reset();
  buildTableTimer.reset();
  readMatTimer.reset();
  evalOrbTimer.reset();
  evalOrb1Timer.reset();
  evalWTimer.reset();
  ExtraStuffTimer.reset();
  buildTableGradTimer.reset();
  readMatGradTimer.reset();
}

void MultiDiracDeterminant::buildOptVariables(std::vector<size_t>& C2node)
{
  if (!Optimizable)
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
    }

  Phi->buildOptVariables(m_act_rot_inds);
}

int MultiDiracDeterminant::build_occ_vec(const std::vector<int>& data,
                                         const size_t nel,
                                         const size_t nmo,
                                         std::vector<int>& occ_vec)
{
  auto it   = data.begin();
  int count = 0; //number of determinants
  while (it != data.end())
  {
    int k = *it; // number of excitations with respect to the reference matrix
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
        occ_vec[*(it + 1 + i)]++;
        occ_vec[*(it + 1 + k + i)]++;
      }
      it += 3 * k + 1;
      count++;
    }
  }
  return count;
}


} // namespace qmcplusplus
