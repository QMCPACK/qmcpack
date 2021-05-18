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
/** set the index of the first particle in the determinant and reset the size of the determinant
 *@param first index of first particle
 *@param nel number of particles in the determinant
 */
void MultiDiracDeterminant::set(int first, int nel)
{
  APP_ABORT(ClassName + "set(int first, int nel) is disabled. \n");
}

void MultiDiracDeterminant::set(int first, int nel, int norb)
{
  FirstIndex = first;
  DetCalculator.resize(nel);
  resize(nel, norb);
  // mmorales; remove later
  //    testDets();
}

void MultiDiracDeterminant::createDetData(ci_configuration2& ref,
                                          std::vector<int>& data,
                                          std::vector<std::pair<int, int>>& pairs,
                                          std::vector<RealType>& sign)
{
  const auto& confgList = *ciConfigList;

  size_t nci = confgList.size(), nex;
  std::vector<size_t> pos(NumPtcls);
  std::vector<size_t> ocp(NumPtcls);
  std::vector<size_t> uno(NumPtcls);
  data.clear();
  sign.resize(nci);
  pairs.clear();
  for (size_t i = 0; i < nci; i++)
  {
    sign[i] = ref.calculateExcitations(confgList[i], nex, pos, ocp, uno);
    data.push_back(nex);
    for (int k = 0; k < nex; k++)
      data.push_back(pos[k]);
    for (int k = 0; k < nex; k++)
      data.push_back(uno[k]);
    for (int k = 0; k < nex; k++)
      data.push_back(ocp[k]);
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
  /*
       std::cout <<"ref: " <<ref << std::endl;
       std::cout <<"list: " << std::endl;
       for(int i=0; i<confgList.size(); i++)
         std::cout <<confgList[i] << std::endl;

       std::cout <<"pairs: " << std::endl;
       for(int i=0; i<pairs.size(); i++)
         std::cout <<pairs[i].first <<"   " <<pairs[i].second << std::endl;
  */
}

//erase
void out1(int n, std::string str = "NULL") {}
//{ std::cout <<"MDD: " <<str <<"  " <<n << std::endl; std::cout.flush(); }


void MultiDiracDeterminant::evaluateForWalkerMove(const ParticleSet& P, bool fromScratch)
{
  evalWTimer.start();
  if (fromScratch)
  {
    Phi->evaluate_notranspose_spin(P, FirstIndex, LastIndex, psiM, dpsiM, d2psiM, dspin_psiM);
  }
  if (NumPtcls == 1)
  {
    //APP_ABORT("Evaluate Log with 1 particle in MultiDiracDeterminant is potentially dangerous. Fix later");
    std::vector<ci_configuration2>::const_iterator it(ciConfigList->begin());
    std::vector<ci_configuration2>::const_iterator last(ciConfigList->end());
    ValueVector_t::iterator det(detValues.begin());
    ValueMatrix_t::iterator lap(lapls.begin());
    GradMatrix_t::iterator grad(grads.begin());
    ValueMatrix_t::iterator spingrad(spingrads.begin());
    while (it != last)
    {
      int orb       = (it++)->occup[0];
      *(det++)      = psiM(0, orb);
      *(lap++)      = d2psiM(0, orb);
      *(grad++)     = dpsiM(0, orb);
      *(spingrad++) = dspin_psiM(0, orb);
    }
  }
  else
  {
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
    InverseTimer.stop();
    const RealType detsign          = (*DetSigns)[ReferenceDeterminant];
    const ValueType det0            = LogToValue<ValueType>::convert(logValueRef);
    detValues[ReferenceDeterminant] = det0;
    BuildDotProductsAndCalculateRatios(ReferenceDeterminant, 0, detValues, psiMinv, TpsiM, dotProducts, *detData,
                                       *uniquePairs, *DetSigns);
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
      grads(ReferenceDeterminant, iat)     = det0 * gradRatio;
      lapls(ReferenceDeterminant, iat)     = det0 * ratioLapl;
      spingrads(ReferenceDeterminant, iat) = det0 * spingradRatio;
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
        BuildDotProductsAndCalculateRatios(ReferenceDeterminant, iat, grads, dpsiMinv, TpsiM, dotProducts, *detData,
                                           *uniquePairs, *DetSigns, idim);
      }
      dpsiMinv = psiMinv;
      it       = confgList[ReferenceDeterminant].occup.begin();
      for (size_t i = 0; i < NumPtcls; i++)
        psiV_temp[i] = d2psiM(iat, *(it++));
      InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, iat, ratioLapl);
      //MultiDiracDeterminant::InverseUpdateByColumn(dpsiMinv,d2psiM,workV1,workV2,iat,ratioLapl,confgList[ReferenceDeterminant].occup.begin());
      for (size_t i = 0; i < NumOrbitals; i++)
        TpsiM(i, iat) = d2psiM(iat, i);
      BuildDotProductsAndCalculateRatios(ReferenceDeterminant, iat, lapls, dpsiMinv, TpsiM, dotProducts, *detData,
                                         *uniquePairs, *DetSigns);

      dpsiMinv = psiMinv;
      it       = confgList[ReferenceDeterminant].occup.begin();
      for (size_t i = 0; i < NumPtcls; i++)
        psiV_temp[i] = dspin_psiM(iat, *(it++));
      InverseUpdateByColumn(dpsiMinv, psiV_temp, workV1, workV2, iat, spingradRatio);
      for (size_t i = 0; i < NumOrbitals; i++)
        TpsiM(i, iat) = dspin_psiM(iat, i);
      BuildDotProductsAndCalculateRatios(ReferenceDeterminant, iat, spingrads, dpsiMinv, TpsiM, dotProducts, *detData,
                                         *uniquePairs, *DetSigns);
      // restore matrix
      for (size_t i = 0; i < NumOrbitals; i++)
        TpsiM(i, iat) = psiM(iat, i);
    }
  } // NumPtcls==1
  psiMinv_temp = psiMinv;
  evalWTimer.stop();
}


MultiDiracDeterminant::LogValueType MultiDiracDeterminant::updateBuffer(ParticleSet& P,
                                                                        WFBufferType& buf,
                                                                        bool fromscratch)
{
  evaluateForWalkerMove(P, (fromscratch || UpdateMode == ORB_PBYP_RATIO));
  buf.put(psiM.first_address(), psiM.last_address());
  buf.put(FirstAddressOfdpsiM, LastAddressOfdpsiM);
  buf.put(d2psiM.first_address(), d2psiM.last_address());
  buf.put(psiMinv.first_address(), psiMinv.last_address());
  buf.put(detValues.first_address(), detValues.last_address());
  buf.put(FirstAddressOfGrads, LastAddressOfGrads);
  buf.put(lapls.first_address(), lapls.last_address());
  buf.put(dspin_psiM.first_address(), dspin_psiM.last_address());
  buf.put(spingrads.first_address(), spingrads.last_address());
  return 1.0;
}

void MultiDiracDeterminant::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  buf.get(psiM.first_address(), psiM.last_address());
  buf.get(FirstAddressOfdpsiM, LastAddressOfdpsiM);
  buf.get(d2psiM.first_address(), d2psiM.last_address());
  buf.get(psiMinv.first_address(), psiMinv.last_address());
  buf.get(detValues.first_address(), detValues.last_address());
  buf.get(FirstAddressOfGrads, LastAddressOfGrads);
  buf.get(lapls.first_address(), lapls.last_address());
  buf.get(dspin_psiM.first_address(), dspin_psiM.last_address());
  buf.get(spingrads.first_address(), spingrads.last_address());
  // only used with ORB_PBYP_ALL,
  //psiM_temp = psiM;
  //dpsiM_temp = dpsiM;
  //d2psiM_temp = d2psiM;
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
  switch (UpdateMode)
  {
  case ORB_PBYP_RATIO:
    psiMinv = psiMinv_temp;
    for (int i = 0; i < NumOrbitals; i++)
      TpsiM(i, WorkingIndex) = psiV[i];
    std::copy(psiV.begin(), psiV.end(), psiM[iat - FirstIndex]);
    std::copy(new_detValues.begin(), new_detValues.end(), detValues.begin());
    break;
  case ORB_PBYP_PARTIAL:
    psiMinv = psiMinv_temp;
    for (int i = 0; i < NumOrbitals; i++)
      TpsiM(i, WorkingIndex) = psiV[i];
    std::copy(new_detValues.begin(), new_detValues.end(), detValues.begin());
    // no use saving these
    //        for(int i=0; i<NumDets; i++)
    //          grads(i,WorkingIndex) = new_grads(i,WorkingIndex);
    std::copy(psiV.begin(), psiV.end(), psiM[WorkingIndex]);
    std::copy(dpsiV.begin(), dpsiV.end(), dpsiM[WorkingIndex]);
    std::copy(d2psiV.begin(), d2psiV.end(), d2psiM[WorkingIndex]);
    std::copy(dspin_psiV.begin(), dspin_psiV.end(), dspin_psiM[WorkingIndex]);
    break;
  default:
    psiMinv = psiMinv_temp;
    for (int i = 0; i < NumOrbitals; i++)
      TpsiM(i, WorkingIndex) = psiV[i];
    std::copy(new_detValues.begin(), new_detValues.end(), detValues.begin());
    std::copy(new_grads.begin(), new_grads.end(), grads.begin());
    std::copy(new_lapls.begin(), new_lapls.end(), lapls.begin());
    std::copy(new_spingrads.begin(), new_spingrads.end(), spingrads.begin());
    std::copy(psiV.begin(), psiV.end(), psiM[WorkingIndex]);
    std::copy(dpsiV.begin(), dpsiV.end(), dpsiM[WorkingIndex]);
    std::copy(d2psiV.begin(), d2psiV.end(), d2psiM[WorkingIndex]);
    std::copy(dspin_psiV.begin(), dspin_psiV.end(), dspin_psiM[WorkingIndex]);
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
      InverseTimer(*timer_manager.createTimer(ClassName + "::inverse")),
      buildTableTimer(*timer_manager.createTimer(ClassName + "::buildTable")),
      readMatTimer(*timer_manager.createTimer(ClassName + "::readMat")),
      evalWTimer(*timer_manager.createTimer(ClassName + "::evalW")),
      evalOrbTimer(*timer_manager.createTimer(ClassName + "::evalOrb")),
      evalOrb1Timer(*timer_manager.createTimer(ClassName + "::evalOrbGrad")),
      readMatGradTimer(*timer_manager.createTimer(ClassName + "::readMatGrad")),
      buildTableGradTimer(*timer_manager.createTimer(ClassName + "::buildTableGrad")),
      ExtraStuffTimer(*timer_manager.createTimer(ClassName + "::ExtraStuff")),
      NP(0),
      FirstIndex(s.FirstIndex),
      ciConfigList(nullptr)
{
  IsCloned = true;

  ReferenceDeterminant = s.ReferenceDeterminant;
  ciConfigList         = s.ciConfigList;
  NumDets              = s.NumDets;
  detData              = s.detData;
  uniquePairs          = s.uniquePairs;
  DetSigns             = s.DetSigns;
  Optimizable          = s.Optimizable;

  registerTimers();
  Phi.reset(s.Phi->makeClone());
  this->resize(s.NumPtcls, s.NumOrbitals);
  this->DetCalculator.resize(s.NumPtcls);
}

SPOSetPtr MultiDiracDeterminant::clonePhi() const { return Phi->makeClone(); }

WaveFunctionComponentPtr MultiDiracDeterminant::makeClone(ParticleSet& tqp) const
{
  APP_ABORT(" Illegal action. Cannot use MultiDiracDeterminant::makeClone");
  return 0;
}

/** constructor
 *@param spos the single-particle orbital set
 *@param first index of the first particle
 */
MultiDiracDeterminant::MultiDiracDeterminant(std::unique_ptr<SPOSet>&& spos, int first)
    : WaveFunctionComponent("MultiDiracDeterminant"),
      UpdateTimer(*timer_manager.createTimer(ClassName + "::update")),
      RatioTimer(*timer_manager.createTimer(ClassName + "::ratio")),
      InverseTimer(*timer_manager.createTimer(ClassName + "::inverse")),
      buildTableTimer(*timer_manager.createTimer(ClassName + "::buildTable")),
      readMatTimer(*timer_manager.createTimer(ClassName + "::readMat")),
      evalWTimer(*timer_manager.createTimer(ClassName + "::evalW")),
      evalOrbTimer(*timer_manager.createTimer(ClassName + "::evalOrb")),
      evalOrb1Timer(*timer_manager.createTimer(ClassName + "::evalOrbGrad")),
      readMatGradTimer(*timer_manager.createTimer(ClassName + "::readMatGrad")),
      buildTableGradTimer(*timer_manager.createTimer(ClassName + "::buildTableGrad")),
      ExtraStuffTimer(*timer_manager.createTimer(ClassName + "::ExtraStuff")),
      NP(0),
      FirstIndex(first),
      Phi(std::move(spos)),
      ciConfigList(nullptr),
      ReferenceDeterminant(0)
{
  (Phi->isOptimizable() == true) ? Optimizable = true : Optimizable = false;

  IsCloned = false;

  ciConfigList = new std::vector<ci_configuration2>;
  detData      = new std::vector<int>;
  uniquePairs  = new std::vector<std::pair<int, int>>;
  DetSigns     = new std::vector<RealType>;

  registerTimers();
}

///default destructor
MultiDiracDeterminant::~MultiDiracDeterminant() {}

MultiDiracDeterminant& MultiDiracDeterminant::operator=(const MultiDiracDeterminant& s)
{
  if (this == &s)
    return *this;

  NP                   = 0;
  IsCloned             = true;
  ReferenceDeterminant = s.ReferenceDeterminant;
  ciConfigList         = s.ciConfigList;
  NumDets              = s.NumDets;

  detData     = s.detData;
  uniquePairs = s.uniquePairs;
  FirstIndex  = s.FirstIndex;
  DetSigns    = s.DetSigns;

  resize(s.NumPtcls, s.NumOrbitals);
  this->DetCalculator.resize(s.NumPtcls);

  return *this;
}

void MultiDiracDeterminant::registerData(ParticleSet& P, WFBufferType& buf)
{
  if (NP == 0)
  //first time, allocate once
  {
    //int norb = cols();
    NP                  = P.getTotalNum();
    FirstAddressOfGrads = &(grads(0, 0)[0]);
    LastAddressOfGrads  = FirstAddressOfGrads + NumPtcls * DIM * NumDets;
    FirstAddressOfdpsiM = &(dpsiM(0, 0)[0]); //(*dpsiM.begin())[0]);
    LastAddressOfdpsiM  = FirstAddressOfdpsiM + NumPtcls * NumOrbitals * DIM;
  }
  evaluateForWalkerMove(P, true);
  //add the data:
  buf.add(psiM.first_address(), psiM.last_address());
  buf.add(FirstAddressOfdpsiM, LastAddressOfdpsiM);
  buf.add(d2psiM.first_address(), d2psiM.last_address());
  buf.add(psiMinv.first_address(), psiMinv.last_address());
  buf.add(detValues.first_address(), detValues.last_address());
  buf.add(FirstAddressOfGrads, LastAddressOfGrads);
  buf.add(lapls.first_address(), lapls.last_address());
  buf.add(dspin_psiM.first_address(), dspin_psiM.last_address());
  buf.add(spingrads.first_address(), spingrads.last_address());
}


void MultiDiracDeterminant::setDetInfo(int ref, std::vector<ci_configuration2>* list)
{
  ReferenceDeterminant = ref;
  ciConfigList         = list;
  NumDets              = list->size();
}

///reset the size: with the number of particles and number of orbtials
/// morb is the total number of orbitals, including virtual
void MultiDiracDeterminant::resize(int nel, int morb)
{
  if (nel <= 0 || morb <= 0)
  {
    APP_ABORT(" ERROR: MultiDiracDeterminant::resize arguments equal to zero. \n");
  }
  if (NumDets == 0 || NumDets != ciConfigList->size())
  {
    APP_ABORT(" ERROR: MultiDiracDeterminant::resize problems with NumDets. \n");
  }

  NumPtcls    = nel;
  NumOrbitals = morb;
  LastIndex   = FirstIndex + nel;
  psiV_temp.resize(nel);
  psiV.resize(NumOrbitals);
  dpsiV.resize(NumOrbitals);
  d2psiV.resize(NumOrbitals);
  dspin_psiV.resize(NumOrbitals);
  psiM.resize(nel, morb);
  dpsiM.resize(nel, morb);
  d2psiM.resize(nel, morb);
  dspin_psiM.resize(nel, morb);
  //psiM_temp.resize(nel,morb);
  //dpsiM_temp.resize(nel,morb);
  //d2psiM_temp.resize(nel,morb);
  TpsiM.resize(morb, nel);
  psiMinv.resize(nel, nel);
  dpsiMinv.resize(nel, nel);
  psiMinv_temp.resize(nel, nel);
  //scratch spaces: stateless
  WorkSpace.resize(std::max(nel, NumDets));
  Pivot.resize(nel);
  workV1.resize(nel);
  workV2.resize(nel);
  detValues.resize(NumDets);
  new_detValues.resize(NumDets);
  grads.resize(NumDets, nel);
  new_grads.resize(NumDets, nel);
  lapls.resize(NumDets, nel);
  new_lapls.resize(NumDets, nel);
  spingrads.resize(NumDets, nel);
  new_spingrads.resize(NumDets, nel);
  dotProducts.resize(morb, morb);

  //if(ciConfigList==nullptr)
  //{
  //  APP_ABORT("ciConfigList was not properly initialized.\n");
  //}
  if (!IsCloned)
    createDetData((*ciConfigList)[ReferenceDeterminant], *detData, *uniquePairs, *DetSigns);
}

void MultiDiracDeterminant::registerTimers()
{
  UpdateTimer.reset();
  RatioTimer.reset();
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
