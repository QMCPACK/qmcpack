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
    
    

#include "QMCWaveFunctions/Fermion/MultiDiracDeterminantBase.h"
//#include "QMCWaveFunctions/Fermion/MultiDiracDeterminantBase_help.h"
#include "QMCWaveFunctions/Fermion/ci_configuration2.h"
#include "Message/Communicate.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/OhmmsBlas.h"
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
void MultiDiracDeterminantBase::set(int first, int nel)
{
  APP_ABORT("MultiDiracDeterminantBase::set(int first, int nel) is disabled. \n");
}

void MultiDiracDeterminantBase::set(int first, int nel,int norb)
{
  FirstIndex = first;
  DetCalculator.resize(nel);
  resize(nel,norb);
// mmorales; remove later
//    testDets();
}

void MultiDiracDeterminantBase::createDetData(ci_configuration2& ref, std::vector<int>& data,
    std::vector<std::pair<int,int> >& pairs, std::vector<RealType>& sign)
{
  const auto& confgList=*ciConfigList;

  size_t nci = confgList.size(), nex;
  std::vector<size_t> pos(NumPtcls);
  std::vector<size_t> ocp(NumPtcls);
  std::vector<size_t> uno(NumPtcls);
  data.clear();
  sign.resize(nci);
  pairs.clear();
  for(size_t i=0; i<nci; i++)
  {
    sign[i] = ref.calculateExcitations(confgList[i],nex,pos,ocp,uno);
    data.push_back(nex);
    for(int k=0; k<nex; k++)
      data.push_back(pos[k]);
    for(int k=0; k<nex; k++)
      data.push_back(uno[k]);
    for(int k=0; k<nex; k++)
      data.push_back(ocp[k]);
    // determine unique pairs, to avoid redundant calculation of matrix elements
    // if storing the entire MOxMO matrix is too much, then make an array and a mapping to it.
    // is there an easier way??
    for(int k1=0; k1<nex; k1++)
      for(int k2=0; k2<nex; k2++)
      {
//           std::pair<int,int> temp(ocp[k1],uno[k2]);
        std::pair<int,int> temp(pos[k1],uno[k2]);
        if(find(pairs.begin(),pairs.end(),temp) == pairs.end()) //pair is new
          pairs.push_back(temp);
      }
  }
  app_log()<<"Number of terms in pairs array: " <<pairs.size() << std::endl;
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
void out1(int n, std::string str="NULL") {}
//{ std::cout <<"MDD: " <<str <<"  " <<n << std::endl; std::cout.flush(); }


void MultiDiracDeterminantBase::evaluateForWalkerMove(ParticleSet& P, bool fromScratch)
{
  evalWTimer.start();
  if(fromScratch)
    Phi->evaluate_notranspose(P,FirstIndex,LastIndex,psiM,dpsiM,d2psiM);
  if(NumPtcls==1)
  {
    //APP_ABORT("Evaluate Log with 1 particle in MultiDiracDeterminantBase is potentially dangerous. Fix later");
    std::vector<ci_configuration2>::const_iterator it(ciConfigList->begin());
    std::vector<ci_configuration2>::const_iterator last(ciConfigList->end());
    ValueVector_t::iterator det(detValues.begin());
    ValueMatrix_t::iterator lap(lapls.begin());
    GradMatrix_t::iterator grad(grads.begin());
    while(it != last)
    {
      int orb = (it++)->occup[0];
      *(det++) = psiM(0,orb);
      *(lap++) = d2psiM(0,orb);
      *(grad++) = dpsiM(0,orb);
    }
  }
  else
  {
    InverseTimer.start();

    const auto& confgList=*ciConfigList;

    //std::vector<int>::iterator it(confgList[ReferenceDeterminant].occup.begin());
    auto it(confgList[ReferenceDeterminant].occup.begin());
    for(size_t i=0; i<NumPtcls; i++)
    {
      for(size_t j=0; j<NumPtcls; j++)
        psiMinv(j,i) = psiM(j,*it);
      it++;
    }
    for(size_t i=0; i<NumPtcls; i++)
    {
      for(size_t j=0; j<NumOrbitals; j++)
        TpsiM(j,i) = psiM(i,j);
    }

    RealType phaseValueRef;
    RealType logValueRef=InvertWithLog(psiMinv.data(),NumPtcls,NumPtcls,WorkSpace.data(),Pivot.data(),phaseValueRef);
    InverseTimer.stop();
    const RealType detsign=(*DetSigns)[ReferenceDeterminant];
#if defined(QMC_COMPLEX)
    const RealType ratioMag = detsign*std::exp(logValueRef);
    //ValueType det0 = DetSigns[ReferenceDeterminant]*std::complex<OHMMS_PRECISION>(std::cos(phaseValueRef)*ratioMag,std::sin(phaseValueRef)*ratioMag);
    const ValueType det0 = std::complex<OHMMS_PRECISION>(std::cos(phaseValueRef)*ratioMag,std::sin(phaseValueRef)*ratioMag);
#else
    //ValueType det0 = DetSigns[ReferenceDeterminant]*std::exp(logValueRef)*std::cos(std::abs(phaseValueRef));
    const ValueType det0 = detsign*std::exp(logValueRef)*std::cos(std::abs(phaseValueRef));
#endif
    detValues[ReferenceDeterminant] = det0;
    BuildDotProductsAndCalculateRatios(ReferenceDeterminant,0,
        detValues,psiMinv,TpsiM,dotProducts,*detData,*uniquePairs,*DetSigns);
    for(size_t iat=0; iat<NumPtcls; iat++)
    {
      it = confgList[ReferenceDeterminant].occup.begin();
      GradType gradRatio;
      ValueType ratioLapl = 0.0;
      for(size_t i=0; i<NumPtcls; i++)
      {
        gradRatio += psiMinv(i,iat)*dpsiM(iat,*it);
        ratioLapl += psiMinv(i,iat)*d2psiM(iat,*it);
        it++;
      }
      grads(ReferenceDeterminant,iat) = det0*gradRatio;
      lapls(ReferenceDeterminant,iat) = det0*ratioLapl;
      for(size_t idim=0; idim<OHMMS_DIM; idim++)
      {
        dpsiMinv = psiMinv;
        it = confgList[ReferenceDeterminant].occup.begin();
        for(size_t i=0; i<NumPtcls; i++)
          psiV_temp[i] = dpsiM(iat,*(it++))[idim];
        InverseUpdateByColumn(dpsiMinv,psiV_temp,workV1,workV2,iat,gradRatio[idim]);
        //MultiDiracDeterminantBase::InverseUpdateByColumn_GRAD(dpsiMinv,dpsiV,workV1,workV2,iat,gradRatio[idim],idim);
        for(size_t i=0; i<NumOrbitals; i++)
          TpsiM(i,iat) = dpsiM(iat,i)[idim];
        BuildDotProductsAndCalculateRatios(ReferenceDeterminant,iat,
            grads,dpsiMinv,TpsiM,dotProducts,*detData,*uniquePairs,*DetSigns,idim);
      }
      dpsiMinv = psiMinv;
      it = confgList[ReferenceDeterminant].occup.begin();
      for(size_t i=0; i<NumPtcls; i++)
        psiV_temp[i] = d2psiM(iat,*(it++));
      InverseUpdateByColumn(dpsiMinv,psiV_temp,workV1,workV2,iat,ratioLapl);
      //MultiDiracDeterminantBase::InverseUpdateByColumn(dpsiMinv,d2psiM,workV1,workV2,iat,ratioLapl,confgList[ReferenceDeterminant].occup.begin());
      for(size_t i=0; i<NumOrbitals; i++)
        TpsiM(i,iat) = d2psiM(iat,i);
      BuildDotProductsAndCalculateRatios(ReferenceDeterminant,iat,
          lapls,dpsiMinv,TpsiM,dotProducts,*detData,*uniquePairs,*DetSigns);
// restore matrix
      for(size_t i=0; i<NumOrbitals; i++)
        TpsiM(i,iat) = psiM(iat,i);
    }
  } // NumPtcls==1
  psiMinv_temp = psiMinv;
  evalWTimer.stop();
}


MultiDiracDeterminantBase::RealType MultiDiracDeterminantBase::updateBuffer(ParticleSet& P,
    WFBufferType& buf, bool fromscratch)
{
  evaluateForWalkerMove(P,(fromscratch || UpdateMode == ORB_PBYP_RATIO) );
  buf.put(psiM.first_address(),psiM.last_address());
  buf.put(FirstAddressOfdpsiM,LastAddressOfdpsiM);
  buf.put(d2psiM.first_address(),d2psiM.last_address());
  buf.put(psiMinv.first_address(),psiMinv.last_address());
  buf.put(detValues.first_address(), detValues.last_address());
  buf.put(FirstAddressOfGrads,LastAddressOfGrads);
  buf.put(lapls.first_address(), lapls.last_address());
  return 1.0;
}

void MultiDiracDeterminantBase::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  buf.get(psiM.first_address(),psiM.last_address());
  buf.get(FirstAddressOfdpsiM,LastAddressOfdpsiM);
  buf.get(d2psiM.first_address(),d2psiM.last_address());
  buf.get(psiMinv.first_address(),psiMinv.last_address());
  buf.get(detValues.first_address(), detValues.last_address());
  buf.get(FirstAddressOfGrads,LastAddressOfGrads);
  buf.get(lapls.first_address(), lapls.last_address());
  // only used with ORB_PBYP_ALL,
  //psiM_temp = psiM;
  //dpsiM_temp = dpsiM;
  //d2psiM_temp = d2psiM;
  psiMinv_temp = psiMinv;
  int n1 = psiM.extent(0);
  int n2 = psiM.extent(1);
  for(int i=0; i<n1; i++)
    for(int j=0; j<n2; j++)
      TpsiM(j,i) = psiM(i,j);
}

/** move was accepted, update the real container
*/
void MultiDiracDeterminantBase::acceptMove(ParticleSet& P, int iat)
{
  WorkingIndex = iat-FirstIndex;
  switch(UpdateMode)
  {
  case ORB_PBYP_RATIO:
    psiMinv = psiMinv_temp;
    for(int i=0; i<NumOrbitals; i++)
      TpsiM(i,WorkingIndex) = psiV[i];
    std::copy(psiV.begin(),psiV.end(),psiM[iat-FirstIndex]);
    std::copy(new_detValues.begin(),new_detValues.end(),detValues.begin());
    break;
  case ORB_PBYP_PARTIAL:
    psiMinv = psiMinv_temp;
    for(int i=0; i<NumOrbitals; i++)
      TpsiM(i,WorkingIndex) = psiV[i];
    std::copy(new_detValues.begin(),new_detValues.end(),detValues.begin());
// no use saving these
//        for(int i=0; i<NumDets; i++)
//          grads(i,WorkingIndex) = new_grads(i,WorkingIndex);
    std::copy(psiV.begin(),psiV.end(),psiM[WorkingIndex]);
    std::copy(dpsiV.begin(),dpsiV.end(),dpsiM[WorkingIndex]);
    std::copy(d2psiV.begin(),d2psiV.end(),d2psiM[WorkingIndex]);
    break;
  default:
    psiMinv = psiMinv_temp;
    for(int i=0; i<NumOrbitals; i++)
      TpsiM(i,WorkingIndex) = psiV[i];
    std::copy(new_detValues.begin(),new_detValues.end(),detValues.begin());
    std::copy(new_grads.begin(),new_grads.end(),grads.begin());
    std::copy(new_lapls.begin(),new_lapls.end(),lapls.begin());
    std::copy(psiV.begin(),psiV.end(),psiM[WorkingIndex]);
    std::copy(dpsiV.begin(),dpsiV.end(),dpsiM[WorkingIndex]);
    std::copy(d2psiV.begin(),d2psiV.end(),d2psiM[WorkingIndex]);
    break;
  }
}

/** move was rejected. copy the real container to the temporary to move on
*/
void MultiDiracDeterminantBase::restore(int iat)
{
  WorkingIndex = iat-FirstIndex;
  psiMinv_temp = psiMinv;
  for(int i=0; i<NumOrbitals; i++)
    TpsiM(i,WorkingIndex) = psiM(WorkingIndex,i);
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
MultiDiracDeterminantBase::MultiDiracDeterminantBase(const MultiDiracDeterminantBase& s):
  OrbitalBase(s), NP(0), FirstIndex(s.FirstIndex),ciConfigList(nullptr),
  UpdateTimer("MultiDiracDeterminantBase::update"),
  RatioTimer("MultiDiracDeterminantBase::ratio"),
  InverseTimer("MultiDiracDeterminantBase::inverse"),
  buildTableTimer("MultiDiracDeterminantBase::buildTable"),
  evalWTimer("MultiDiracDeterminantBase::evalW"),
  evalOrbTimer("MultiDiracDeterminantBase::evalOrb"),
  evalOrb1Timer("MultiDiracDeterminantBase::evalOrbGrad"),
  readMatTimer("MultiDiracDeterminantBase::readMat"),
  readMatGradTimer("MultiDiracDeterminantBase::readMatGrad"),
  buildTableGradTimer("MultiDiracDeterminantBase::buildTableGrad"),
  ExtraStuffTimer("MultiDiracDeterminantBase::ExtraStuff")
{
  IsCloned=true;

  ReferenceDeterminant = s.ReferenceDeterminant;
  ciConfigList=s.ciConfigList;
  NumDets = s.NumDets;
  detData = s.detData;
  uniquePairs = s.uniquePairs;
  DetSigns = s.DetSigns;

  registerTimers();
  Phi = (s.Phi->makeClone());
  this->resize(s.NumPtcls,s.NumOrbitals);
  this->DetCalculator.resize(s.NumPtcls);
}

SPOSetBasePtr  MultiDiracDeterminantBase::clonePhi() const
{
  return Phi->makeClone();
}

OrbitalBasePtr MultiDiracDeterminantBase::makeClone(ParticleSet& tqp) const
{
  APP_ABORT(" Illegal action. Cannot use MultiDiracDeterminantBase::makeClone");
  return 0;
}

/** constructor
 *@param spos the single-particle orbital set
 *@param first index of the first particle
 */
MultiDiracDeterminantBase::MultiDiracDeterminantBase(SPOSetBasePtr const &spos, int first):
  NP(0),Phi(spos),FirstIndex(first),ReferenceDeterminant(0), ciConfigList(nullptr),
  UpdateTimer("MultiDiracDeterminantBase::update"),
  RatioTimer("MultiDiracDeterminantBase::ratio"),
  InverseTimer("MultiDiracDeterminantBase::inverse"),
  buildTableTimer("MultiDiracDeterminantBase::buildTable"),
  evalWTimer("MultiDiracDeterminantBase::evalW"),
  evalOrbTimer("MultiDiracDeterminantBase::evalOrb"),
  evalOrb1Timer("MultiDiracDeterminantBase::evalOrbGrad"),
  readMatTimer("MultiDiracDeterminantBase::readMat"),
  readMatGradTimer("MultiDiracDeterminantBase::readMatGrad"),
  buildTableGradTimer("MultiDiracDeterminantBase::buildTableGrad"),
  ExtraStuffTimer("MultiDiracDeterminantBase::ExtraStuff")
{
  Optimizable=true;
  OrbitalName="MultiDiracDeterminantBase";

  IsCloned=false;

  ciConfigList=new std::vector<ci_configuration2>;
  detData=new std::vector<int>;
  uniquePairs=new std::vector<std::pair<int,int> >;
  DetSigns=new std::vector<RealType>;

  registerTimers();
}

///default destructor
MultiDiracDeterminantBase::~MultiDiracDeterminantBase() {}

MultiDiracDeterminantBase& MultiDiracDeterminantBase::operator=(const MultiDiracDeterminantBase& s)
{
  if(this == & s) return *this;

  NP=0;
  IsCloned=true;
  ReferenceDeterminant = s.ReferenceDeterminant;
  ciConfigList=s.ciConfigList;
  NumDets = s.NumDets;

  detData=s.detData;
  uniquePairs=s.uniquePairs;
  FirstIndex=s.FirstIndex;
  DetSigns=s.DetSigns;

  resize(s.NumPtcls, s.NumOrbitals);
  this->DetCalculator.resize(s.NumPtcls);

  return *this;
}

void
MultiDiracDeterminantBase::registerData(ParticleSet& P, WFBufferType& buf)
{
  if(NP == 0)
    //first time, allocate once
  {
    //int norb = cols();
    NP=P.getTotalNum();
    FirstAddressOfGrads = &(grads(0,0)[0]);
    LastAddressOfGrads = FirstAddressOfGrads + NumPtcls*DIM*NumDets;
    FirstAddressOfdpsiM = &(dpsiM(0,0)[0]); //(*dpsiM.begin())[0]);
    LastAddressOfdpsiM = FirstAddressOfdpsiM + NumPtcls*NumOrbitals*DIM;
  }
  evaluateForWalkerMove(P,true);
  //add the data:
  buf.add(psiM.first_address(),psiM.last_address());
  buf.add(FirstAddressOfdpsiM,LastAddressOfdpsiM);
  buf.add(d2psiM.first_address(),d2psiM.last_address());
  buf.add(psiMinv.first_address(),psiMinv.last_address());
  buf.add(detValues.first_address(), detValues.last_address());
  buf.add(FirstAddressOfGrads,LastAddressOfGrads);
  buf.add(lapls.first_address(), lapls.last_address());
}


void MultiDiracDeterminantBase::setDetInfo(int ref, std::vector<ci_configuration2>* list)
{
  ReferenceDeterminant = ref;
  ciConfigList=list;
  NumDets = list->size();
}

///reset the size: with the number of particles and number of orbtials
/// morb is the total number of orbitals, including virtual
void MultiDiracDeterminantBase::resize(int nel, int morb)
{
  if(nel <= 0 || morb <= 0)
  {
    APP_ABORT(" ERROR: MultiDiracDeterminantBase::resize arguments equal to zero. \n");
  }
  if(NumDets == 0 || NumDets != ciConfigList->size())
  {
    APP_ABORT(" ERROR: MultiDiracDeterminantBase::resize problems with NumDets. \n");
  }

  NumPtcls=nel;
  NumOrbitals=morb;
  LastIndex = FirstIndex + nel;
  psiV_temp.resize(nel);
  psiV.resize(NumOrbitals);
  dpsiV.resize(NumOrbitals);
  d2psiV.resize(NumOrbitals);
  psiM.resize(nel,morb);
  dpsiM.resize(nel,morb);
  d2psiM.resize(nel,morb);
  //psiM_temp.resize(nel,morb);
  //dpsiM_temp.resize(nel,morb);
  //d2psiM_temp.resize(nel,morb);
  TpsiM.resize(morb,nel);
  psiMinv.resize(nel,nel);
  dpsiMinv.resize(nel,nel);
  psiMinv_temp.resize(nel,nel);
  //scratch spaces: stateless
  WorkSpace.resize(std::max(nel,NumDets)); 
  Pivot.resize(nel);
  workV1.resize(nel);
  workV2.resize(nel);
  detValues.resize(NumDets);
  new_detValues.resize(NumDets);
  grads.resize(NumDets,nel);
  new_grads.resize(NumDets,nel);
  lapls.resize(NumDets,nel);
  new_lapls.resize(NumDets,nel);
  dotProducts.resize(morb,morb);

  //if(ciConfigList==nullptr)
  //{
  //  APP_ABORT("ciConfigList was not properly initialized.\n");
  //}
  if(!IsCloned)
    createDetData((*ciConfigList)[ReferenceDeterminant], *detData,*uniquePairs,*DetSigns);
}

void MultiDiracDeterminantBase::registerTimers()
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
  TimerManager.addTimer (&UpdateTimer);
  TimerManager.addTimer (&RatioTimer);
  TimerManager.addTimer (&InverseTimer);
  TimerManager.addTimer (&buildTableTimer);
  TimerManager.addTimer (&readMatTimer);
  TimerManager.addTimer (&evalWTimer);
  TimerManager.addTimer (&evalOrbTimer);
  TimerManager.addTimer (&evalOrb1Timer);
  TimerManager.addTimer (&buildTableGradTimer);
  TimerManager.addTimer (&readMatGradTimer);
  TimerManager.addTimer (&ExtraStuffTimer);

}




}
