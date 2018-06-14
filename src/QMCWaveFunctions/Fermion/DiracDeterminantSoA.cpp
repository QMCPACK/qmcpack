//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: 
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
    
    
/**@file DiracDeterminantSoA.h
 */
#include "QMCWaveFunctions/Fermion/DiracDeterminantSoA.h"
#include "QMCWaveFunctions/Fermion/dirac_computeGL.h"

namespace qmcplusplus
{
  DiracDeterminantSoA::DiracDeterminantSoA(SPOSetBasePtr const &spos, int first): 
    DiracDeterminantBase(spos,first)
  { 
    Need2Compute4PbyP=false; 
  }

  /** nothing useful yet */
  DiracDeterminantSoA::~DiracDeterminantSoA() {}
  DiracDeterminantSoA::DiracDeterminantSoA(const DiracDeterminantSoA& s):DiracDeterminantBase(s){}

  DiracDeterminantBase* 
    DiracDeterminantSoA::makeCopy(SPOSetBasePtr spo) const
    {
      DiracDeterminantSoA* dclone= new DiracDeterminantSoA(spo,FirstIndex);
      dclone->resize(NumPtcls,NumOrbitals);
      return dclone;
    }

  //void DiracDeterminantSoA::set(int first, int nel)
  //{
  //  FirstIndex = first;
  //  resize(nel,nel);
  //}

  void DiracDeterminantSoA::resize(int nel, int morb)
  {
    size_t norb=morb;
    if(norb <= 0) norb = nel; // for morb == -1 (default)
    //compute the sizes to ensure alignment
    NumPtcls=nel;
    NumOrbitals=norb;
    NorbPad=getAlignedSize<ValueType>(norb);
    LastIndex = FirstIndex + nel;

    /** use aligned sized */
    psiM.resize(nel,NorbPad);
    psiM_temp.resize(nel,NorbPad);
#ifdef MIXED_PRECISION
    psiM_hp.resize(nel,NorbPad);
#endif
    psiV.resize(NorbPad); 

    BlockSize=NorbPad*(OHMMS_DIM+1);
    memoryPool.resize(nel*BlockSize);
    mGL.resize(nel);
    //map mGL[i] to a block of memoryPool
    for(size_t i=0; i<nel; ++i)
      mGL[i].attachReference(norb,NorbPad,memoryPool.data()+i*BlockSize);
    vVGL.resize(norb);
  }

  DiracDeterminantBase::GradType
    DiracDeterminantSoA::evalGrad(ParticleSet& P, int iat)
    {
      WorkingIndex = iat-FirstIndex;
      return computeG(psiM[WorkingIndex],mGL[WorkingIndex]);
    }

  DiracDeterminantSoA::ValueType
    DiracDeterminantSoA::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
    {
      Phi->evaluateVGL(P, iat, vVGL,true); //use the new position
      WorkingIndex = iat-FirstIndex;

      UpdateMode=ORB_PBYP_PARTIAL;
      curRatio=simd::dot(psiM[WorkingIndex],vVGL.data(0),NumOrbitals);

      ValueType f=RealType(1)/curRatio;
      GradType rv= computeG(psiM[WorkingIndex],vVGL);
      grad_iat += f * rv;
      return curRatio;
    }

  /** move was accepted, update the real container
  */
  void DiracDeterminantSoA::acceptMove(ParticleSet& P, int iat)
  {
    PhaseValue += evaluatePhase(curRatio);
    LogValue +=std::log(std::abs(curRatio));

    if(UpdateMode == ORB_PBYP_PARTIAL)
    {
      detEng.updateRow(psiM,vVGL.data(0),WorkingIndex,curRatio);
      simd::copy_n(vVGL.data(1), BlockSize, mGL[WorkingIndex].data());
    }
    else
    {
      detEng.updateRow(psiM,psiV.data(),WorkingIndex,curRatio);
    }

    curRatio=ValueType(1);
  }

  /** done PbyP update, prepare for the measurements */
  void DiracDeterminantSoA::updateAfterSweep(ParticleSet& P,
      ParticleSet::ParticleGradient_t& G,
      ParticleSet::ParticleLaplacian_t& L)
  {
    if(UpdateMode == ORB_PBYP_RATIO) 
    {//ratio only method need to compute mGL
      bool newp=false;
      for(size_t i=0; i<NumPtcls; ++i)
      {
        Phi->evaluateVGL(P, i+FirstIndex, vVGL, newp); 
        simd::copy_n(vVGL.data(1), BlockSize, mGL[i].data());
      }
    }

    mGradType grad;
    mValueType lap;
    for(size_t i=0,iat=FirstIndex; i<NumPtcls; ++i,++iat)
    {
      computeGL(psiM[i],mGL[i],grad,lap);
      G[iat]+=grad;
      L[iat]+=lap-dot(grad,grad);
    }
  }

  void
    DiracDeterminantSoA::registerData(ParticleSet& P, WFBufferType& buf)
    {
      //add the data: determinant, inverse, gradient and laplacians
      buf.add(psiM.first_address(),psiM.last_address());
      buf.add(memoryPool.data(),memoryPool.data()+memoryPool.size());
      buf.add(LogValue);
      buf.add(PhaseValue);
    }

  DiracDeterminantSoA::RealType 
    DiracDeterminantSoA::updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch)
    {
      if(fromscratch)
        LogValue=evaluateLog(P,P.G,P.L);
      else
        updateAfterSweep(P,P.G,P.L);

      buf.put(psiM.first_address(),psiM.last_address());
      buf.put(memoryPool.data(),memoryPool.data()+memoryPool.size());
      buf.put(LogValue);
      buf.put(PhaseValue);
      return LogValue;
    }

  void DiracDeterminantSoA::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
  {
    buf.get(psiM.first_address(),psiM.last_address());
    buf.get(memoryPool.data(),memoryPool.data()+memoryPool.size());
    buf.get(LogValue);
    buf.get(PhaseValue);
  }

  DiracDeterminantBase::ValueType 
    DiracDeterminantSoA::ratio(ParticleSet& P, int iat)
    {
      UpdateMode=ORB_PBYP_RATIO;
      WorkingIndex = iat-FirstIndex;
      Phi->evaluate(P, iat, psiV);
      curRatio=simd::dot(psiM[WorkingIndex],psiV.data(),NumOrbitals);
      return curRatio;
    }

  DiracDeterminantSoA::RealType
    DiracDeterminantSoA::evaluateLog(ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
    {
      recompute(P);
      mGradType rv;
      mValueType lap;
      for(size_t i=0,iat=FirstIndex; i<NumPtcls; ++i,++iat)
      {
        computeGL(psiM[i],mGL[i],rv,lap);
        G[iat]+=rv;
        L[iat]+=lap-dot(rv,rv);
      }
      return LogValue;
    }

  /** recompute the inverse
   */
  void DiracDeterminantSoA::recompute(ParticleSet& P)
  { 
    bool curpos=false;
    for(size_t i=0,iat=FirstIndex; i<NumPtcls; ++i,++iat)
    {
      Phi->evaluateVGL(P, iat, vVGL, curpos); 
      simd::copy_n(vVGL.data(0), NumOrbitals, psiM_temp[i]);
      simd::copy_n(vVGL.data(1), BlockSize,  mGL[i].data());
    }
#ifdef MIXED_PRECISION
    simd::transpose(psiM_temp.data(), NumOrbitals, NorbPad, psiM_hp.data(), NumOrbitals, psiM_hp.cols());
    detEng_hp.invert(psiM_hp,true);
    LogValue  =static_cast<RealType>(detEng_hp.LogDet);
    PhaseValue=static_cast<RealType>(detEng_hp.Phase);
    psiM=psiM_hp;
#else
    simd::transpose(psiM_temp.data(), NumOrbitals, NorbPad, psiM.data(), NumOrbitals, psiM.cols());
    detEng.invert(psiM,true);
    LogValue  =detEng.LogDet;
    PhaseValue=detEng.Phase;
#endif
  }




}
