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
    BufferMode=1;
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
    BlockSize=NorbPad*(OHMMS_DIM+1);

    /** use aligned sized */
    Phi->t_logpsi.resize(nel,NorbPad);
    psiM.resize(nel,NorbPad);
    psiV.resize(NorbPad);

#ifdef MIXED_PRECISION
    size_t nh=getAlignedSize<mValueType>(norb);
    psiM_hp.resize(nel,nh);
#endif
    LastIndex = FirstIndex + nel;

    memoryPool.resize(nel*BlockSize);
    mGL.resize(nel);
    
    if(mGL.capacity()!= NorbPad)
    {
      APP_ABORT("DiracDeterminantSoA::resize failed due to size mismatch");
    }

    for(size_t i=0; i<nel; ++i)
      mGL[i].resetByRef(norb,NorbPad,memoryPool.data()+i*BlockSize);
    vGL.resize(norb);
  }

  DiracDeterminantBase::GradType
    DiracDeterminantSoA::evalGrad(ParticleSet& P, int iat)
    {
      WorkingIndex = iat-FirstIndex;
      if(Need2Compute4PbyP)//need to compute mGL because
      {
        VectorViewer<ValueType> arow(psiV.data(),NumOrbitals);
        Phi->evaluateVGL(P, iat, arow, mGL[WorkingIndex],false); 
      }
      return computeG(psiM[WorkingIndex],mGL[WorkingIndex]);
    }

  DiracDeterminantSoA::ValueType
    DiracDeterminantSoA::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
    {
      VectorViewer<ValueType> arow(psiV.data(),NumOrbitals);
      Phi->evaluateVGL(P, iat, arow, vGL,true); //use the new position
      WorkingIndex = iat-FirstIndex;

      UpdateMode=ORB_PBYP_PARTIAL;
      curRatio=simd::dot(psiM[WorkingIndex],psiV.data(),NumOrbitals);

      ValueType f=RealType(1)/curRatio;
      GradType rv= computeG(psiM[WorkingIndex],vGL);
      grad_iat += f * rv;
      return curRatio;
    }

  /** move was accepted, update the real container
  */
  void DiracDeterminantSoA::acceptMove(ParticleSet& P, int iat)
  {
    PhaseValue += evaluatePhase(curRatio);
    LogValue +=std::log(std::abs(curRatio));

    detEng.updateRow(psiM,psiV.data(),WorkingIndex,curRatio);

    if(UpdateMode==ORB_PBYP_PARTIAL)
      simd::copy_n(vGL.data(), BlockSize, mGL[WorkingIndex].data());

    curRatio=ValueType(1);
  }

  /** done PbyP update, prepare for the measurements */
  void DiracDeterminantSoA::updateAfterSweep(ParticleSet& P,
      ParticleSet::ParticleGradient_t& G,
      ParticleSet::ParticleLaplacian_t& L)
  {
    if(UpdateMode == ORB_PBYP_RATIO) //ratio only method does not compute GL
    {
      bool newp=false;
      for(size_t i=0,iat=FirstIndex; i<NumPtcls; ++i,++iat)
      {
        VectorViewer<ValueType> arow(Phi->t_logpsi[i],NumOrbitals);
        Phi->evaluateVGL(P, iat, arow, mGL[i],newp); 
      }
    }

    GradType grad;
    ValueType lap;
    for(size_t i=0,iat=FirstIndex; i<NumPtcls; ++i,++iat)
    {
      computeGL(psiM[i],mGL[i],grad,lap);
      G[iat]+=grad;
      L[iat]+=lap-dot(grad,grad);
    }
  }

  DiracDeterminantSoA::RealType
    DiracDeterminantSoA::registerData(ParticleSet& P, PooledData<RealType>& buf)
    {
      LogValue=evaluateLog(P,P.G,P.L);
      //add the data: determinant, inverse, gradient and laplacians
      buf.add(psiM.first_address(),psiM.last_address());
      if(BufferMode)
        buf.add(memoryPool.data(),memoryPool.data()+memoryPool.size());
      buf.add(LogValue);
      buf.add(PhaseValue);
      return LogValue;
    }

  DiracDeterminantSoA::RealType
    DiracDeterminantSoA::evaluateLog(ParticleSet& P, PooledData<RealType>& buf)
    {
      buf.put(psiM.first_address(),psiM.last_address());

      if(BufferMode)
        buf.put(memoryPool.data(),memoryPool.data()+memoryPool.size());
      buf.put(LogValue);
      buf.put(PhaseValue);
      return LogValue;
    }


  DiracDeterminantSoA::RealType 
    DiracDeterminantSoA::updateBuffer(ParticleSet& P, PooledData<RealType>& buf, bool fromscratch)
    {
      if(fromscratch)
        LogValue=evaluateLog(P,P.G,P.L);
      else
        updateAfterSweep(P,P.G,P.L);

      buf.put(psiM.first_address(),psiM.last_address());
      if(BufferMode)
        buf.put(memoryPool.data(),memoryPool.data()+memoryPool.size());
      buf.put(LogValue);
      buf.put(PhaseValue);
      return LogValue;
    }

  void DiracDeterminantSoA::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf)
  {
    buf.get(psiM.first_address(),psiM.last_address());
    if(BufferMode)
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
      GradType rv;
      ValueType lap;
      for(size_t i=0,iat=FirstIndex; i<NumPtcls; ++i,++iat)
      {
        computeGL(psiM[i],mGL[i],rv,lap);
        G[iat]+=rv;
        L[iat]+=lap-dot(rv,rv);
      }
      return LogValue;
    }

  void DiracDeterminantSoA::recompute(ParticleSet& P)
  {
    bool newp=false;
    for(size_t i=0,iat=FirstIndex; i<NumPtcls; ++i,++iat)
    {
      VectorViewer<ValueType> arow(Phi->t_logpsi[i],NumOrbitals);
      Phi->evaluateVGL(P, iat, arow, mGL[i],newp); 
    }

#ifdef MIXED_PRECISION
    simd::transpose(Phi->t_logpsi.data(), NumOrbitals, NorbPad, psiM_hp.data(), NumOrbitals, psiM_hp.cols());
    detEng_hp.invert(psiM_hp,true);
    LogValue  =static_cast<RealType>(detEng_hp.LogDet);
    PhaseValue=static_cast<RealType>(detEng_hp.Phase);
    psiM=psiM_hp;
#else
    simd::transpose(Phi->t_logpsi.data(), NumOrbitals, NorbPad, psiM.data(), NumOrbitals, psiM.cols());
    detEng.invert(psiM,true);
    LogValue  =detEng.LogDet;
    PhaseValue=detEng.Phase;
#endif
  }




}
