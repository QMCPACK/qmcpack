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
namespace qmcplusplus
{
  DiracDeterminantSoA::DiracDeterminantSoA(SPOSetBasePtr const &spos, int first): 
    DiracDeterminantBase(spos,first)
  { }

  /** nothing useful yet */
  DiracDeterminantSoA::~DiracDeterminantSoA() {}
  DiracDeterminantSoA::DiracDeterminantSoA(const DiracDeterminantSoA& s):DiracDeterminantBase(s){}

  //void DiracDeterminantSoA::set(int first, int nel)
  //{
  //  FirstIndex = first;
  //  resize(nel,nel);
  //}

  void DiracDeterminantSoA::resize(int nel, int morb)
  {
    int norb=morb;
    if(norb <= 0) norb = nel; // for morb == -1 (default)

    //compute the sizes to ensure alignment
    NumPtcls=nel;
    NumOrbitals=norb;
    NorbPad=getAlignedSize<ValueType>(norb);
    BlockSize=NorbPad*(OHMMS_DIM+1);

    /** use aligned sized */
    Phi->t_logpsi.resize(nel,NorbPad);
    psiM.resize(nel,NorbPad);
    psiV.resize(norb);

#ifdef MIXED_PRECISION
    psiM_hp.resize(nel,NorbPad);
#endif
    LastIndex = FirstIndex + nel;

    memoryPool.resize(nel*BlockSize);
    mGL.resize(nel);
    for(size_t i=0; i<nel; ++i)
      mGL[i].resetByRef(norb,NorbPad,memoryPool.data()+i*BlockSize);
    vGL.resize(norb);
  }

  template<typename T>
  inline void computeGL(T* row, VectorSoaContainer<T,4>& gl_v, TinyVector<T,3>& grad, T& lap)
  {
    constexpr T czero(0);
    constexpr T cone(1);
    int four=4;
    int na=gl_v.size();
    int lda=gl_v.capacity();
    T y[]={czero,czero,czero,czero};
    BLAS::gemv('T',na,four,cone,gl_v.data(),lda,row,1,czero,y,1);
    grad[0]=y[0]; grad[1]=y[1]; grad[2]=y[2]; lap=y[3];
  }

  template<typename T>
  inline TinyVector<T,3> computeG(T* row, VectorSoaContainer<T,4>& gl_v)
  {
    constexpr T czero(0);
    constexpr T cone(1);
    int three=3;
    int na=gl_v.size();
    int lda=gl_v.capacity();
    T y[]={czero,czero,czero,czero};
    BLAS::gemv('T',na,three,cone,gl_v.data(),lda,row,1,czero,y,1);
    return TinyVector<T,3>(y[0],y[1],y[2]);
  }


  DiracDeterminantSoA::RealType
    DiracDeterminantSoA::evaluateLog(ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
    {
      if(NumPtcls==1)
      {
        APP_ABORT("Cannot do a single particle with DiracDeterminantSoA");
        ////CurrentDet=psiM(0,0);
        //ValueType det=psiM(0,0);
        //ValueType y=(RealType)1.0/det;
        //psiM(0,0)=y;
        //G(FirstIndex) += rv;
        //L(FirstIndex) += y*d2psiM(0,0) - dot(rv,rv);
        //LogValue = evaluateLogAndPhase(det,PhaseValue);
      }
      else
      {
        recompute(P);
        GradType rv;
        ValueType lap;
        for(size_t i=0,iat=FirstIndex; i<NumPtcls; ++i,++iat)
        {
          computeGL(psiM[i],mGL[i],rv,lap);
          P.G[iat]+=rv;
          P.L[iat]+=lap-dot(rv,rv);
        }

        RatioTimer.stop();
      }
      psiM_temp = psiM;
      return LogValue;
    }

  DiracDeterminantBase::ValueType 
    DiracDeterminantSoA::ratio(ParticleSet& P, int iat)
    {
      UpdateMode=ORB_PBYP_RATIO;
      WorkingIndex = iat-FirstIndex;
      SPOVTimer.start();
      Phi->evaluate(P, iat, psiV);
      SPOVTimer.stop();
      RatioTimer.start();
      curRatio=simd::dot(psiM[WorkingIndex],psiV.data(),NumOrbitals);
      RatioTimer.stop();
      return curRatio;
    }

  DiracDeterminantBase::GradType
    DiracDeterminantSoA::evalGrad(ParticleSet& P, int iat)
    {
      WorkingIndex = iat-FirstIndex;
      if(ReCompute)
      {//need to compute mGL
        VectorViewer<ValueType> arow(psiV.data(),NumOrbitals);
        Phi->evaluateVGL(P, iat, arow, mGL[WorkingIndex],false); 
      }
      return computeG(psiM[WorkingIndex],mGL[WorkingIndex]);
    }

  DiracDeterminantSoA::ValueType
    DiracDeterminantSoA::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
    {
      SPOVGLTimer.start();
      VectorViewer<ValueType> arow(psiV.data(),NumOrbitals);
      Phi->evaluateVGL(P, iat, arow, vGL,true); //use the new position
      SPOVGLTimer.stop();
      RatioTimer.start();
      WorkingIndex = iat-FirstIndex;
      UpdateMode=ORB_PBYP_PARTIAL;
      curRatio=simd::dot(psiM[WorkingIndex],psiV.data(),NumOrbitals);

      constexpr RealType cone(1);

      //GradType rv=simd::dot(psiM[WorkingIndex],dpsiV.data(),NumOrbitals);
      GradType rv= computeG(psiM[WorkingIndex],vGL);
      grad_iat += (cone/curRatio) * rv;
      RatioTimer.stop();
      return curRatio;
    }

  void DiracDeterminantSoA::restore(int iat)
  {
    curRatio=1.0;
  }

  /** move was accepted, update the real container
  */
  void DiracDeterminantSoA::acceptMove(ParticleSet& P, int iat)
  {
    PhaseValue += evaluatePhase(curRatio);
    LogValue +=std::log(std::abs(curRatio));
    UpdateTimer.start();

    detEng.updateRow(psiM,psiV.data(),WorkingIndex,curRatio);

    if(UpdateMode==ORB_PBYP_PARTIAL)
      simd::copy_n(vGL.data(), BlockSize, mGL[WorkingIndex].data());

    UpdateTimer.stop();
    curRatio=1.0;
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
    simd::transpose(Phi->t_logpsi.data(), NumOrbitals, NorbPad, psiM_hp.data(), NumOrbitals, NorbPad);
    detEng_hp.invert(psiM_hp,true);
    LogValue  =static_cast<RealType>(detEng_hp.LogDet);
    PhaseValue=static_cast<RealType>(detEng_hp.Phase);
    psiM=psiM_hp;
#else
    simd::transpose(Phi->t_logpsi.data(), NumOrbitals, NorbPad, psiM.data(), NumOrbitals, NorbPad);
    detEng.invert(psiM,true);
    LogValue  =detEng.LogDet;
    PhaseValue=detEng.Phase;
#endif
  }


  DiracDeterminantSoA::RealType
    DiracDeterminantSoA::registerData(ParticleSet& P, PooledData<RealType>& buf)
    {
      NP=P.getTotalNum();
      LogValue=evaluateLog(P,P.G,P.L);
      //add the data: determinant, inverse, gradient and laplacians
      buf.add(psiM.first_address(),psiM.last_address());
      if(!ReCompute) buf.add(memoryPool.begin(),memoryPool.end());
      buf.add(LogValue);
      buf.add(PhaseValue);
      return LogValue;
    }

  DiracDeterminantSoA::RealType 
    DiracDeterminantSoA::updateBuffer(ParticleSet& P, PooledData<RealType>& buf, bool fromscratch)
    {
      if(fromscratch)
      {
        LogValue=evaluateLog(P,P.G,P.L);
      }
      else
      {
        GradType grad;
        ValueType lap;
        if(UpdateMode == ORB_PBYP_RATIO) //ratio only method does not compute GL
        {
          for(size_t i=0,iat=FirstIndex; i<NumPtcls; ++i,++iat)
          {
            VectorViewer<ValueType> arow(psiV.data(),NumOrbitals);
            Phi->evaluateVGL(P, i, arow, mGL[i], true); 
            computeGL(psiM[i],mGL[i],grad,lap);
            P.G[iat]+=grad;
            P.L[iat]+=lap-dot(grad,grad);
          }
        }
        else
        {
          for(size_t i=0,iat=FirstIndex; i<NumPtcls; ++i,++iat)
          {
            computeGL(psiM[i],mGL[i],grad,lap);
            P.G[iat]+=grad;
            P.L[iat]+=lap-dot(grad,grad);
          }
        }
      }

      BufferTimer.start();
      buf.put(psiM.first_address(),psiM.last_address());
      if(!ReCompute) buf.put(memoryPool.begin(),memoryPool.end());
      buf.put(LogValue);
      buf.put(PhaseValue);
      BufferTimer.stop();
      return LogValue;
    }

  void DiracDeterminantSoA::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf)
  {
    buf.get(psiM.first_address(),psiM.last_address());
    if(!ReCompute) buf.get(memoryPool.begin(),memoryPool.end());
    buf.get(LogValue);
    buf.get(PhaseValue);
  }


}
