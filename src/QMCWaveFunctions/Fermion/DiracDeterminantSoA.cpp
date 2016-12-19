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
  DiracDeterminantSoA(SPOSetBasePtr const &spos, int first): DiracDeterminantBase(spo,first)
  { }

  /** nothing useful yet */
  DiracDeterminantSoA::~DiracDeterminantSoA() {}
  DiracDeterminantSoA(const DiracDeterminantSoA& s){}
  DiracDeterminantSoA& DiracDeterminantSoA::operator=(const DiracDeterminantSoA& s) {}

  void DiracDeterminantSoA::set(int first, int nel)
  {
    FirstIndex = first;
    resize(nel,nel);
  }

  void DiracDeterminantBase::resize(int nel, int morb)
  {
    int norb=morb;
    if(norb <= 0) norb = nel; // for morb == -1 (default)
    psiM.resize(nel,norb);
    psiV.resize(norb);
#ifdef MIXED_PRECISION
    psiM_hp.resize(nel,norb);
    WorkSpace_hp.resize(nel);
#endif
    WorkSpace.resize(nel);
    Pivot.resize(nel);
    LastIndex = FirstIndex + nel;
    NumPtcls=nel;
    NumOrbitals=norb;

    /** allocate GL combo */
    NorbPad=getAlignedSize<ValueType>(norb);
    BlockSize=NorbPad*(OHMMS_DIM+1);
    memoryPool.resize(nel*BlockSize);
    mGL.resize(nel);
    for(size_t i=0; i<nel; ++i)
      mGL[i].resetByRef(Norb, NorbPad,memoryPool.data()+i*BlockSize);
    vGL.resize(Norb);
  }

  template<typename T>
  inline void computeGL(const VectorViwer<T>& row, const VectorSoaContainer<T>& gl_v, TinyVector<T,3>& grad, T& lap)
  {
    constexpr char trans='T';
    constexpr T czero(0);
    constexpr T cone(1);
    int four=4;
    int na=gl_v.size();
    int lda=gl_v.capacity();
    T y[]={czero,czero,czero,czero};
    BLAS::gemv(transa,na,four,cone,gl_v.data(),lda,row.data(),1,czero,y,1);
    grad[0]=y[0]; grad[1]=y[1]; grad[2]=y[2]; lap=y[3];
  }

  template<typename T>
  inline TinyVector<T,3> computeG(const VectorViwer<T>& row, const VectorSoaContainer<T>& gl_v)
  {
    TinyVector<T,3> grad;
    constexpr char trans='T';
    constexpr T czero(0);
    constexpr T cone(1);
    int three=3;
    int na=gl_v.size();
    int lda=gl_v.capacity();
    BLAS::gemv(transa,na,three,cone,gl_v.data(),lda,row.data(),1,czero,grad.data(),1);
    return grad;
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
        SPOVGLTimer.start();
        for(size_t i=0,iat=FirstIndex; i<NumPtcls; ++i,++iat)
          Phi->evaluate(P, iat, VectorViewer<ValueType>(psiM[i]), mGL[i]); 
        SPOVGLTimer.stop();
        InverseTimer.start();
#ifdef MIXED_PRECISION
        psiM_hp = psiM;
        ParticleSet::Scalar_t PhaseValue_hp;
        LogValue = InvertWithLog(psiM_hp.data(),NumPtcls,NumOrbitals,WorkSpace_hp.data(),Pivot.data(),PhaseValue_hp);
        psiM = psiM_hp;
        PhaseValue = PhaseValue_hp;
#else
        LogValue = InvertWithLog(psiM.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);
#endif
        InverseTimer.stop();
        RatioTimer.start();

        mGradType rv;
        mValueType lap;
        for(size_t i=0,iat=FirstIndex; i<NumPtcls; ++i,++iat)
        {
          computeGL(VectorView<ValueType>(psiM[i]),mGL[i],rv,lap);
          P.G[iat]+=rv;
          P.L[iat]+=lap-dot(rv,rv);
        }

        RatioTimer.stop();
      }
      psiM_temp = psiM;
      return LogValue;
    }

  DiracDeterminantBase::GradType
    DiracDeterminantBase::evalGrad(ParticleSet& P, int iat)
    {
      WorkingIndex = iat-FirstIndex;
      return computeG(VectorView<ValueType>(psiM[WorkingIndex]),mGL[WorkingIndex]);
    }

  DiracDeterminantSoA::ValueType
    DiracDeterminantSoA::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
    {
      SPOVGLTimer.start();
      Phi->evaluate(P, iat, VectorViewer<ValueType>(psiV.data()), vGL);
      SPOVGLTimer.stop();
      RatioTimer.start();
      WorkingIndex = iat-FirstIndex;
      UpdateMode=ORB_PBYP_PARTIAL;
      curRatio=simd::dot(psiM[WorkingIndex],psiV.data(),NumOrbitals);

      constexpr RealType cone(1);

      //GradType rv=simd::dot(psiM[WorkingIndex],dpsiV.data(),NumOrbitals);
      GradType rv= computeG(VectorView<ValueType>(psiM[WorkingIndex]),vGL);
      grad_iat += (cone/curRatio) * rv;
      RatioTimer.stop();
      return curRatio;
    }

  void DiracDeterminantBase::restore(int iat)
  {
    curRatio=1.0;
  }

  /** move was accepted, update the real container
  */
  void DiracDeterminantBase::acceptMove(ParticleSet& P, int iat)
  {
    PhaseValue += evaluatePhase(curRatio);
    LogValue +=std::log(std::abs(curRatio));
    UpdateTimer.start();

    InverseUpdateByRow(psiM,psiV,workV1,workV2,WorkingIndex,curRatio);

    if(UpdateMode==ORB_PBYP_PARTIAL)
      simd::copy_n(vGL.data(), BlockSize, mGL[WorkingIndex].data());

    UpdateTimer.stop();
    curRatio=1.0;
  }

  void DiracDeterminantSoA::recompute(ParticleSet& P)
  {
    for(size_t i=0,iat=FirstIndex; i<NumPtcls; ++i,++iat)
      Phi->evaluate(P, iat, VectorViewer<ValueType>(psiM[i]), mGL[i]); 
#ifdef MIXED_PRECISION
    psiM_hp = psiM;
    ParticleSet::Scalar_t PhaseValue_hp;
    LogValue = InvertWithLog(psiM_hp.data(),NumPtcls,NumOrbitals,WorkSpace_hp.data(),Pivot.data(),PhaseValue_hp);
    psiM = psiM_hp;
    PhaseValue = PhaseValue_hp;
#else
    LogValue = InvertWithLog(psiM.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);
#endif
  }


  DiracDeterminantSoA::RealType
    DiracDeterminantSoA::registerData(ParticleSet& P, PooledData<RealType>& buf)
    {
      if(NP == 0) 
      {
        workV1.resize(NumOrbitals);
        workV2.resize(NumOrbitals);
        NP=P.getTotalNum();
      }
      LogValue=evaluateLog(P,P.G,P.L);
      //add the data: determinant, inverse, gradient and laplacians
      buf.add(psiM.first_address(),psiM.last_address());
      buf.add(mGL.first_address(),mGL.last_address());
      buf.add(LogValue);
      buf.add(PhaseValue);
      return LogValue;
    }

  DiracDeterminantSoA::RealType DiracDeterminantBase::updateBuffer(ParticleSet& P,
      PooledData<RealType>& buf, bool fromscratch)
  {
    if(fromscratch)
    {
      LogValue=evaluateLog(P,P.G,P.L);
    }
    else
    {
      mGradType grad;
      mValueType lap;
      if(UpdateMode == ORB_PBYP_RATIO) //ratio only method does not compute GL
      {
        for(size_t i=0; i<NP; ++i)
        {
          VectorView arow(psiM[i]);
          Phi->evaluate(P, i, arow, mGL[i]); 
          computeGL(arow,mGL[i],grad,lap);
          P.G[iat]+=grad;
          P.L[iat]+=lap-dot(grad,grad);
        }
      }
      else
      {
        for(size_t i=0,iat=FirstIndex; i<NumPtcls; ++i,++iat)
        {
          computeGL(VectorView<ValueType>(psiM[i]),mGL[i],grad,lap);
          P.G[iat]+=grad;
          P.L[iat]+=lap-dot(grad,grad);
        }
      }
    }

    BufferTimer.start();
    buf.put(psiM.first_address(),psiM.last_address());
    buf.put(mGL.first_address(),mGL.last_address());
    buf.put(LogValue);
    buf.put(PhaseValue);
    BufferTimer.stop();
    return LogValue;
  }

  void DiracDeterminantSoA::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf)
  {
    buf.get(psiM.first_address(),psiM.last_address());
    buf.get(mGL.first_address(),mGL.last_address());
    buf.get(LogValue);
    buf.get(PhaseValue);
  }


}
