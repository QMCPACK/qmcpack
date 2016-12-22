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
  { 
    Need2Compute4PbyP=true;
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
    psiV.resize(norb);

#ifdef MIXED_PRECISION
    size_t nh=getAlignedSize<mValueType>(norb);
    psiM_hp.resize(nel,nh);
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
    const T* restrict gx_p=gl_v.data(0);
    const T* restrict gy_p=gl_v.data(1);
    const T* restrict gz_p=gl_v.data(2);
    const T* restrict l_p=gl_v.data(3);
    lap=czero;
    T gx=czero, gy=czero, gz=czero,l=czero;
    const int n=gl_v.size();
//#pragma omp simd reduction(+:gx,gy,gz,l)
    for(size_t i=0; i<n; ++i)
    {
      gx +=row[i]*gx_p[i];
      gy +=row[i]*gy_p[i];
      gz +=row[i]*gz_p[i];
      l+=row[i]*l_p[i];
    }
    grad[0]=gx;
    grad[1]=gy;
    grad[2]=gz;
    lap=l;
    //int four=4;
    //int na=gl_v.size();
    //int lda=gl_v.capacity();
    //T y[]={czero,czero,czero,czero};
    //BLAS::gemv('T',na,four,cone,gl_v.data(),lda,row,1,czero,y,1);
    //grad[0]=y[0]; grad[1]=y[1]; grad[2]=y[2]; lap=y[3];
  }

  template<typename T>
  inline TinyVector<T,3> computeG(T* row, VectorSoaContainer<T,4>& gl_v)
  {
    constexpr T czero(0);
    constexpr T cone(1);
    const T* restrict gx_p=gl_v.data(0);
    const T* restrict gy_p=gl_v.data(1);
    const T* restrict gz_p=gl_v.data(2);
    T gx=czero, gy=czero, gz=czero;
    const int n=gl_v.size();
//#pragma omp simd reduction(+:gx,gy,gz)
    for(size_t i=0; i<n; ++i)
    {
      gx+=row[i]*gx_p[i];
      gy+=row[i]*gy_p[i];
      gz+=row[i]*gz_p[i];
    }
    return TinyVector<T,3>(gx,gy,gz);
    //grad[0]=gx;
    //grad[1]=gy;
    //grad[2]=gz;
    //int three=3;
    //int na=gl_v.size();
    //int lda=gl_v.capacity();
    //T y[]={czero,czero,czero,czero};
    //BLAS::gemv('T',na,three,cone,gl_v.data(),lda,row,1,czero,y,1);
    //return TinyVector<T,3>(y[0],y[1],y[2]);
  }
#if 0
  template<typename S, typename T>
  inline void computeGL(S* row, VectorSoaContainer<S,4>& gl_v, TinyVector<T,3>& grad, T& lap)
  {
    CONSTEXPR T czero(0);
    CONSTEXPR T cone(1);
    int four=4;
    int na=gl_v.size();
    int lda=gl_v.capacity();
    T y[]={czero,czero,czero,czero};
    BLAS::gemv('T',na,four,cone,gl_v.data(),lda,row,1,czero,y,1);
    grad[0]=y[0]; grad[1]=y[1]; grad[2]=y[2]; lap=y[3];
  }

  template<typename S, typename T>
  inline void computeG(S* row, VectorSoaContainer<S,4>& gl_v, TinyVector<T,3>& grad)
  {
    constexpr T czero(0);
    constexpr T cone(1);
    int three=3;
    int na=gl_v.size();
    int lda=gl_v.capacity();
    T y[]={czero,czero,czero,czero};
    BLAS::gemv('T',na,three,cone,gl_v.data(),lda,row,1,czero,y,1);
    grad[0]=y[0]; grad[1]=y[1]; grad[2]=y[2];
  }
#endif

  DiracDeterminantSoA::RealType
    DiracDeterminantSoA::evaluateLog(ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
    {
      if(NumPtcls==1)
      { //waste of good code
        VectorViewer<ValueType> arow(psiM[0],1);
        Phi->evaluateVGL(P, FirstIndex, arow, mGL[0],false); 
        ValueType det=psiM(0,0);
        ValueType y=(RealType)1.0/det;
        psiM(0,0)=y;
        GradType rv = y*GradType(mGL[0].data(0),OHMMS_DIM+1);
        ValueType lap= mGL[0][0][OHMMS_DIM];
        G(FirstIndex) += rv;
        L(FirstIndex) += y*lap - dot(rv,rv);

        LogValue = evaluateLogAndPhase(det,PhaseValue);
      }
      else
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
      }
      psiM_temp = psiM;
      return LogValue;
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

  DiracDeterminantBase::GradType
    DiracDeterminantSoA::evalGrad(ParticleSet& P, int iat)
    {
      WorkingIndex = iat-FirstIndex;
      if(Need2Compute4PbyP)//need to compute mGL
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

      //GradType rv=simd::dot(psiM[WorkingIndex],dpsiV.data(),NumOrbitals);
      GradType rv= computeG(psiM[WorkingIndex],vGL);
      grad_iat += f * rv;
      //grad_iat += (cone/curRatio) * rv;
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

    detEng.updateRow(psiM,psiV.data(),WorkingIndex,curRatio);

    if(UpdateMode==ORB_PBYP_PARTIAL)
      simd::copy_n(vGL.data(), BlockSize, mGL[WorkingIndex].data());

    curRatio=ValueType(1);
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


  DiracDeterminantSoA::RealType
    DiracDeterminantSoA::registerData(ParticleSet& P, PooledData<RealType>& buf)
    {
      NP=P.getTotalNum();
      LogValue=evaluateLog(P,P.G,P.L);

      //add the data: determinant, inverse, gradient and laplacians
      buf.add(psiM.first_address(),psiM.last_address());
      //if(!Need2Compute4PbyP) buf.add(memoryPool.data(),memoryPool.data()+memoryPool.size());
      buf.add(LogValue);
      buf.add(PhaseValue);
      return LogValue;
    }

  /** done PbyP update, prepare for the measurements */
  void DiracDeterminantSoA::updateAfterSweep(ParticleSet& P,
      ParticleSet::ParticleGradient_t& G,
      ParticleSet::ParticleLaplacian_t& L)
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


  DiracDeterminantSoA::RealType 
    DiracDeterminantSoA::updateBuffer(ParticleSet& P, PooledData<RealType>& buf, bool fromscratch)
    {
      //if(fromscratch)
        LogValue=evaluateLog(P,P.G,P.L);
      //else
      //  updateAfterSweep(P,P.G,P.L);

      buf.put(psiM.first_address(),psiM.last_address());
      //if(!Need2Compute4PbyP) buf.put(memoryPool.data(),memoryPool.data()+memoryPool.size());
      buf.put(LogValue);
      buf.put(PhaseValue);
      return LogValue;
    }

  void DiracDeterminantSoA::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf)
  {
    buf.get(psiM.first_address(),psiM.last_address());
    //if(!Need2Compute4PbyP) buf.get(memoryPool.data(),memoryPool.data()+memoryPool.size());
    buf.get(LogValue);
    buf.get(PhaseValue);
  }


}
