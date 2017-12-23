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
    
    
#include "QMCWaveFunctions/lcao/LCAOrbitalSet.h"
#include <Numerics/MatrixOperators.h>

namespace qmcplusplus
{
  LCAOrbitalSet::LCAOrbitalSet(basis_type* bs,int rl): 
    myBasisSet(nullptr),ReportLevel(rl)
  {
    if(bs != nullptr) setBasisSet(bs);
  }

  LCAOrbitalSet::~LCAOrbitalSet() {}

  void LCAOrbitalSet::setBasisSet(basis_type* bs)
  {
    myBasisSet=bs;
    BasisSetSize=myBasisSet->getBasisSetSize();
    Temp.resize(BasisSetSize);
  }

  SPOSetBase* LCAOrbitalSet::makeClone() const
  {
    LCAOrbitalSet* myclone = new LCAOrbitalSet(*this);
    myclone->myBasisSet = myBasisSet->makeClone();
    myclone->IsCloned=true;
    return myclone;
  }

  void LCAOrbitalSet::evaluate(const ParticleSet& P, 
      int iat, ValueVector_t& psi)
  {
    if(Identity)
    { //PAY ATTENTION TO COMPLEX
      myBasisSet->evaluateV(P,iat,psi.data());
    }
    else
    {
      VectorViewer<ValueType> vTemp(Temp.data(0),BasisSetSize);
      myBasisSet->evaluateV(P,iat,vTemp.data());
      simd::gemv(*C,Temp.data(0),psi.data());
    }
  }

  /** Find a better place for other user classes, Matrix should be padded as well */
  template<typename T,unsigned D>
    inline void Product_ABt(const VectorSoaContainer<T,D>& A, const Matrix<T>& B, VectorSoaContainer<T,D>& C)
    {
      CONSTEXPR char transa = 't';
      CONSTEXPR char transb = 'n';
      CONSTEXPR T zone(1);
      CONSTEXPR T zero(0);
      BLAS::gemm(transa, transb, B.rows(), D, B.cols(),
          zone, B.data(), B.cols(), A.data(), A.capacity(),
          zero, C.data(), C.capacity());
    }

  inline void LCAOrbitalSet::evaluate_vgl_impl(const vgl_type& temp,
      ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi) const
  {
    simd::copy_n(temp.data(0),OrbitalSetSize,psi.data());
    const ValueType* restrict gx=temp.data(1);
    const ValueType* restrict gy=temp.data(2);
    const ValueType* restrict gz=temp.data(3);
    for(size_t j=0; j<OrbitalSetSize; j++)
    {
      dpsi[j][0]=gx[j];
      dpsi[j][1]=gy[j];
      dpsi[j][2]=gz[j];
    }
    simd::copy_n(temp.data(4),OrbitalSetSize,d2psi.data());
  }

  void LCAOrbitalSet::evaluate(const ParticleSet& P, int iat,
      ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
    {
      //TAKE CARE OF IDENTITY
      myBasisSet->evaluateVGL(P,iat,Temp);
      if(Identity) 
        evaluate_vgl_impl(Temp,psi,dpsi,d2psi);
      else
      {
        Product_ABt(Temp,*C,Tempv);
        evaluate_vgl_impl(Tempv,psi,dpsi,d2psi);
      }
    }

  void LCAOrbitalSet::evaluateVGL(const ParticleSet& P, int iat, 
      VGLVector_t vgl)
  {
    if(Identity)
      myBasisSet->evaluateVGL(P,iat,vgl);
    else
    {
      myBasisSet->evaluateVGL(P,iat,Temp);
      Product_ABt(Temp,*C,vgl);
    }
  }

  void LCAOrbitalSet::evaluate(const ParticleSet& P, int iat,
        ValueVector_t& psi, GradVector_t& dpsi,
        HessVector_t& grad_grad_psi)
    {
#if 0
        myBasisSet->evaluateForPtclMoveWithHessian(P,iat);
        simd::gemv(C,myBasisSet->Phi.data(),psi.data());
        simd::gemv(C,myBasisSet->dPhi.data(),dpsi.data());
        simd::gemv(C,myBasisSet->grad_grad_Phi.data(),grad_grad_psi.data());
#endif
      }

  /* implement using gemm algorithm */
  inline void LCAOrbitalSet::evaluate_vgl_impl(const vgl_type& temp, int i,
      ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet) const
  {
    simd::copy_n(temp.data(0),OrbitalSetSize,logdet[i]);
    const ValueType* restrict gx=temp.data(1);
    const ValueType* restrict gy=temp.data(2);
    const ValueType* restrict gz=temp.data(3);
    for(size_t j=0; j<OrbitalSetSize; j++)
    {
      dlogdet[i][j][0]=gx[j];
      dlogdet[i][j][1]=gy[j];
      dlogdet[i][j][2]=gz[j];
    }
    simd::copy_n(temp.data(4),OrbitalSetSize,d2logdet[i]);
  }

  void LCAOrbitalSet::evaluate_notranspose(const ParticleSet& P, int first, int last,
      ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
  {
    if(Identity)
    {
      for(size_t i=0, iat=first; iat<last; i++,iat++)
      {
        myBasisSet->evaluateVGL(P,iat,Temp);
        evaluate_vgl_impl(Temp,i,logdet,dlogdet,d2logdet);
      }
    }
    else
    {
      for(size_t i=0, iat=first; iat<last; i++,iat++)
      {
        myBasisSet->evaluateVGL(P,iat,Temp);
        Product_ABt(Temp,*C,Tempv);
        evaluate_vgl_impl(Tempv,i,logdet,dlogdet,d2logdet);
      }
    }
  }

  void LCAOrbitalSet::evaluate_notranspose(const ParticleSet& P, int first, int last,
      ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
#if 0
    const ValueType* restrict cptr=C.data();
#pragma ivdep
    for(int i=0,ij=0, iat=first; iat<last; i++,iat++)
    {
      myBasisSet->evaluateWithHessian(P,iat);
      MatrixOperators::product(C,myBasisSet->Phi,logdet[i]);
      const typename BS::GradType* restrict dptr=myBasisSet->dPhi.data();
      const typename BS::HessType* restrict d2ptr=myBasisSet->grad_grad_Phi.data();
      for(int j=0,jk=0; j<OrbitalSetSize; j++)
      {
        register GradType dres;
        register HessType d2res;
        for(int b=0; b<BasisSetSize; ++b,++jk)
        {
          dres +=  cptr[jk]*dptr[b];
          d2res +=  cptr[jk]*d2ptr[b];
        }
        dlogdet(ij)=dres;
        grad_grad_logdet(ij)=d2res;
        ++ij;
      }
    }
#endif
  }

  void LCAOrbitalSet::evaluate_notranspose(const ParticleSet& P, int first, int last
      , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet)
  {
#if 0
    const ValueType* restrict cptr=C.data();
#pragma ivdep
    for(int i=0,ij=0, iat=first; iat<last; i++,iat++)
    {
      myBasisSet->evaluateWithThirdDeriv(P,iat);
      MatrixOperators::product(C,myBasisSet->Phi,logdet[i]);
      const typename BS::GradType* restrict dptr=myBasisSet->dPhi.data();
      const typename BS::HessType* restrict d2ptr=myBasisSet->grad_grad_Phi.data();
      const typename BS::GGGType* restrict gggptr=myBasisSet->grad_grad_grad_Phi.data();
      for(int j=0,jk=0; j<OrbitalSetSize; j++)
      {
        register GradType dres;
        register HessType d2res;
        register GGGType gggres;
        for(int b=0; b<BasisSetSize; ++b)
        {
          dres +=  cptr[jk]*dptr[b];
          d2res +=  cptr[jk]*d2ptr[b];
          gggres[0] +=  cptr[jk]*(gggptr[b])[0];
          gggres[1] +=  cptr[jk]*(gggptr[b])[1];
          gggres[2] +=  cptr[jk++]*(gggptr[b])[2];
        }
        dlogdet(ij)=dres;
        grad_grad_logdet(ij)=d2res;
        grad_grad_grad_logdet(ij)=gggres;
        ++ij;
      }
    }
#endif
  }

  void LCAOrbitalSet::evaluateThirdDeriv(const ParticleSet& P, int first, int last
      , GGGMatrix_t& grad_grad_grad_logdet)
  {
#if 0
    const ValueType* restrict cptr=C.data();
#pragma ivdep
    for(int i=0,ij=0, iat=first; iat<last; i++,iat++)
    {
      myBasisSet->evaluateThirdDerivOnly(P,iat);
      const typename BS::GGGType* restrict gggptr=myBasisSet->grad_grad_grad_Phi.data();
      for(int j=0,jk=0; j<OrbitalSetSize; j++)
      {
        register GGGType gggres;
        for(int b=0; b<BasisSetSize; ++b)
        {
          gggres[0] +=  cptr[jk]*(gggptr[b])[0];
          gggres[1] +=  cptr[jk]*(gggptr[b])[1];
          gggres[2] +=  cptr[jk++]*(gggptr[b])[2];
        }
        grad_grad_grad_logdet(ij)=gggres;
        ++ij;
      }
    }
#endif
  }
}
