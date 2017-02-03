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
    
    
#ifndef QMCPLUSPLUS_SOA_LINEARCOMIBINATIONORBITALSET_TEMP_H
#define QMCPLUSPLUS_SOA_LINEARCOMIBINATIONORBITALSET_TEMP_H

#include "QMCWaveFunctions/SPOSetBase.h"
#include <Numerics/MatrixOperators.h>

namespace qmcplusplus
{

#if 0
  struct NoCuspCorrection
  {
    static const size_t BasisSetSize=0;

    inline void setBasisSetSize(size_t nbs) { }

    inline NoCuspCorrecton* makeClone() const 
    {
      return new NoCorrection;
    }

    template<typename VGL>
      inline void evaluateVGL(const ParticleSet& P, int iat, VGL& vgl, bool trialMove) {}

    template<typename VA>
      inline void evaluateV(const ParticleSet& P, int iat, VA& vals) {}
  };
#endif

  /** class to handle linear combinations of basis orbitals used to evaluate the Dirac determinants.
   *
   * LCOrbitalSet stands for (L)inear(C)ombinationOrbitals
   * Any single-particle orbital \f$ \psi_j \f$ that can be represented by
   * \f[
   * \psi_j ({\bf r}_i) = \sum_k C_{jk} \phi_{k}({\bf r}_i),
   * \f]
   * where \f$ \phi_{k} \f$ is the k-th basis.
   * A templated version is LCOrbitals.
   */
  template<typename BS>
    struct LCAOrbitalSet: public SPOSetBase
  {
    ///level of printing
    int ReportLevel;
    ///pointer to the basis set
    BS* myBasisSet;
    ///Temp(2+DIM,BasisSetSize) : Row index=V,L,Gx,Gy,Gz
    ValueMatrix_t Temp; 
    ///Tempv(2+DIM,OrbitalSetSize) Tempv=C*Temp
    ValueMatrix_t Tempv;
    /** constructor
     * @param bs pointer to the BasisSet
     * @param id identifier of this LCOrbitalSet
     */
    LCAOrbitalSet(BS* bs=nullptr,int rl=0): myBasisSet(nullptr),ReportLevel(rl)
    {
      if(bs != nullptr) setBasisSet(bs);
    }
    /** destructor
     *
     * BasisSet is deleted by the object with ID == 0
     */
    ~LCAOrbitalSet() {}

    SPOSetBase* makeClone() const
    {
      SoaLCOrbitalSet<BS>* myclone = new SoaLCOrbitalSet<BS>(*this);
      myclone->myBasisSet = myBasisSet->makeClone();
      return myclone;
    }

    ///reset
    void resetParameters(const opt_variables_type& active)
    {
      myBasisSet->resetParameters(active);
    }

    ///reset the target particleset
    void resetTargetParticleSet(ParticleSet& P)
    {
      myBasisSet->resetTargetParticleSet(P);
    }

    /** set the OrbitalSetSize
    */
    void setOrbitalSetSize(int norbs)
    {
      OrbitalSetSize=norbs;
      Tempv.resize(OHMMS_DIM+2,OrbitalSetSize);
    }

    /** set the basis set
    */
    void setBasisSet(BS* bs)
    {
      myBasisSet=bs;
      BasisSetSize=myBasisSet->getBasisSetSize();
      Temp.resize(OHMMS_DIM+2,BasisSetSize);
    }

    /** return the size of the basis set
    */
    inline int getBasisSetSize() const
    {
      return (myBasisSet==nullptr)? 0: myBasisSet->getBasisSetSize();
    }

    inline void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
    {
      VectorViewer<value_type> vTemp(Temp[0],BasisSetSize);
      myBasisSet->evaluateV(P,iat,vTemp);
      simd::gemv(C,Temp[0],psi.data());
    }

    inline void
      evaluate(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
      {
        const bool trialmove=true;
        BS->evaluateVGL(P,iat,Temp,trialmove);
        MatrixOperators::product_ABt(Temp,C,Tempv);
        simd::copy_n(Tempv[0],OrbitalSetSize,psi.data());
        simd::copy_n(Tempv[1],OrbitalSetSize,d2psi.data());
        for(size_t j=0; j<OrbitalSetSize; j++)
        {
          dpsi[j][0]=Tempv[2][j];
          dpsi[j][1]=Tempv[3][j];
          dpsi[j][2]=Tempv[4][j];
        }
      }

    inline void
      evaluate(const ParticleSet& P, int iat,
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

    void evaluate_notranspose(const ParticleSet& P, int first, int last,
        ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
    {
      const bool curpos=false;
      for(size_t i=0, iat=first; iat<last; i++,iat++)
      {
        BS->evaluateVGL(P,iat,Temp,curpos);
        MatrixOperators::product_ABt(Temp,C,Tempv);
        simd::copy_n(Tempv[0],OrbitalSetSize,logdet[i]);
        simd::copy_n(Tempv[1],OrbitalSetSize,d2logdet[i]);
        for(int j=0; j<OrbitalSetSize; j++)
        {
          dlogdet[i][j][0]=Tempv[2][j];
          dlogdet[i][j][1]=Tempv[3][j];
          dlogdet[i][j][2]=Tempv[4][j];
        }
      }
    }

    void evaluate_notranspose(const ParticleSet& P, int first, int last,
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

    void evaluate_notranspose(const ParticleSet& P, int first, int last
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

    void evaluateThirdDeriv(const ParticleSet& P, int first, int last
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

  };
}
#endif
