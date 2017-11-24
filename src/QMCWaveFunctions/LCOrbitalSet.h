//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_LINEARCOMIBINATIONORBITALSET_TEMP_H
#define QMCPLUSPLUS_LINEARCOMIBINATIONORBITALSET_TEMP_H

#include "QMCWaveFunctions/SPOSetBase.h"
#include <Numerics/MatrixOperators.h>

namespace qmcplusplus
{

/** decalaration of generic class to handle a linear-combination of basis function*/
template<class BS, bool IDENTITY>
class LCOrbitalSet: public SPOSetBase
{
};

template<class BS>
class LCOrbitalSet<BS,true>: public SPOSetBase
{

public:

  ///level of printing
  int ReportLevel;
  ///pointer to the basis set
  BS* myBasisSet;
  ValueMatrix_t Temp;
  ValueMatrix_t Tempv;

  /** constructor
   * @param bs pointer to the BasisSet
   * @param id identifier of this LCOrbitalSet
   */
  LCOrbitalSet(BS* bs=0,int rl=0): myBasisSet(0), ReportLevel(rl)
  {
    if(bs)
      setBasisSet(bs);
  }

  /** destructor
   *
   * BasisSet is deleted by the object with ID == 0
   */
  ~LCOrbitalSet() {}

  SPOSetBase* makeClone() const
  {
    LCOrbitalSet<BS,true>* myclone = new LCOrbitalSet<BS,true>(*this);
    myclone->myBasisSet = myBasisSet->makeClone();
    myclone->IsCloned=true;
    return myclone;
  }

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
    Tempv.resize(OrbitalSetSize,5);
  }

  /** set the basis set
   */
  void setBasisSet(BS* bs)
  {
    myBasisSet=bs;
    BasisSetSize=myBasisSet->getBasisSetSize();
    Temp.resize(BasisSetSize,5);
  }

  /** return the size of the basis set
   */
  inline int getBasisSetSize() const
  {
    return (myBasisSet==0)? 0: myBasisSet->getBasisSetSize();
  }

  inline void
  evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
  {
    myBasisSet->evaluateForPtclMove(P,iat);
    for(int j=0 ; j<OrbitalSetSize; j++)
      psi[j] = myBasisSet->Phi[j];
  }

  inline void
  evaluate(const ParticleSet& P, int iat,
           ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
  {
    myBasisSet->evaluateAllForPtclMove(P,iat);
    for(int j=0; j<OrbitalSetSize; j++)
      psi[j]=myBasisSet->Phi[j];
    for(int j=0; j<OrbitalSetSize; j++)
      dpsi[j]=myBasisSet->dPhi[j];
    for(int j=0; j<OrbitalSetSize; j++)
      d2psi[j]=myBasisSet->d2Phi[j];
  }

  inline void
  evaluate(const ParticleSet& P, int iat,
           ValueVector_t& psi, GradVector_t& dpsi,
           HessVector_t& grad_grad_psi)
  {
    myBasisSet->evaluateForPtclMoveWithHessian(P,iat);
    for(int j=0; j<OrbitalSetSize; j++)
      psi[j]=myBasisSet->Phi[j];
    for(int j=0; j<OrbitalSetSize; j++)
      dpsi[j]=myBasisSet->dPhi[j];
    for(int j=0; j<OrbitalSetSize; j++)
      grad_grad_psi[j]=myBasisSet->grad_grad_Phi[j];
  }

  void evaluate_notranspose(const ParticleSet& P, int first, int last,
                            ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
  {
    for(int i=0, iat=first; iat<last; i++,iat++)
    {
      myBasisSet->evaluateForWalkerMove(P,iat);
      std::copy(myBasisSet->Phi.data(),myBasisSet->Phi.data()+OrbitalSetSize,logdet[i]);
      std::copy(myBasisSet->dPhi.data(),myBasisSet->dPhi.data()+OrbitalSetSize,dlogdet[i]);
      std::copy(myBasisSet->d2Phi.data(),myBasisSet->d2Phi.data()+OrbitalSetSize,d2logdet[i]);
    }
  }

  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    APP_ABORT("Need specialization of LCOrbitalSet<BS,true>::evaluate_notranspose() for grad_grad_logdet. \n");
  }

  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet)
  {
    APP_ABORT("Need specialization of LCOrbitalSet<BS,true>::evaluate_notranspose() for grad_grad_grad_logdet. \n");
  }

};

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
template<class BS>
class LCOrbitalSet<BS,false>: public SPOSetBase
{

public:

  ///level of printing
  int ReportLevel;
  ///pointer to the basis set
  BS* myBasisSet;

  ///temporary matricies
  ValueMatrix_t Temp;
  ValueMatrix_t Tempv;

  ///algorithm switch
  int Algo;

  /** constructor
   * @param bs pointer to the BasisSet
   * @param id identifier of this LCOrbitalSet
   */
  LCOrbitalSet(BS* bs=0,int rl=0, std::string algorithm=""): myBasisSet(0), ReportLevel(rl)
  {
    if(algorithm=="legacy_gemv")
    {
      Algo=0;
      /// input sample: <coefficient size="57" id="updetC" algorithm="legacy_gemv">
      app_log() << "SPO coefficients are applied using the legacy algorithm" << std::endl;
    }
    else
    {
      Algo=1;
      app_log() << "SPO coefficients are applied by a single dgemm (BLAS3)" << std::endl;
    }
    if(bs)
      setBasisSet(bs);
  }

  /** destructor
   *
   * BasisSet is deleted by the object with ID == 0
   */
  ~LCOrbitalSet() {}

  SPOSetBase* makeClone() const
  {
    LCOrbitalSet<BS,false>* myclone = new LCOrbitalSet<BS,false>(*this);
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
    Tempv.resize(5,OrbitalSetSize);
  }

  /** set the basis set
   */
  void setBasisSet(BS* bs)
  {
    myBasisSet=bs;
    BasisSetSize=myBasisSet->getBasisSetSize();
    Temp.resize(5,BasisSetSize);
  }

  /** return the size of the basis set
   */
  inline int getBasisSetSize() const
  {
    return (myBasisSet==0)? 0: myBasisSet->getBasisSetSize();
  }

  inline void
  evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
  {
    myBasisSet->evaluateForPtclMove(P,iat);
    simd::gemv(*C,myBasisSet->Phi.data(),psi.data());
  }

  inline void
  evaluate(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
  {
    myBasisSet->evaluateAllForPtclMove(P,iat);

    if(Algo==1)
    {
      // Fast algorithm, using a single dgemm
      // Temp stores SOA BasisSet
      simd::copy(Temp.data(), myBasisSet->Phi.data(), BasisSetSize);
      simd::copy(Temp.data() + BasisSetSize, myBasisSet->d2Phi.data(), BasisSetSize);
      for(int i=0; i<BasisSetSize; i++)
      {
        Temp[2][i]=myBasisSet->dPhi[i][0];
        Temp[3][i]=myBasisSet->dPhi[i][1];
        Temp[4][i]=myBasisSet->dPhi[i][2];
      }
    
      MatrixOperators::product_ABt(Temp,*C,Tempv);
    
      // Tempv stores SOA SPOset
      simd::copy(psi.data(), Tempv.data(), OrbitalSetSize);
      simd::copy(d2psi.data(), Tempv.data() + OrbitalSetSize, OrbitalSetSize);
      for(int i=0; i<OrbitalSetSize; i++)
      {
        dpsi[i][0]=Tempv[2][i];
        dpsi[i][1]=Tempv[3][i];
        dpsi[i][2]=Tempv[4][i];
      }
    }
    else
    {
      // legacy algorithm
      simd::gemv(*C,myBasisSet->Phi.data(),psi.data());
      simd::gemv(*C,myBasisSet->dPhi.data(),dpsi.data());
      simd::gemv(*C,myBasisSet->d2Phi.data(),d2psi.data());
    }
  }

  inline void
  evaluate(const ParticleSet& P, int iat,
           ValueVector_t& psi, GradVector_t& dpsi,
           HessVector_t& grad_grad_psi)
  {
    myBasisSet->evaluateForPtclMoveWithHessian(P,iat);
    simd::gemv(*C,myBasisSet->Phi.data(),psi.data());
    simd::gemv(*C,myBasisSet->dPhi.data(),dpsi.data());
    simd::gemv(*C,myBasisSet->grad_grad_Phi.data(),grad_grad_psi.data());
//#if defined(USE_BLAS2)
//      MatrixOperators::product(C,myBasisSet->Phi.data(),psi.data());
//      MatrixOperators::product(C,myBasisSet->dPhi.data(),dpsi.data());
//      MatrixOperators::product(C,myBasisSet->grad_grad_Phi.data(),grad_grad_psi.data());
//#else
//      const ValueType* restrict cptr=C->data();
//      const typename BS::ValueType* restrict pptr=myBasisSet->Phi.data();
//      const typename BS::HessType* restrict d2ptr=myBasisSet->grad_grad_Phi.data();
//      const typename BS::GradType* restrict dptr=myBasisSet->dPhi.data();
//#pragma ivdep
//      for(int j=0,kk=0; j<OrbitalSetSize; j++) {
//        register ValueType res=0.0;
//        register GradType dres;
//        register HessType hess;
//        for(int b=0; b<BasisSetSize; b++,kk++)
//        {
//          res += cptr[kk]*pptr[b];
//          hess += cptr[kk]*d2ptr[b];
//          dres += cptr[kk]*dptr[b];
//        }
//        psi[j]=res; dpsi[j]=dres; grad_grad_psi[j]=hess;
//      }
//#endif
  }

  void evaluate_notranspose(const ParticleSet& P, int first, int last,
                            ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
  {
    const ValueType* restrict cptr=C->data();
    for(int i=0,ij=0, iat=first; iat<last; i++,iat++)
    {
      myBasisSet->evaluateForWalkerMove(P,iat);

      if(Algo==1)
      {
        // Fast algorithm, using a single dgemm
        // Temp stores SOA BasisSet
        simd::copy(Temp.data(), myBasisSet->Phi.data(), BasisSetSize);
        simd::copy(Temp.data() + BasisSetSize, myBasisSet->d2Phi.data(), BasisSetSize);
        for(int k=0; k<BasisSetSize; k++)
        {
          Temp[2][k]=myBasisSet->dPhi[k][0];
          Temp[3][k]=myBasisSet->dPhi[k][1];
          Temp[4][k]=myBasisSet->dPhi[k][2];
        }

        MatrixOperators::product_ABt(Temp,*C,Tempv);

        // Tempv stores SOA SPOset
        simd::copy(logdet[i], Tempv.data(), OrbitalSetSize);
        simd::copy(d2logdet[i], Tempv.data() + OrbitalSetSize, OrbitalSetSize);
        for(int j=0; j<OrbitalSetSize; j++)
        {
          dlogdet[i][j][0]=Tempv[2][j];
          dlogdet[i][j][1]=Tempv[3][j];
          dlogdet[i][j][2]=Tempv[4][j];
        }
      }
      else
      {
        // legacy algorithm
        MatrixOperators::product(*C,myBasisSet->Phi,logdet[i]);
        MatrixOperators::product(*C,myBasisSet->d2Phi,d2logdet[i]);
        const typename BS::GradType* restrict dptr=myBasisSet->dPhi.data();
        for(int j=0,jk=0; j<OrbitalSetSize; j++)
        {
          register GradType dres;
          for(int b=0; b<BasisSetSize; ++b)
            dres +=  cptr[jk++]*dptr[b];
          dlogdet(ij)=dres;
          ++ij;
        }
      }
    }
  }

  void evaluate_notranspose(const ParticleSet& P, int first, int last,
                            ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    const ValueType* restrict cptr=C->data();
#pragma ivdep
    for(int i=0,ij=0, iat=first; iat<last; i++,iat++)
    {
      myBasisSet->evaluateWithHessian(P,iat);
      MatrixOperators::product(*C,myBasisSet->Phi,logdet[i]);
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
  }

  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet)
  {
    const ValueType* restrict cptr=C->data();
#pragma ivdep
    for(int i=0,ij=0, iat=first; iat<last; i++,iat++)
    {
      myBasisSet->evaluateWithThirdDeriv(P,iat);
      MatrixOperators::product(*C,myBasisSet->Phi,logdet[i]);
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
  }

  /** evaluate psiM for virtual moves
   *
   * For the i-th virtual move and the j-th orbital,
   * \f$ psiM(i,j)= \sum_k phiM(i,k)*C(j,k) \f$
   */
  void evaluateValues(const ParticleSet& P, ValueMatrix_t& psiM)
  {
    ValueMatrix_t phiM(P.getTotalNum(),BasisSetSize);
    myBasisSet->evaluateValues(P,phiM);
    MatrixOperators::product_ABt(phiM,*C,psiM);
    //for(int i=0; i<psiM.rows(); ++i)
    //  for(int j=0; j<psiM.cols(); ++j)
    //    psiM(i,j)=simd::dot(C[j],phiM[i],BasisSetSize);
  }

  void evaluateThirdDeriv(const ParticleSet& P, int first, int last
                          , GGGMatrix_t& grad_grad_grad_logdet)
  {
    const ValueType* restrict cptr=C->data();
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
  }

};
}
#endif


