//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_LINEARCOMIBINATIONORBITALSETWITHCORRECTION_TEMP_H
#define QMCPLUSPLUS_LINEARCOMIBINATIONORBITALSETWITHCORRECTION_TEMP_H

#include <vector>
#include <cstdlib>
#include "QMCWaveFunctions/SPOSetBase.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/MatrixOperators.h"
#include "Message/Communicate.h"
#include "QMCWaveFunctions/Experimental/CorrectingOrbitalBasisSet.h"
#include "QMCWaveFunctions/MolecularOrbitals/NGOBuilder.h"
#include "QMCWaveFunctions/LCOrbitalSet.h"
#include "QMCWaveFunctions/Experimental/CuspCorr.h"
#include "Numerics/OptimizableFunctorBase.h"
#include "Numerics/OneDimQuinticSpline.h"
#include "QMCWaveFunctions/SphericalBasisSet.h"

namespace qmcplusplus
{

// mmorales:
/* A class designed for cusp corrected LCAO-type molecular orbitals.
 * It is basically a copy of LCOrbitalSet with a few minor modifications.
 */

/** decalaration of generic class to handle a linear-combination of basis function*/
template<class BS, bool IDENTITY>
class LCOrbitalSetWithCorrection: public SPOSetBase
{
};

template<class BS>
class LCOrbitalSetWithCorrection<BS,true>: public SPOSetBase
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
   * @param id identifier of thisLCOrbitalSetWithCorrection
   */
  LCOrbitalSetWithCorrection(BS* bs=0,int rl=0): myBasisSet(0),ReportLevel(rl)
  {
    if(bs)
      setBasisSet(bs);
  }

  /** destructor
   *
   * BasisSet is deleted by the object with ID == 0
   */
  ~LCOrbitalSetWithCorrection() {}

  SPOSetBase* makeClone() const
  {
    LCOrbitalSetWithCorrection<BS,true>* myclone = new LCOrbitalSetWithCorrection<BS,true>(*this);
    myclone->setBasisSet(myBasisSet->makeClone());
    myclone->setOrbitalSetSize(OrbitalSetSize);
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


  ///** evaluate everything for the walker move
  // *
  // * Using gemm can improve the performance for a larger problem
  // */
  //void evaluate(const ParticleSet& P, int first, int last,
  //    ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet) {
  //  for(int i=0, iat=first; iat<last; i++,iat++){
  //    myBasisSet->evaluateForWalkerMove(P,iat);
  //    for(int j=0; j<OrbitalSetSize; j++) logdet(j,i)=myBasisSet->Phi[j];
  //    for(int j=0; j<OrbitalSetSize; j++) dlogdet(i,j)=myBasisSet->dPhi[j];
  //    for(int j=0; j<OrbitalSetSize; j++) d2logdet(i,j)=myBasisSet->d2Phi[j];
  //  }
  //}

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

  void evaluate(const ParticleSet& P, int first, int last,
                ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
  {
    APP_ABORT("LCOrbitalSet:evaluate() not implemented yet. \n");
  }

  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    APP_ABORT("Need specialization of LCOrbitalSetWithCorrection<BS,true>::evaluate_notranspose() for grad_grad_logdet. \n");
  }

  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet)
  {
    APP_ABORT("Need specialization of LCOrbitalSetWithCorrection::evaluate_notranspose() for grad_grad_grad_logdet. \n");
  }

  void evaluateThirdDeriv(const ParticleSet& P, int first, int last
                          , GGGMatrix_t& grad_grad_grad_logdet)
  {
    APP_ABORT("Need specialization of LCOrbitalSetWithCorrection::evaluateThirdDeriv()\n");
  }

};

/** class to handle linear combinations of basis orbitals used to evaluate the Dirac determinants.
 *
 * LCOrbitalSetWithCorrection stands for (L)inear(C)ombinationOrbitals wit correction
 * Any single-particle orbital \f$ \psi_j \f$ that can be represented by
 * \f[
 * \psi_j ({\bf r}_i) = \sum_k C_{jk} \phi_{k}({\bf r}_i) + \sum_k' C'_{jk'} \phi^{'}_{k'}({\bf r}_i),
 * \f]
 * where \f$ \phi_{k} \f$ is the k-th basis and \f$ \phi^{'}_{k'} \f$ is the k'-th basis which contains all corrected sperically symmetric functions.
 * A templated version is LCOrbitals.
 */
template<class BS>
class LCOrbitalSetWithCorrection<BS,false>: public SPOSetBase
{

public:

  typedef typename NGOBuilder::CenteredOrbitalType COT;
  typedef typename NGOBuilder::GridType GridType_;

  ///level of printing
  int ReportLevel;
  ///pointer to the basis set
  BS* myBasisSet;
  ///pointer to the basis set of corrected s functions
  CorrectingOrbitalBasisSet<COT>* corrBasisSet;
  /// vector that defines whether a center is corrected
  vector<bool> corrCenter;

  ///target ParticleSet
  ParticleSet* targetPtcl;
  ///source ParticleSet
  ParticleSet* sourcePtcl;

  ValueMatrix_t Temp;
  ValueMatrix_t Tempv;
  /** constructor
   * @param bs pointer to the BasisSet
   * @param id identifier of thisLCOrbitalSetWithCorrection
   */
  LCOrbitalSetWithCorrection(BS* bs=0, ParticleSet* els=0, ParticleSet* ions=0, int rl=0, double rc=-1.0, string cusp_info=""): myBasisSet(0),corrBasisSet(0),targetPtcl(els),sourcePtcl(ions),ReportLevel(rl),Rcut(rc),cuspInfoFile(cusp_info)
  {
// call APP_ABORT if els==0, ions==0
    if(bs)
      setBasisSet(bs);
  }

  /** destructor
   *
   * BasisSet is deleted by the object with ID == 0
   */
  ~LCOrbitalSetWithCorrection() {}

  SPOSetBase* makeClone() const
  {
    LCOrbitalSetWithCorrection<BS,false>* myclone = new LCOrbitalSetWithCorrection<BS,false>(*this);
//       myclone->myBasisSet = myBasisSet->makeClone();
    myclone->setBasisSet(myBasisSet->makeClone());
    myclone->setOrbitalSetSize(OrbitalSetSize);
    myclone->corrBasisSet = corrBasisSet->makeClone();
    return myclone;
  }

  ///reset
  void resetParameters(const opt_variables_type& active)
  {
    myBasisSet->resetParameters(active);
    corrBasisSet->resetParameters(active);
  }

  ///reset the target particleset
  void resetTargetParticleSet(ParticleSet& P)
  {
    myBasisSet->resetTargetParticleSet(P);
    corrBasisSet->resetTargetParticleSet(P);
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
    corrBasisSet->evaluateForPtclMove(P,iat);
    int numCenters=corrBasisSet->NumCenters;
    for(int j=0 ; j<OrbitalSetSize; j++)
    {
      psi[j] = simd::dot(C[j],myBasisSet->Phi.data(),BasisSetSize);
      for(int k=0 ; k<numCenters; k++)
        psi[j] += corrBasisSet->Phi[k*OrbitalSetSize+j];
    }
    //overhead of blas::gemv is big
    //MatrixOperators::product(C,myBasisSet->Phi,psi.data());
    //overhead of blas::dot is too big, better to use the inline function
    //for((int j=0 ; j<OrbitalSetSize; j++)
    //  psi[j] = BLAS::dot(BasisSetSize,C[j],myBasisSet->Phi.data());
  }

  inline void
  evaluate(const ParticleSet& P, int iat,
           ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
  {
    myBasisSet->evaluateAllForPtclMove(P,iat);
    corrBasisSet->evaluateAllForPtclMove(P,iat);
    //optimal on tungsten
    const ValueType* restrict cptr=C.data();
    int numCenters=corrBasisSet->NumCenters;
    const typename BS::ValueType* restrict pptr=myBasisSet->Phi.data();
    const typename BS::ValueType* restrict d2ptr=myBasisSet->d2Phi.data();
    const typename BS::GradType* restrict dptr=myBasisSet->dPhi.data();
#pragma ivdep
    for(int j=0; j<OrbitalSetSize; j++)
    {
      register ValueType res=0.0, d2res=0.0;
      register GradType dres;
      for(int b=0; b<BasisSetSize; b++,cptr++)
      {
        res += *cptr*pptr[b];
        d2res += *cptr*d2ptr[b];
        dres += *cptr*dptr[b];
        //res += *cptr * (*pptr++);
        //d2res += *cptr* (*d2ptr++);
        //dres += *cptr* (*dptr++);
      }
      for(int k=0 ; k<numCenters; k++)
      {
        res += corrBasisSet->Phi[k*OrbitalSetSize+j];
        d2res += corrBasisSet->d2Phi[k*OrbitalSetSize+j];
        dres += corrBasisSet->dPhi[k*OrbitalSetSize+j];
      }
      psi[j]=res;
      dpsi[j]=dres;
      d2psi[j]=d2res;
    }
    //blasI is not too good
    //for(int j=0 ; j<OrbitalSetSize; j++) {
    //  psi[j]   = dot(C[j],myBasisSet->Phi.data(),  BasisSetSize);
    //  dpsi[j]  = dot(C[j],myBasisSet->dPhi.data(), BasisSetSize);
    //  d2psi[j] = dot(C[j],myBasisSet->d2Phi.data(),BasisSetSize);
    //  //psi[j]   = dot(C[j],myBasisSet->Y[0],  BasisSetSize);
    //  //dpsi[j]  = dot(C[j],myBasisSet->dY[0], BasisSetSize);
    //  //d2psi[j] = dot(C[j],myBasisSet->d2Y[0],BasisSetSize);
    //}
  }

  inline void
  evaluate(const ParticleSet& P, int iat,
           ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& d2psi)
  {
    myBasisSet->evaluateForPtclMoveWithHessian(P,iat);
    corrBasisSet->evaluateForPtclMoveWithHessian(P,iat);
    //optimal on tungsten
    const ValueType* restrict cptr=C.data();
    int numCenters=corrBasisSet->NumCenters;
    const typename BS::ValueType* restrict pptr=myBasisSet->Phi.data();
    const typename BS::HessType* restrict d2ptr=myBasisSet->grad_grad_Phi.data();
    const typename BS::GradType* restrict dptr=myBasisSet->dPhi.data();
#pragma ivdep
    for(int j=0; j<OrbitalSetSize; j++)
    {
      register ValueType res=0.0;
      register HessType d2res=0.0;
      register GradType dres;
      for(int b=0; b<BasisSetSize; b++,cptr++)
      {
        res += *cptr*pptr[b];
        d2res += *cptr*d2ptr[b];
        dres += *cptr*dptr[b];
      }
      for(int k=0 ; k<numCenters; k++)
      {
        res += corrBasisSet->Phi[k*OrbitalSetSize+j];
        d2res += corrBasisSet->grad_grad_Phi[k*OrbitalSetSize+j];
        dres += corrBasisSet->dPhi[k*OrbitalSetSize+j];
      }
      psi[j]=res;
      dpsi[j]=dres;
      d2psi[j]=d2res;
    }
  }


  /** evaluate everything for the walker move
   *
   * Using gemm can improve the performance for a larger problem
   */
//    void evaluate(const ParticleSet& P, int first, int last,
//        ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet) {
//      //unroll myself: rely on the compilers
//      //optimal on xeon
//#pragma ivdep
//      for(int i=0, iat=first; iat<last; i++,iat++){
//        myBasisSet->evaluateForWalkerMove(P,iat);
//        const ValueType* restrict cptr=C.data();
//        const typename BS::ValueType* restrict pptr=myBasisSet->Phi.data();
//        const typename BS::ValueType* restrict d2ptr=myBasisSet->d2Phi.data();
//        const typename BS::GradType* restrict dptr=myBasisSet->dPhi.data();
//        for(int j=0; j<OrbitalSetSize; j++) {
//          register ValueType res=0.0, d2res=0.0;
//          register GradType dres;
//          for(int b=0; b<BasisSetSize; b++,cptr++) {
//            res += *cptr*pptr[b];
//            d2res += *cptr*d2ptr[b];
//            dres += *cptr*dptr[b];
//          }
//          logdet(j,i)=res;dlogdet(i,j)=dres;d2logdet(i,j)=d2res;
//        }
//      }
//    }

  void evaluate_notranspose(const ParticleSet& P, int first, int last,
                            ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
  {
    const ValueType* restrict cptr=C.data();
    int numCenters=corrBasisSet->NumCenters;
#pragma ivdep
    for(int i=0,ij=0, iat=first; iat<last; i++,iat++)
    {
      myBasisSet->evaluateForWalkerMove(P,iat);
      corrBasisSet->evaluateForWalkerMove(P,iat);
      MatrixOperators::product(C,myBasisSet->Phi,logdet[i]);
      MatrixOperators::product(C,myBasisSet->d2Phi,d2logdet[i]);
      const typename BS::GradType* restrict dptr=myBasisSet->dPhi.data();
      for(int j=0,jk=0; j<OrbitalSetSize; j++)
      {
        register GradType dres;
        for(int k=0 ; k<numCenters; k++)
        {
          logdet(i,j) += corrBasisSet->Phi[k*OrbitalSetSize+j];
          d2logdet(i,j) += corrBasisSet->d2Phi[k*OrbitalSetSize+j];
          dres += corrBasisSet->dPhi[k*OrbitalSetSize+j];
        }
        for(int b=0; b<BasisSetSize; ++b)
        {
          dres +=  cptr[jk++]*dptr[b];
        }
        dlogdet(ij)=dres;
        ++ij;
      }
    }
  }

  void evaluate(const ParticleSet& P, int first, int last,
                ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
  {
    evaluate_notranspose(P,first,last,logdet,dlogdet,d2logdet);
    MatrixOperators::transpose(logdet);
  }

  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    const ValueType* restrict cptr=C.data();
    int numCenters=corrBasisSet->NumCenters;
#pragma ivdep
    for(int i=0,ij=0, iat=first; iat<last; i++,iat++)
    {
      myBasisSet->evaluateWithHessian(P,iat);
      corrBasisSet->evaluateWithHessian(P,iat);
      MatrixOperators::product(C,myBasisSet->Phi,logdet[i]);
      const typename BS::GradType* restrict dptr=myBasisSet->dPhi.data();
      const typename BS::HessType* restrict d2ptr=myBasisSet->grad_grad_Phi.data();
      for(int j=0,jk=0; j<OrbitalSetSize; j++)
      {
        register GradType dres;
        register HessType d2res;
        for(int k=0 ; k<numCenters; k++)
        {
          logdet(i,j) += corrBasisSet->Phi[k*OrbitalSetSize+j];
          d2res += corrBasisSet->grad_grad_Phi[k*OrbitalSetSize+j];
          dres += corrBasisSet->dPhi[k*OrbitalSetSize+j];
        }
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
    const ValueType* restrict cptr=C.data();
    int numCenters=corrBasisSet->NumCenters;
#pragma ivdep
    for(int i=0,ij=0, iat=first; iat<last; i++,iat++)
    {
      myBasisSet->evaluateWithThirdDeriv(P,iat);
      corrBasisSet->evaluateWithThirdDeriv(P,iat);
      MatrixOperators::product(C,myBasisSet->Phi,logdet[i]);
      const typename BS::GradType* restrict dptr=myBasisSet->dPhi.data();
      const typename BS::HessType* restrict d2ptr=myBasisSet->grad_grad_Phi.data();
      const typename BS::GGGType* restrict gggptr=myBasisSet->grad_grad_grad_Phi.data();
      for(int j=0,jk=0; j<OrbitalSetSize; j++)
      {
        register GradType dres;
        register HessType d2res;
        register GGGType gggres;
        for(int k=0 ; k<numCenters; k++)
        {
          logdet(i,j) += corrBasisSet->Phi[k*OrbitalSetSize+j];
          d2res += corrBasisSet->grad_grad_Phi[k*OrbitalSetSize+j];
          dres += corrBasisSet->dPhi[k*OrbitalSetSize+j];
          gggres[0] +=  (corrBasisSet->grad_grad_grad_Phi[k*OrbitalSetSize+j])[0];
          gggres[1] +=  (corrBasisSet->grad_grad_grad_Phi[k*OrbitalSetSize+j])[1];
          gggres[2] +=  (corrBasisSet->grad_grad_grad_Phi[k*OrbitalSetSize+j])[2];
        }
        for(int b=0; b<BasisSetSize; ++b,++jk)
        {
          dres +=  cptr[jk]*dptr[b];
          d2res +=  cptr[jk]*d2ptr[b];
          gggres[0] +=  cptr[jk]*(gggptr[b])[0];
          gggres[1] +=  cptr[jk]*(gggptr[b])[1];
          gggres[2] +=  cptr[jk]*(gggptr[b])[2];
        }
        dlogdet(ij)=dres;
        grad_grad_logdet(ij)=d2res;
        grad_grad_grad_logdet(ij)=gggres;
        ++ij;
      }
    }
  }

  void evaluateThirdDeriv(const ParticleSet& P, int first, int last
                          , GGGMatrix_t& grad_grad_grad_logdet)
  {
    const ValueType* restrict cptr=C.data();
    int numCenters=corrBasisSet->NumCenters;
#pragma ivdep
    for(int i=0,ij=0, iat=first; iat<last; i++,iat++)
    {
      myBasisSet->evaluateThirdDerivOnly(P,iat);
      corrBasisSet->evaluateThirdDerivOnly(P,iat);
      const typename BS::GGGType* restrict gggptr=myBasisSet->grad_grad_grad_Phi.data();
      for(int j=0,jk=0; j<OrbitalSetSize; j++)
      {
        register GGGType gggres;
        for(int k=0 ; k<numCenters; k++)
        {
          gggres[0] +=  (corrBasisSet->grad_grad_grad_Phi[k*OrbitalSetSize+j])[0];
          gggres[1] +=  (corrBasisSet->grad_grad_grad_Phi[k*OrbitalSetSize+j])[1];
          gggres[2] +=  (corrBasisSet->grad_grad_grad_Phi[k*OrbitalSetSize+j])[2];
        }
        for(int b=0; b<BasisSetSize; ++b,++jk)
        {
          gggres[0] +=  cptr[jk]*(gggptr[b])[0];
          gggres[1] +=  cptr[jk]*(gggptr[b])[1];
          gggres[2] +=  cptr[jk]*(gggptr[b])[2];
        }
        grad_grad_grad_logdet(ij)=gggres;
        ++ij;
      }
    }
  }


  LCOrbitalSet<BS,false>* clone2LCOrbitalSet();

  bool transformSPOSet();

  bool readCuspInfo(Matrix<TinyVector<RealType,9> > &);

private:

  BS* extractHighYLM(vector<bool> &rmv);
  LCOrbitalSet<BS,false>* originalSPOSet;
  double Rcut;
  vector<double> Z;
  string cuspInfoFile;

  void createLCOSets(int centr, LCOrbitalSet<BS,false>* Phi, LCOrbitalSet<BS,false>* Eta);

};

}

#include "QMCWaveFunctions/Experimental/LCOrbitalSetWithCorrection.cpp"

#endif

/***************************************************************************
 * $RCSfile$   $Author: jeongnim.kim $
 * $Revision: 4351 $   $Date: 2009-11-01 15:40:52 -0600 (Sun, 01 Nov 2009) $
 * $Id: LCOrbitalSetWithCorrection.h 4351 2010-01-11 21:40:52Z jeongnim.kim $
 ***************************************************************************/

