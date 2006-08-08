//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_LINEARCOMIBINATIONORBITALSET_TEMP_H
#define QMCPLUSPLUS_LINEARCOMIBINATIONORBITALSET_TEMP_H

#include "QMCWaveFunctions/SPOSetBase.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/MatrixOperators.h"

namespace qmcplusplus {

  /** class to handle linear combinations of basis orbitals used to evaluate the Dirac determinants.
   *
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
  class LCOrbitalSet: public SPOSetBase {

  public:

    ///pointer to the basis set
    BS* myBasisSet;
    /** constructor
     * @param bs pointer to the BasisSet
     * @param id identifier of this LCOrbitalSet
     */
    LCOrbitalSet(BS* bs=0): myBasisSet(0) {
      if(bs) setBasisSet(bs);
    }

    /** destructor
     *
     * BasisSet is deleted by the object with ID == 0
     */
    ~LCOrbitalSet() {}

    ///reset
    void reset() {
      myBasisSet->reset();
    }

    ///reset the target particleset
    void resetTargetParticleSet(ParticleSet& P) {
      myBasisSet->resetTargetParticleSet(P);
    }

    /** set the OrbitalSetSize
     */
    void setOrbitalSetSize(int norbs) {
      OrbitalSetSize=norbs;
    }

    /** set the basis set
     */
    void setBasisSet(BS* bs) {
      myBasisSet=bs;
      BasisSetSize=myBasisSet->getBasisSetSize();
    }

    /** return the size of the basis set
     */
    inline int getBasisSetSize() const { 
      return (myBasisSet==0)? 0: myBasisSet->getBasisSetSize();
    }

    inline void 
    evaluate(const ParticleSet& P, int iat, ValueVector_t& psi) {
      myBasisSet->evaluate(P,iat);
      for(int j=0 ; j<OrbitalSetSize; j++) 
        psi[j] = dot(C[j],myBasisSet->Phi.data(),BasisSetSize);
        //psi[j] = BLAS::dot(BasisSetSize,C[j],myBasisSet->y(0));
      //MatrixOperators::product(C,myBasisSet->Phi,psi.data());
    }

    inline void 
    evaluate(const ParticleSet& P, int iat, 
        ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi) {
      myBasisSet->evaluateAll(P,iat);
      for(int j=0 ; j<OrbitalSetSize; j++) {
        psi[j]   = dot(C[j],myBasisSet->Phi.data(),  BasisSetSize);
        dpsi[j]  = dot(C[j],myBasisSet->dPhi.data(), BasisSetSize);
        d2psi[j] = dot(C[j],myBasisSet->d2Phi.data(),BasisSetSize);
        //psi[j]   = dot(C[j],myBasisSet->Y[0],  BasisSetSize);
        //dpsi[j]  = dot(C[j],myBasisSet->dY[0], BasisSetSize);
        //d2psi[j] = dot(C[j],myBasisSet->d2Y[0],BasisSetSize);
      }
    }

    void evaluate(const ParticleSet& P, int first, int last,
        ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet) {
      //if(first ==0) myBasisSet->evaluate(P);
      //for(int i=0, iat=first; iat<last; i++,iat++){
      //  for(int j=0 ; j<OrbitalSetSize; j++) {
      //    //logdet(j,i) = \f$\sum_k^{nb} C(j,k) Y(i+first,k)\f$
      //    logdet(j,i)   = dot(C[j],myBasisSet->Y[iat],  BasisSetSize);
      //    dlogdet(i,j)  = dot(C[j],myBasisSet->dY[iat], BasisSetSize);
      //    d2logdet(i,j) = dot(C[j],myBasisSet->d2Y[iat],BasisSetSize);
      //  }
      //}
     
      for(int i=0, iat=first; iat<last; i++,iat++){
        myBasisSet->evaluateBasis(P,iat);
        for(int j=0 ; j<OrbitalSetSize; j++) {
          //logdet(j,i) = \f$\sum_k^{nb} C(j,k) Y(i+first,k)\f$
          logdet(j,i)   = dot(C[j],myBasisSet->Phi.data(),  BasisSetSize);
          dlogdet(i,j)  = dot(C[j],myBasisSet->dPhi.data(), BasisSetSize);
          d2logdet(i,j) = dot(C[j],myBasisSet->d2Phi.data(),BasisSetSize);
        }
      }
    }
  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

