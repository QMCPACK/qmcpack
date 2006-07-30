//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
/** @file BasisSetBase.h
 * @brief Declaration of a base class of BasisSet
 */
#ifndef QMCPLUSPLUS_BASISSETBASE_H
#define QMCPLUSPLUS_BASISSETBASE_H

#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "Particle/ParticleSet.h"

namespace qmcplusplus {

  /** base class for a basis set
   *
   * BasisSetBase can function as a special case of LCOrbitalSet where the coefficient is an 
   * Identity matrix.
   */
  struct BasisSetBase: public OrbitalSetTraits {

    ///size of the basis set
    IndexType BasisSetSize;
    ///phi[i] the value of the i-th basis set 
    ValueVector_t Phi;
    ///dphi[i] the gradient of the i-th basis set 
    GradVector_t  dPhi;
    ///d2phi[i] the laplacian of the i-th basis set 
    ValueVector_t d2Phi;
    ///container to store value, laplacian and gradient
    ValueMatrix_t Temp;

    ///default constructor
    BasisSetBase();
    ///virtual destructor
    virtual ~BasisSetBase();

    /** resize the container */
    void resize();

    /** return the basis set size */
    inline IndexType getBasisSetSize() const {
      return BasisSetSize;
    }

    /** fill in Phi
     * @param iat index of the particle for \f$\Phi_{i} (r_{iat})\f$
     */
    void evaluate(const ParticleSet& P, int iat) {
      evaluateBasis(P,iat,Phi.data());
    }

    /** fill in Phi, dPhi and d2Phi
     * @param iat index of the particle for \f$\Phi_{i} (r_{iat})\f$
     *
     * Evaluate the gradients and laplacians are evaluated as well.
     */
    void evaluateAll(const ParticleSet& P, int iat) {
      evaluateBasis(P,iat,Phi.data(),dPhi.data(),d2Phi.data());
    }

    ///reset the basis set
    virtual void reset() = 0;
    ///resize the basis set
    virtual void setBasisSetSize(int nbs) = 0;
    ///reset the target particle set
    virtual void resetTargetParticleSet(ParticleSet& P)=0;
    ///evaluate the value of basis functions for the iath-particle
    virtual void evaluateBasis(const ParticleSet& P, int iat, ValueType* restrict psi)=0;
    ///evaluate the value, gradient and laplacian of basis functions for the iath-particle
    virtual void evaluateBasis(const ParticleSet& P, int iat, ValueType* restrict psi, 
        GradType* restrict dpsi, ValueType* restrict d2psi)=0;
    ///evaluate the value, gradient and laplacian of basis functions for the iath-particle
    virtual void evaluateBasis(const ParticleSet& P, int iat, ValueMatrix_t& temp)=0;

  };

  /** base class for the BasisSet builder
   */
  struct BasisSetBuilder {
    BasisSetBase* myBasisSet;
    BasisSetBuilder():myBasisSet(0) {}
    virtual ~BasisSetBuilder(){}
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
