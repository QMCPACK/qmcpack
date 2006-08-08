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

  class SPOSetBase;

  /** base class for a basis set
   *
   * Define a common storage for the derived classes and 
   * provides  a minimal set of interfaces to get/set BasisSetSize.
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

    ValueMatrix_t Y;
    GradMatrix_t dY;
    ValueMatrix_t d2Y;

    ///default constructor
    BasisSetBase();
    ///virtual destructor
    virtual ~BasisSetBase();

    /** resize the container */
    void resize(int ntargets);

    /** return the basis set size */
    inline IndexType getBasisSetSize() const {
      return BasisSetSize;
    }

    ///reset the basis set
    virtual void reset() = 0;
    ///resize the basis set
    virtual void setBasisSetSize(int nbs) = 0;
    ///reset the target particle set
    virtual void resetTargetParticleSet(ParticleSet& P)=0;
  };


  /** base class for the BasisSet builder
   */
  struct BasisSetBuilder {
    BasisSetBase* myBasisSet;
    BasisSetBuilder():myBasisSet(0) {}
    virtual ~BasisSetBuilder(){}
    virtual bool put(xmlNodePtr cur)=0;
    virtual SPOSetBase* createSPOSet(xmlNodePtr cur)=0;
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
