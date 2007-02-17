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
#ifndef QMCPLUSPLUS_SINGLEPARTICLEORBITALSETBASE_H
#define QMCPLUSPLUS_SINGLEPARTICLEORBITALSETBASE_H

#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"
//#define ENABLE_SMARTPOINTER
#if defined(ENABLE_SMARTPOINTER)
#include <boost/shared_ptr.hpp>
#endif

namespace qmcplusplus {

  /** base class for Single-particle orbital sets
   *
   * SPOSetBase stands for S(ingle)P(article)O(rbital)SetBase which contains
   * a number of single-particle orbitals with capabilities of evaluating \f$ \psi_j({\bf r}_i)\f$ 
   */
  class SPOSetBase: public OrbitalSetTraits {

  public:

    ///true if C is an identity matrix
    bool Identity;
    ///number of Single-particle orbtials
    IndexType OrbitalSetSize;
    ///number of Single-particle orbtials
    int BasisSetSize;
    ///matrix containing the coefficients
    ValueMatrix_t C;
    ///occupation number
    vector<RealType> Occ;

    /** constructor
     */
    SPOSetBase():Identity(false),OrbitalSetSize(0),BasisSetSize(0) {}

    /** destructor
     */
    virtual ~SPOSetBase() {}

    /** return the size of the orbital set
     */
    inline int getOrbitalSetSize() const { 
      return OrbitalSetSize;
    }

    inline int getBasisSetSize() const { 
      return BasisSetSize;
    }

    bool setIdentity(bool useIdentity) {
      Identity=useIdentity;
      C.resize(OrbitalSetSize,BasisSetSize);
      for(int i=0; i<OrbitalSetSize; i++) C(i,i)=1.0;
      return true;
    }

    void checkObject();

    ///get C and Occ
    bool put(xmlNodePtr cur);

    ///reset
    virtual void resetParameters(OptimizableSetType& optVariables)=0;
    ///reset the target particleset
    virtual void resetTargetParticleSet(ParticleSet& P)=0;
    /** set the OrbitalSetSize
     * @param norbs number of single-particle orbitals
     */
    virtual void setOrbitalSetSize(int norbs)=0;

    /** evaluate the values of this single-particle orbital set
     * @param P current ParticleSet
     * @param iat active particle
     * @param psi values of the SPO
     */
    virtual void 
    evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)=0;

    /** evaluate the values, gradients and laplacians of this single-particle orbital set
     * @param P current ParticleSet
     * @param iat active particle
     * @param psi values of the SPO
     */
    virtual void 
    evaluate(const ParticleSet& P, int iat, 
        ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)=0;

    /** evaluate the values, gradients and laplacians of this single-particle orbital for [first,last)particles
     * @param P current ParticleSet
     * @param first starting index of the particles
     * @param last ending index of the particles
     * @param logdet determinant matrix to be inverted
     * @param dlogdet gradients
     * @param d2logdet laplacians
     */
    virtual void evaluate(const ParticleSet& P, int first, int last,
        ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)=0;

protected:
    bool putOccupation(xmlNodePtr occ_ptr);
    bool putFromXML(xmlNodePtr coeff_ptr);
    bool putFromH5(const char* fname, xmlNodePtr coeff_ptr);
  };

#if defined(ENABLE_SMARTPOINTER)
  typedef boost::shared_ptr<SPOSetBase> SPOSetBasePtr;
#else
  typedef SPOSetBase*                   SPOSetBasePtr;
#endif

}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

