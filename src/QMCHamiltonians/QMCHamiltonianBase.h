//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
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
/**@file QMCHamiltonianBase.h
 *@brief Declaration of QMCHamiltonianBase
 */
#ifndef OHMMS_QMC_HAMILTONIANBASE_H
#define OHMMS_QMC_HAMILTONIANBASE_H

#include "Particle/ParticleSet.h"
#include "Utilities/PooledData.h"

namespace ohmmsqmc {

  class DistanceTableData;

  /** An abstract class for Local Energy operators 
   *
   * Return_t is defined as RealTye. 
   * The types should be checked when using complex wave functions.
   */ 
  struct QMCHamiltonianBase: public QMCTraits {
    
    /**\typedef Return value of the functions
     */
    typedef RealType Return_t;
 
    RealType Tau;
    RealType Value;
   
    ///constructor
    QMCHamiltonianBase():Tau(0.0),Value(0.0){}

    ///virtual destructor
    virtual ~QMCHamiltonianBase() { }

    /** Evaluate the local energies of an N-particle configuration
     *@param P input configuration containing N particles
     *@return the value of the Hamiltonian
     */
    virtual Return_t evaluate(ParticleSet& P) = 0; 

    /** Evaluate the local energies of an N-particle configuration
     *@param P input configuration containing N particles
     *@param x the sum of local energies
     *@return the value the Local Energy
    */
    virtual Return_t evaluate(ParticleSet& P, RealType& x) = 0;
    
    /** read the input parameter
     * @param cur xml node for a QMCHamiltonianBase object
     */
    virtual bool put(xmlNodePtr cur)=0;
    
    /** set Tau
     * @param tau time step
     */
    inline void setTau(RealType tau) { Tau = tau;}
  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

