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

  class WalkerSetRef;
  class DistanceTableData;

  /** An abstract class for Local Energy operators 
   *
   * Use of ValueType is questionable. The types should be checked when using
   * complex wave functions.
   */ 
  struct QMCHamiltonianBase: public QMCTraits {

    RealType Tau;
    RealType Value;

    typedef ParticleAttrib<ValueType>  ValueVectorType;

    ///constructor
    QMCHamiltonianBase():Tau(0.0),Value(0.0){}

    ///virtual destructor
    virtual ~QMCHamiltonianBase() { }

    /** Evaluate the local energies of an N-particle configuration
     *@param P input configuration containing N particles
     *@return the value of the Hamiltonian
     */
    virtual ValueType evaluate(ParticleSet& P) = 0; 

    /** Evaluate the local energies of an N-particle configuration
     *@param P input configuration containing N particles
     *@param x the sum of local energies
     *@return the value the Local Energy
    */
    virtual ValueType evaluate(ParticleSet& P, RealType& x) = 0;

    /** Evaluate the local energies of the entire walkers
     *@param W a set of walkers (N-particle configurations)
     *@param LE return a vector containing the value
     */
    virtual 
    void evaluate(WalkerSetRef& W, ValueVectorType& LE) = 0;

    inline void setTau(RealType tau) { Tau = tau;}
  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

