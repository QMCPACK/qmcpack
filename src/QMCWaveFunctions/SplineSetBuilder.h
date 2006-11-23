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
#ifndef QMCPLUSPLUS_SPLINESETBUILDER_H
#define QMCPLUSPLUS_SPLINESETBUILDER_H

#include "Numerics/TriCubicSplineT.h"
#include "QMCWaveFunctions/BasisSetBase.h"
#include "QMCWaveFunctions/SingleParticleOrbitalSet.h"

namespace qmcplusplus {

  /**@ingroup WFSBuilder
   * A builder class for a set of Spline functions
   */
  class SplineSetBuilder: public BasisSetBuilder {

  public:

    typedef TriCubicSplineT<ValueType,RealType>           SPOType;
    typedef TriCubicSplineT<ValueType,RealType>::GridType GridType;
    typedef SingleParticleOrbitalSet<SPOType>             SPOSetType;
    typedef map<string,ParticleSet*> PtclPoolType;

    /** constructor
     * @param p target ParticleSet
     * @param psets a set of ParticleSet objects
     */
    SplineSetBuilder(ParticleSet& p, PtclPoolType& psets);

    bool put(xmlNodePtr cur);

    /** initialize the Antisymmetric wave function for electrons
     *@param cur the current xml node
     */
    SPOSetBase* createSPOSet(xmlNodePtr cur);

  private:
    ///target ParticleSet
    ParticleSet& targetPtcl;
    ///reference to a ParticleSetPool
    PtclPoolType& ptclPool;
    ///global GridType*, should be generalized for any number of grids
    GridType* GridXYZ;
    ///set of SPOType*
    map<string,SPOType* > NumericalOrbitals;
    ///set of SPOSetType*
    map<string,SPOSetType*> SPOSet;
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
