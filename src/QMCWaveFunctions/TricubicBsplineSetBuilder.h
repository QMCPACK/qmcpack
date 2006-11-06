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
#ifndef QMCPLUSPLUS_TRICUBIC_BSPLINESETBUILDER_H
#define QMCPLUSPLUS_TRICUBIC_BSPLINESETBUILDER_H

#include "QMCWaveFunctions/BasisSetBase.h"
#include "QMCWaveFunctions/GroupedOrbitalSet.h"
#include "Numerics/TricubicBsplineSet.h"

namespace qmcplusplus {

  /**@ingroup WFSBuilder
   * A builder class for a set of Spline functions
   */
  class TricubicBsplineSetBuilder: public BasisSetBuilder {

  public:

    typedef TricubicBsplineSet<ValueType>              OrbitalGroupType;      
    typedef TricubicBsplineSet<ValueType>::StorageType StorageType;
    typedef GroupedOrbitalSet<OrbitalGroupType>        SPOSetType;             
    typedef map<string,ParticleSet*>                   PtclPoolType;

    /** constructor
     * @param p target ParticleSet
     * @param psets a set of ParticleSet objects
     */
    TricubicBsplineSetBuilder(ParticleSet& p, PtclPoolType& psets);

    bool put(xmlNodePtr cur);

    /** initialize the Antisymmetric wave function for electrons
     *@param cur the current xml node
     */
    SPOSetBase* createSPOSet(xmlNodePtr cur);

  private:
    bool DebugWithEG;
    ///target ParticleSet
    ParticleSet& targetPtcl;
    ///reference to a ParticleSetPool
    PtclPoolType& ptclPool;

    PosType LowerBox;
    PosType UpperBox;
    TinyVector<IndexType,DIM> BoxGrid;
    ///three-dimnesional grid
    //GridType* GridXYZ;
    ///set of StorageType*
    map<string,StorageType*> BigDataSet;
    ///set of WFSetType*
    map<string,OrbitalGroupType*> myBasis;
    ///single-particle orbital sets
    map<string,SPOSetType*> mySPOSet;

    SPOSetBase* createSPOSetWithEG();
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
