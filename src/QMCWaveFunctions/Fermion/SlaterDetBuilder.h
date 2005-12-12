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
#ifndef QMCPLUSPLUS_SLATERDETSET_BUILDER_WITH_BASISSET_H
#define QMCPLUSPLUS_SLATERDETSET_BUILDER_WITH_BASISSET_H

#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"

namespace qmcplusplus {

  struct SlaterDetBuilder: public OrbitalBuilderBase {
  
    typedef DiracDeterminantBase Det_t;
    typedef SlaterDet            SlaterDeterminant_t;

    /** index to map SPO and builder */
    int SPOSetID;
    /** reference to OrbitalBuilderBase:PtclPoolType */
    PtclPoolType& ptclPool;
    /** pointer to a SPO Builder */
    OrbitalBuilderBase* spoBuilder;
    /** pointer to a BasisSetBase 
     */
    BasisSetBase* BasisSet;

    /** name of the basis set*/
    string basisName;

    /** constructor
     *@param p particleset whose positions defines the wave function
     *@param psi trial wavefuntion to which determinant terms are added
     *@param abuilder a BasisBuilderT object, provides addBasisSet and typedefs
     */
    SlaterDetBuilder(ParticleSet& p, TrialWaveFunction& psi,
        PtclPoolType& psets);
 
    /** process the current xml node to create single-particle orbital
     *@param cur xmlNodePtr to be processed
     *@return true when the size of determinants is positive.
     */
    bool put(xmlNodePtr cur);


    /** create a BasisSet
     */
    void createBasisSet(xmlNodePtr cur);

    /** create a Determinant
     */
    Det_t* createDeterminant(xmlNodePtr cur, int is, int first);

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
