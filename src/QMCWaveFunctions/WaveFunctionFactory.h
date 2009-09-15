/////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
/**@file WaveFunctionFactory.h
 *@brief Declaration of a WaveFunctionFactory 
 */
#ifndef QMCPLUSPLUS_TRIALWAVEFUNCTION_FACTORY_H
#define QMCPLUSPLUS_TRIALWAVEFUNCTION_FACTORY_H

#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "Message/MPIObjectBase.h"
namespace qmcplusplus {

  /** Factory class to build a many-body wavefunction 
   */
  struct WaveFunctionFactory: public MPIObjectBase
  {

    typedef map<string,ParticleSet*> PtclPoolType;
    ///target ParticleSet
    ParticleSet* targetPtcl;
    ///many-body wavefunction object
    TrialWaveFunction* targetPsi;
    ///reference to the PtclPoolType
    PtclPoolType&  ptclPool;
    ///input node for a many-body wavefunction
    xmlNodePtr myNode;
    ///builder tree
    std::vector<OrbitalBuilderBase*> psiBuilder;

    /** constructor
     * @param qp quantum particleset
     * @param pset pool of particlesets
     * @param c  communicator
     */
    WaveFunctionFactory(ParticleSet* qp, PtclPoolType& pset, Communicate* c);

    ~WaveFunctionFactory();

    void setPsi(TrialWaveFunction* psi);

    ///read from xmlNode
    bool put(xmlNodePtr cur);

    ///reset member data
    void reset();

    /** process xmlNode to populate targetPsi
     */
    bool build(xmlNodePtr cur, bool buildtree=true);

    /** add Fermion wavefunction term */
    bool addFermionTerm(xmlNodePtr cur);

    /** add an OrbitalBuilder and the matching xml node
     * @param b OrbitalBuilderBase*
     * @oaram cur xmlNode for b
     * @return true if successful
     */
    bool addNode(OrbitalBuilderBase* b, xmlNodePtr cur);

    void setCloneSize(int np);

    WaveFunctionFactory* clone(ParticleSet* qp, int ip, const string& aname);

    vector<WaveFunctionFactory*> myClones;
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
