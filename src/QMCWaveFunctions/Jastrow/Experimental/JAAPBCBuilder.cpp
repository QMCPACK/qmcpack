////////////////////////////////////////////////////////////////// // (c) Copyright 2003-  by Jeongnim Kim
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
#include "Utilities/OhmmsInfo.h"
#include "QMCWaveFunctions/Jastrow/JAAPBCBuilder.h"
#include "QMCWaveFunctions/Jastrow/LRTwoBodyJastrow.h"
#include "LongRange/LRJastrowSingleton.h"

namespace qmcplusplus {

  /** constructor
   * @param p target ParticleSet whose wave function is to be initialized
   *@param psi the wavefunction
   *@param psets a vector containing the ParticleSets
   * 
   *Jastrow wavefunctions chose the needed distance tables and the 
   *DistanceTableData objects are initialized based on the source 
   *and target particle sets.
   */
  JAAPBCBuilder::JAAPBCBuilder(ParticleSet& p, TrialWaveFunction& psi):
    OrbitalBuilderBase(p,psi)
  { }

  bool JAAPBCBuilder::put(xmlNodePtr cur) {
    typedef LRJastrowSingleton::LRHandlerType HandlerType;
    HandlerType* handler = LRJastrowSingleton::getHandler(targetPtcl);
    LRTwoBodyJastrow *J2 = new LRTwoBodyJastrow(targetPtcl, handler);
    bool success = J2->put(cur, targetPsi.VarList);
    J2->setOptimizable(true);
    targetPsi.addOrbital(J2);
    return success;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
