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
#ifndef QMCPLUSPLUS_NUMERICAL_PBC_JASTROW_AA_BUILDER_H
#define QMCPLUSPLUS_NUMERICAL_PBC_JASTROW_AA_BUILDER_H
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/Jastrow/NumericalJastrowFunctor.h"

namespace qmcplusplus {

  //forward declaration
  class ParticleSet;

  /**@ingroup WFSBuilder
   *A builder class to add a numerical two-body Jastrow function to a TrialWaveFunction
   *
   *A xml node with OrbtialBuilderBase::jastrow_tag is parsed recursively.
   */
  struct JAAPBCBuilder: public OrbitalBuilderBase {

    JAAPBCBuilder(ParticleSet& p, TrialWaveFunction& psi);

    /**@param cur the current xmlNodePtr to be processed by NumericalJastrowBuilder
     *@return true if succesful
     */
    bool put(xmlNodePtr cur);

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
