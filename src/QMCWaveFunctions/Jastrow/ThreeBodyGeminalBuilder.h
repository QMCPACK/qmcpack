//////////////////////////////////////////////////////////////////
// (c) Copyright 2005-  by Jeongnim Kim
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
/**@file ThreeBodyGeminalBuilder.h
 *@brief declaration of three-body jastrow of Geminal functions
 */
#ifndef QMCPLUSPLUS_THREEBODY_GEMINALBUILDER_H
#define QMCPLUSPLUS_THREEBODY_GEMINALBUILDER_H 
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
namespace qmcplusplus {

  class ThreeBodyGeminal;
  class BasisSetBuilder;

  /**@ingroup WFSBuilder
   * @brief An abstract class for wave function builders
   */
  class ThreeBodyGeminalBuilder: public OrbitalBuilderBase  
  {
    
  public:

    ThreeBodyGeminalBuilder(ParticleSet& els, TrialWaveFunction& wfs,
      ParticleSet& ions);

    /// process a xml node at cur
    bool put(xmlNodePtr cur);

  protected:

    ParticleSet& sourcePtcl;
    BasisSetBuilder* basisBuilder;
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
