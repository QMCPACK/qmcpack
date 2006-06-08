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
/**@file ThreeBodyPadeBuilder.h
 *@brief declaration of three-body jastrow of Pade functions
 */
#ifndef QMCPLUSPLUS_THREEBODY_PADE_BUILDER_H
#define QMCPLUSPLUS_THREEBODY_PADE_BUILDER_H 
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
namespace qmcplusplus {

  //class GTOMolecularOrbitals;
  class ThreeBodyPade;
 // class GridMolecularOrbitals;

  /**@ingroup WFSBuilder
   * @brief An abstract class for wave function builders
   */
  class ThreeBodyPadeBuilder: public OrbitalBuilderBase  {
    
  public:

    ThreeBodyPadeBuilder(ParticleSet& els, TrialWaveFunction& wfs,
      ParticleSet& ions);

    /// process a xml node at cur
    bool put(xmlNodePtr cur);

  protected:

//    GridMolecularOrbitals* gtoBuilder;
    ThreeBodyPade* J3;
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
