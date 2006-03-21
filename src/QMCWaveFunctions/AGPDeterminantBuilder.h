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
#ifndef QMCPLUSPLUS_AGPDETERMINANT_GEMINALBUILDER_H
#define QMCPLUSPLUS_AGPDETERMINANT_GEMINALBUILDER_H 
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
namespace qmcplusplus {

  class AGPDeterminant;

  /**@ingroup WFSBuilder
   * @brief An abstract class for wave function builders
   */
  class AGPDeterminantBuilder: public OrbitalBuilderBase  {
    
  public:

    AGPDeterminantBuilder(ParticleSet& els, TrialWaveFunction& wfs,
      ParticleSet& ions);

    /// process a xml node at cur
    bool put(xmlNodePtr cur);

  protected:

    OrbitalBuilderBase* basisBuilder;
    AGPDeterminant* agpDet;
    ParticleSet& ionRef;

    template <class BasisBuilderT> 
    bool createAGP(BasisBuilderT* abuilder, xmlNodePtr cur);
    
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
