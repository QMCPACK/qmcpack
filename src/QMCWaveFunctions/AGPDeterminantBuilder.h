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
/**@file AGPDeterminantBuilder.h
 *@brief declaration of a builder class for AGPDeterminant
 */
#ifndef QMCPLUSPLUS_AGPDETERMINANT_GEMINALBUILDER_H
#define QMCPLUSPLUS_AGPDETERMINANT_GEMINALBUILDER_H 
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/BasisSetFactory.h"
namespace qmcplusplus {

  class AGPDeterminant;

  /**@ingroup WFSBuilder
   * @brief An abstract class for wave function builders
   */
  class AGPDeterminantBuilder: public OrbitalBuilderBase  
  {
    
  public:

    AGPDeterminantBuilder(ParticleSet& els, TrialWaveFunction& wfs, PtclPoolType& pset);

    /// process a xml node at cur
    bool put(xmlNodePtr cur);

  protected:

    ///reference to a PtclPoolType
    PtclPoolType& ptclPool;
    ///basiset Factory
    BasisSetFactory* myBasisSetFactory;
    ///AGPDeterminant
    AGPDeterminant* agpDet;
    string funcOpt;
    string transformOpt;

    template <typename BasisBuilderT> 
    bool createAGP(BasisBuilderT* abuilder, xmlNodePtr cur);
    
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
