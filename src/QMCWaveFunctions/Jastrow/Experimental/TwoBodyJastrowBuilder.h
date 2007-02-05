//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_TWOBODYJASTROWCOMBO_BUILDER_H
#define QMCPLUSPLUS_TWOBODYJASTROWCOMBO_BUILDER_H
#include "QMCWaveFunctions/OrbitalBuilderBase.h"

namespace qmcplusplus {

  /** TwoBodyJastrow Jastrow Builder with constraints
   */
  class TwoBodyJastrowBuilder: public OrbitalBuilderBase {

    public:

    TwoBodyJastrowBuilder(ParticleSet& p, TrialWaveFunction& psi, PtclPoolType& psets);

    bool put(xmlNodePtr cur);

    private:
    PtclPoolType& ptclPool;
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
