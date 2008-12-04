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
#ifndef QMCPLUSPLUS_ORIGINAL_JASTROW_AA_BUILDER_H
#define QMCPLUSPLUS_ORIGINAL_JASTROW_AA_BUILDER_H
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"

namespace qmcplusplus {

  //forward declaration
  class ParticleSet;

  /** Generic two-body Jastrow builder
   * 
   * Replacement of JastrowBuilder::createTwoBodySpin and JastrowBuilder::createTwoBodyNoSpin
   */
  struct JAABuilder: public OrbitalBuilderBase {

    JAABuilder(ParticleSet& p, TrialWaveFunction& psi);

    bool put(xmlNodePtr cur);

    template <class FN> TwoBodyJastrowOrbital<FN>* createJAA(xmlNodePtr cur, const string& jname);

    bool IgnoreSpin;

    DistanceTableData* d_table;
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
