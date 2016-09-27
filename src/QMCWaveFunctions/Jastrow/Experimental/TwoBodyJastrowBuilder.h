//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//                    Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_TWOBODYJASTROWCOMBO_BUILDER_H
#define QMCPLUSPLUS_TWOBODYJASTROWCOMBO_BUILDER_H
#include "QMCWaveFunctions/OrbitalBuilderBase.h"

namespace qmcplusplus
{

/** TwoBodyJastrow Jastrow Builder with constraints
 */
class TwoBodyJastrowBuilder: public OrbitalBuilderBase
{

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
