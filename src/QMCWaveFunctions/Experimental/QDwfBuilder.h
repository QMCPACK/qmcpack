//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_QDWFBUILDER_H
#define QMCPLUSPLUS_QDWFBUILDER_H

#include "OhmmsData/OhmmsElementBase.h"
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/SingleParticleOrbitalSet.h"
#include "QMCWaveFunctions/QDwf.h"

namespace qmcplusplus
{

class QDwfBuilder: public OrbitalBuilderBase
{

  typedef SingleParticleOrbitalSet<QDwf> SPOSet_t;

public :

  QDwfBuilder(TrialWaveFunction& a) : OrbitalBuilderBase(a) {}

  bool put(xmlNodePtr cur);

};

}
#endif
