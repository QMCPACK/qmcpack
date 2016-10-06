//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include "QMCWaveFunctions/OptimizableSPOBuilder.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

OptimizableSPOBuilder::OptimizableSPOBuilder
(ParticleSet& p, PtclPoolType& psets, xmlNodePtr cur) :
  targetPtcl (&p)
{
}

bool
OptimizableSPOBuilder::put (xmlNodePtr cur)
{
  return true;
}


SPOSetBase*
OptimizableSPOBuilder::createSPOSetFromXML(xmlNodePtr cur)
{
  OptimizableSPOSet *spo =  new OptimizableSPOSet();
  return spo;
}

}
