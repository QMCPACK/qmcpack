//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include "QMCWaveFunctions/AFMSPOBuilder.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

AFMSPOBuilder::AFMSPOBuilder
(ParticleSet& p, PtclPoolType& psets, xmlNodePtr cur) :
  targetPtcl (&p)
{
}

bool
AFMSPOBuilder::put (xmlNodePtr cur)
{
  return true;
}


SPOSetBase*
AFMSPOBuilder::createSPOSetFromXML(xmlNodePtr cur)
{
  AFMSPOSet *spo =  new AFMSPOSet();
  return spo;
}

}
