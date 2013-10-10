//////////////////////////////////////////////////////////////////
// (c) Copyright 2010-  by Ken Esler and Jeongnim Kim
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: esler@uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////

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
