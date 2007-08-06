//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Ken Esler and Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////

#include "QMCWaveFunctions/GroupedOrbitalSet.h"
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "OhmmsData/AttributeSet.h"
#include "Utilities/Timer.h"
#include "Message/Communicate.h"

namespace qmcplusplus {
  EinsplineSetBuilder::EinsplineSetBuilder(ParticleSet& p, PtclPoolType& psets, xmlNodePtr cur) 
  {
  }

  EinsplineSetBuilder::~EinsplineSetBuilder()
  {

  }


  bool 
  EinsplineSetBuilder::put(xmlNodePtr cur) {
    string hdfName;
    OhmmsAttributeSet attribs;
    attribs.add (hdfName, "href");
    attribs.put(cur);
    cerr << "In EinsplineSetBuilder::put.  href=\"" << hdfName << "\".\n";
  }

  SPOSetBase*
  EinsplineSetBuilder::createSPOSet(xmlNodePtr cur) {
    string hdfName;
    OhmmsAttributeSet attribs;
    attribs.add (hdfName, "href");
    attribs.put (cur);
    cerr << "In EinsplineSetBuilder::createSPOSet().  href=\"" << hdfName << "\".\n";
  }
}
