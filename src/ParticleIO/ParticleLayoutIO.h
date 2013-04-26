//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_PARTICLELAYOUT_INPUTOUTPUT_UTILITY_H
#define OHMMS_PARTICLELAYOUT_INPUTOUTPUT_UTILITY_H

#include "OhmmsData/OhmmsElementBase.h"
#include "Configuration.h"

namespace qmcplusplus
{

class LatticeParser
{

  typedef PtclOnLatticeTraits::ParticleLayout_t ParticleLayout_t;
  ParticleLayout_t& ref_;

public:

  LatticeParser(ParticleLayout_t& lat): ref_(lat) { }
  bool put(xmlNodePtr cur);
};


class LatticeXMLWriter
{

  typedef PtclOnLatticeTraits::ParticleLayout_t ParticleLayout_t;
  ParticleLayout_t& ref_;
public:

  LatticeXMLWriter(ParticleLayout_t& lat): ref_(lat) { }
  bool get(ostream& ) const;
  xmlNodePtr createNode();
};


}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
