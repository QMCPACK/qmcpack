//////////////////////////////////////////////////////////////////
// (c) Copyright 2009- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file ESHDFParticleParser.h
 * @brief Definition of the parsers of ions and electrons from a ESHDF file.
 */
#ifndef ESHDF_PARTICLE_PARSER_H
#define ESHDF_PARTICLE_PARSER_H

#include "Particle/ParticleSet.h"
#include "OhmmsData/HDFAttribIO.h"

class Communicate;

namespace qmcplusplus
{

struct ESHDFElectronsParser
{

  ESHDFElectronsParser(ParticleSet& aptcl, hid_t h=-1, Communicate* c=0);

  bool put(xmlNodePtr cur);

  ParticleSet& ref_;
  hid_t hfile_id;
  Communicate* myComm;
};

struct ESHDFIonsParser
{

  ESHDFIonsParser(ParticleSet& aptcl, hid_t h=-1, Communicate* c=0);

  bool put(xmlNodePtr cur);
  void readESHDF();
  //expand the ionic systems
  void expand(Tensor<int,3>& tilematrix);

  ParticleSet& ref_;
  hid_t hfile_id;
  Communicate* myComm;
  string atomic_number_tag;
  string charge_tag;
  string mass_tag;
};

}

#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 757 $   $Date: 2005-10-31 10:10:28 -0600 (Mon, 31 Oct 2005) $
 * $Id: ESHDFParticleParser.h 757 2005-10-31 16:10:28Z jnkim $
 ***************************************************************************/
