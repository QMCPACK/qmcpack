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
#ifndef OHMMS_PARTICLE_INPUTOUTPUT_HDF_UTILITY_H
#define OHMMS_PARTICLE_INPUTOUTPUT_HDF_UTILITY_H

#include "OhmmsData/OhmmsElementBase.h"
#include "OhmmsData/RecordProperty.h"
#include "Particle/ParticleSet.h"
#include "OhmmsData/HDFAttribIO.h"

namespace qmcplusplus
{

class HDFParticleParser: public ParticleTags
{

public:

  typedef ParticleSet Particle_t;

  HDFParticleParser(Particle_t& aptcl):ref_(aptcl) { }

  ///reading from a file
  bool put(const char*);

  bool put(xmlNodePtr cur);

private:
  Particle_t& ref_;
};

class HDFSaveParticle:
  public ParticleTags,
  public RecordProperty
{

public:

  typedef ParticleSet Particle_t;

  HDFSaveParticle(Particle_t& pin): ref_(pin) { }

  ~HDFSaveParticle();

  void reset(const char* fileroot, bool append=false);

  void report(int iter);

  void finalize() { }

  bool put(xmlNodePtr cur);

private:

  Particle_t& ref_;
  string FileRoot;

};
}

#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
