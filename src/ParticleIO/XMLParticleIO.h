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
#ifndef OHMMS_PARTICLE_INPUTOUTPUT_XML_UTILITY_H
#define OHMMS_PARTICLE_INPUTOUTPUT_XML_UTILITY_H

#include "OhmmsData/OhmmsElementBase.h"
#include "OhmmsData/RecordProperty.h"
#include "Particle/ParticleSet.h"

namespace qmcplusplus
{

//   struct XMLParticleBase {
//     static string root_tag;
//     static string attrib_tag;
//     static string datatype_tag;
//     static string condition_tag;
//     static string size_tag;
//   };

class XMLParticleParser: public ParticleTags
{

  typedef ParticleSet                  Particle_t;
  typedef Particle_t::ParticleIndex_t  ParticleIndex_t;
  typedef Particle_t::ParticleScalar_t ParticleScalar_t;
  typedef Particle_t::ParticlePos_t    ParticlePos_t;
  typedef Particle_t::ParticleTensor_t ParticleTensor_t;

  bool AssignmentOnly;
  Particle_t& ref_;
  Tensor<int,OHMMS_DIM>& TileMatrix;

  bool putSpecial(xmlNodePtr cur);

  /** read the data of a particle attribute
   *@param cur the xmlnode
   *@param nat the number of particle attributes to be read
   *@param nloc the current local count to which nat particle attributes are added.
   */
  void getPtclAttrib(xmlNodePtr cur, int nat, int nloc);

public:

  /**constructor
   *@param aptcl the particleset to be initialized
   *@param donotresize if true, only assignment is done
   */
  XMLParticleParser(Particle_t& aptcl, Tensor<int,OHMMS_DIM>& tmat
                    , bool donotresize=false);

  ///reading from a file
  bool put(const string& fname_in, const string& fext_in);

  bool put(xmlNodePtr cur);

  /** reset the properties of a particle set
   */
  bool reset(xmlNodePtr cur);

};

class XMLSaveParticle:
  public ParticleTags,
  public RecordProperty
{

  typedef ParticleSet                  Particle_t;
  typedef Particle_t::ParticleIndex_t  ParticleIndex_t;
  typedef Particle_t::ParticleScalar_t ParticleScalar_t;
  typedef Particle_t::ParticlePos_t    ParticlePos_t;
  typedef Particle_t::ParticleTensor_t ParticleTensor_t;

  Particle_t& ref_;
  string FileRoot;
  vector<string> SpeciesName;

public:

  XMLSaveParticle(Particle_t& pin);

  ~XMLSaveParticle();

  void reset(const char* fileroot, bool append=false);

  void report(int iter);

  void finalize() { }

  bool put(xmlNodePtr cur);

  void get(ostream& os, int olevel) const;

  xmlNodePtr createNode(bool addlattice);

private:

};
}

#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
