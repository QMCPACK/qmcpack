//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef OHMMS_PARTICLE_INPUTOUTPUT_XML_UTILITY_H
#define OHMMS_PARTICLE_INPUTOUTPUT_XML_UTILITY_H

#include "OhmmsData/OhmmsElementBase.h"
#include "OhmmsData/RecordProperty.h"
#include "Particle/ParticleSet.h"

namespace qmcplusplus
{

//   struct XMLParticleBase {
//     static std::string root_tag;
//     static std::string attrib_tag;
//     static std::string datatype_tag;
//     static std::string condition_tag;
//     static std::string size_tag;
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
  bool put(const std::string& fname_in, const std::string& fext_in);

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
  std::string FileRoot;
  std::vector<std::string> SpeciesName;

public:

  XMLSaveParticle(Particle_t& pin);

  ~XMLSaveParticle();

  void reset(const char* fileroot, bool append=false);

  void report(int iter);

  void finalize() { }

  bool put(xmlNodePtr cur);

  void get(std::ostream& os, int olevel) const;

  xmlNodePtr createNode(bool addlattice);

private:

};
}

#endif

