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
class AttribListType : public ParticleTags
{
  std::map<std::string, int> AttribTypeMap;
  std::map<std::string, OhmmsObject*> AttribList;
  ///** objects created by getXYZAttrib(aname) */
  //std::vector<OhmmsObject*>             AllocatedList;

public:
  AttribListType()
  {
    //map[type-name] to enum
    AttribTypeMap[ParticleTags::indextype_tag]  = PA_IndexType;
    AttribTypeMap[ParticleTags::scalartype_tag] = PA_ScalarType;
    AttribTypeMap[ParticleTags::stringtype_tag] = PA_StringType;
    AttribTypeMap[ParticleTags::postype_tag]    = PA_PositionType;
    AttribTypeMap[ParticleTags::tensortype_tag] = PA_TensorType;
  }

  /*
  ~AttribListType()
  {
    for(int i=0; i<AllocatedList.size(); i++)
      delete AllocatedList[i];
  }
  */

  /** add ParticleAttrib<AT>
   * @tparam AT any element type, int, double, float ...
   */
  template<typename AT>
  int add(ParticleAttrib<AT>& pa)
  {
    int oid                                          = AttribList.size();
    std::map<std::string, OhmmsObject*>::iterator it = AttribList.find(pa.objName());
    if (it == AttribList.end())
    {
      AttribList[pa.objName()] = &pa;
      pa.setID(oid);
    }
    else
    {
      oid = (*it).second->id();
    }
    return oid;
  }

  ///return a type id: one of the enum values
  inline int getAttribType(const std::string& tname) { return AttribTypeMap[tname]; }

  /** generic get function attribute function
   * @param tname attribute type name
   * @param oname attribute name
   * @return pointer to the attribute
   */
  template<typename AT>
  ParticleAttrib<AT>* getAttribute(const std::string& tname, const std::string& oname)
  {
    using attrib_type = ParticleAttrib<AT>;
    const auto it     = AttribList.find(oname);
    if (it != AttribList.end())
    {
      OhmmsObject* o = (*it).second;
      return dynamic_cast<attrib_type*>(o);
    }
    else
      throw std::runtime_error("AttribListType::getAttribute Unknown attribute " + oname + "\n");
    return nullptr;
  }
};

class XMLParticleParser : public ParticleTags
{
  using Particle_t     = ParticleSet;
  using ParticleIndex  = Particle_t::ParticleIndex;
  using ParticleScalar = Particle_t::ParticleScalar;
  using ParticlePos    = Particle_t::ParticlePos;
  using ParticleTensor = Particle_t::ParticleTensor;

  Particle_t& ref_;
  AttribListType ref_AttribList;

  /** read the data of a particle attribute
   *@param cur the xmlnode
   *@param in_offset the location offset to read from XML element node body.
   *@param copy_size the number of particle attributes to be read
   *@param out_offset the current local count to which copy_size particle attributes are added.
   */
  void getPtclAttrib(xmlNodePtr cur, int in_offset, int copy_size, int out_offset);

  void checkGrouping(int nat, const std::vector<int>& nat_group) const;

public:
  /**constructor
   *@param aptcl the particleset to be initialized
   */
  XMLParticleParser(Particle_t& aptcl);

  bool readXML(xmlNodePtr cur);

  /** reset the properties of a particle set
   */
  bool reset(xmlNodePtr cur);
};

class XMLSaveParticle : public ParticleTags, public RecordProperty
{
  using Particle_t     = ParticleSet;
  using ParticleIndex  = Particle_t::ParticleIndex;
  using ParticleScalar = Particle_t::ParticleScalar;
  using ParticlePos    = Particle_t::ParticlePos;
  using ParticleTensor = Particle_t::ParticleTensor;

  Particle_t& ref_;
  AttribListType ref_AttribList;
  std::string FileRoot;
  std::vector<std::string> SpeciesName;

public:
  XMLSaveParticle(Particle_t& pin);

  ~XMLSaveParticle() override;

  void reset(const char* fileroot, bool append = false) override;

  void report(int iter) override;

  void finalize() override {}

  bool put(xmlNodePtr cur) override;

  void get(std::ostream& os, int olevel) const;

  xmlNodePtr createNode(bool addlattice);

private:
};
} // namespace qmcplusplus

#endif
