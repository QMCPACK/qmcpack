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

class AttribListType: public ParticleTags
{
  std::map<std::string,int>             AttribTypeMap;
  std::map<std::string,OhmmsObject*>    AttribList;
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
   * @tparm AT any element type, int, double, float ...
   */
  template<typename AT>
    int add(ParticleAttrib<AT>& pa)
  {
    int oid=AttribList.size();
    std::map<std::string,OhmmsObject*>::iterator it= AttribList.find(pa.objName());
    if(it == AttribList.end())
    {
      AttribList[pa.objName()]=&pa;
      pa.setID(oid);
    }
    else
    {
      oid=(*it).second->id();
    }
    return oid;
  }

  ///return a type id: one of the enum values
  inline int getAttribType(const std::string& tname)
  {
    return AttribTypeMap[tname];
  }

  /** generic get function attribute function
   * @param tname attribute type name
   * @param oname attribute name
   * @param pa pointer to ParticleAttrib<AT>*
   * @return pointer to the attribute
   */
  template<typename AT>
    ParticleAttrib<AT>* getAttribute(const std::string& tname, const std::string& oname, ParticleAttrib<AT>*  pa)
  {
    bool foundit=false;
    typedef ParticleAttrib<AT> attrib_type;
    int oid=AttribList.size();
    std::map<std::string,OhmmsObject*>::iterator it= AttribList.find(oname);
    if(it != AttribList.end())
    {
      OhmmsObject* o=(*it).second;
      return dynamic_cast<attrib_type*>(o);
    }
    else
    {
      APP_ABORT("AttribListType::getAttribute Unknown attribute "+oname+"\n");
      /*
      if(pa == nullptr) //only
      {
        pa=new attrib_type(tname,oname);
        pa->resize(LocalNum);
        pa->setID(AttribList.size());

        AllocatedList.push_back(pa);
        AttribList[oname]=pa;
      }
      */
    }
    return pa;
  }

};

class XMLParticleParser: public ParticleTags
{

  typedef ParticleSet                  Particle_t;
  typedef Particle_t::ParticleIndex_t  ParticleIndex_t;
  typedef Particle_t::ParticleScalar_t ParticleScalar_t;
  typedef Particle_t::ParticlePos_t    ParticlePos_t;
  typedef Particle_t::ParticleTensor_t ParticleTensor_t;

  bool AssignmentOnly;
  Particle_t& ref_;
  AttribListType ref_AttribList;
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
  AttribListType ref_AttribList;
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

