//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef OHMMS_PARTICLEBASE_H
#define OHMMS_PARTICLEBASE_H

#include <string>
#include <vector>
#include <map>
#include "ParticleTags.h"

namespace qmcplusplus
{
/**@file ParticleBase.h
 *@brief Declaration of ParticleBase<PT>
 */

/** Container class for particle attributes.
 *
 * The template parameter PT is a trait class which provides typedefs.
 * See, ParticleTraits.h
 */
template<class PT>
class ParticleBase: public ParticleTags
{

public:

  typedef typename PT::ParticleLayout_t      ParticleLayout_t;
  typedef typename PT::Index_t               Index_t;
  typedef typename PT::Scalar_t              Scalar_t;
  typedef typename PT::SingleParticleIndex_t SingleParticleIndex_t;
  typedef typename PT::SingleParticlePos_t   SingleParticlePos_t;
  typedef typename PT::Tensor_t              Tensor_t;
  typedef typename PT::ParticleValue_t       ParticleValue_t;

  //@{containers
  //! Three types of particle attributes are available.
  //! Integer quantity associated with each atom, e.g., ID, GroupID
  typedef typename PT::ParticleIndex_t       ParticleIndex_t;

  //! Scalar quantity associated with each atom, e.g., Local Energy
  typedef typename PT::ParticleScalar_t      ParticleScalar_t;

  //! D-dimensional vector quantity associated with each atom, e.g., R
  typedef typename PT::ParticlePos_t         ParticlePos_t;

  //! D-dimensional tensor quantity associated with each atom, e.g., R
  typedef typename PT::ParticleTensor_t      ParticleTensor_t;

  //! D-dimensional real/complex quantity associated with each atom, e.g., R
  typedef typename PT::ParticleGradient_t     ParticleGradient_t;

  //! D-dimensional real/complex quantity associated with each atom, e.g., R
  typedef typename PT::ParticleLaplacian_t     ParticleLaplacian_t;
  //@}

  //! Mapping between an attribute name and particle attribute
  typedef std::map<std::string,OhmmsObject*>::iterator  PAListIterator;


  /*! \note Lattice and three particle attributes, ID, GroupID, and R are
   * public. Encapsulation of ParticleBase is realized by making
   * a ParticleBase a protected member of PropagatorBase.
   * Modification of particle attributes is only allowed by PropagatorBase
   * and its derived classes.
   */
  //!< ParticleLayout
  ParticleLayout_t  Lattice, PrimitiveLattice;
  //!< unique, persistent ID for each particle
  ParticleIndex_t   ID;
  //!< Species ID
  ParticleIndex_t   GroupID;
  //!< Position
  ParticlePos_t    R;
  //!< Instantaneous Position in unit used by layout/nnlist engines, not added to myAttribList
  //ParticlePos_t    curR;
  //!< Default constructor
  ParticleBase();

  ParticleBase(const ParticleBase& P): Counter(0), LocalNum(0), GlobalNum(0)
  {
    initBase();
    assign(P);
  }

  virtual ~ParticleBase();

  ///return a type id: one of the enum values
  inline int getAttribType(const std::string& tname)
  {
    return AttribTypeMap[tname];
  }
  inline int getNumAttrib() const
  {
    return AttribList.size();
  }

  inline PAListIterator first_attrib()
  {
    return AttribList.begin();
  }
  inline PAListIterator last_attrib()
  {
    return AttribList.end();
  }

  bool hasAttrib(const std::string& attrib_name);

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
      if(pa == nullptr) //only 
      {
        pa=new attrib_type(tname,oname);
        pa->resize(LocalNum);
        pa->setID(AttribList.size());

        AllocatedList.push_back(pa);
        AttribList[oname]=pa;
      }
    }
    return pa;
  }

  void createBase(size_t m);
  void createBase(const std::vector<int>& agroup);

  void resize(size_t m);
  void clear();

  virtual void assign(const ParticleBase<PT>& ptclin)
  {
    size_t nextra=ptclin.AllocatedList.size();
    if(nextra)
    {
      AllocatedList.resize(nextra);
      for(size_t i=0; i<nextra; ++i)
      {
        OhmmsObject* o=ptclin.AllocatedList[i]->makeClone();
        AllocatedList[i]=o;
        AttribList[o->objName()]=o;
      }
    }

    resize(ptclin.getLocalNum());
    Lattice = ptclin.Lattice;
    PrimitiveLattice = ptclin.PrimitiveLattice;
    R.InUnit = ptclin.R.InUnit;
    R = ptclin.R;
    ID = ptclin.ID;
    GroupID = ptclin.GroupID;
    if(ptclin.SubPtcl.size())
    {
      SubPtcl.resize(ptclin.SubPtcl.size());
      SubPtcl =ptclin.SubPtcl;
    }
  }

  inline int getLocalNum() const
  {
    return LocalNum;
  }
  inline int getTotalNum() const
  {
    return GlobalNum;
  }

  ///return the number of groups
  inline int groups() const
  {
    return SubPtcl.size()-1;
  }

  ///return the first index of a group i
  inline int first(int igroup) const
  {
    return SubPtcl[igroup];
  }

  ///return the last index of a group i
  inline int last(int igroup) const
  {
    return SubPtcl[igroup+1];
  }

  /// Returns the current counter
  inline int  current() const
  {
    return Counter;
  }

  /// Increments the counter
  inline void advance()
  {
    Counter++;
  }

  /// Sets the counter
  inline void setCounter(int i = 0)
  {
    Counter = i;
  }

protected:

  void initBase();

  //!< Internal counter
  int Counter;
  int LocalNum;
  int GlobalNum;

  ///array to handle a group of distinct particles per species
  ParticleIndex_t                       SubPtcl;
  std::map<std::string,int>             AttribTypeMap;
  std::map<std::string,OhmmsObject*>    AttribList;
  /** objects created by getXYZAttrib(aname) */
  std::vector<OhmmsObject*>        AllocatedList;

  /** add ParticleAttrib<AT>
   * @tparm AT any element type, int, double, float ...
   */
  template<typename AT> 
    int addAttribute(ParticleAttrib<AT>& pa)
  {
    if(pa.size() < getLocalNum()) pa.resize(getLocalNum());

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

};

/** Default constructor */
template<class PL>
ParticleBase<PL>::ParticleBase():Counter(0), LocalNum(0), GlobalNum(0)
{
  initBase();
}

template<class PL>
void ParticleBase<PL>::initBase()
{
  R.setTypeName(ParticleTags::postype_tag);
  R.setObjName(ParticleTags::position_tag);
  ID.setTypeName(ParticleTags::indextype_tag);
  ID.setObjName(ParticleTags::id_tag);
  GroupID.setTypeName(ParticleTags::indextype_tag);
  GroupID.setObjName(ParticleTags::ionid_tag);
  //map[type-name] to enum
  AttribTypeMap[ParticleTags::indextype_tag]  = PA_IndexType;
  AttribTypeMap[ParticleTags::scalartype_tag] = PA_ScalarType;
  AttribTypeMap[ParticleTags::stringtype_tag] = PA_StringType;
  AttribTypeMap[ParticleTags::postype_tag]    = PA_PositionType;
  AttribTypeMap[ParticleTags::tensortype_tag] = PA_TensorType;
  //add basic attributes
  addAttribute(R);
  addAttribute(ID);
  addAttribute(GroupID);
  //curR is always in unit
  //curR.setUnit(PosUnit::LatticeUnit);
}

template<class PL>
ParticleBase<PL>::~ParticleBase()
{
  for(int i=0; i<AllocatedList.size(); i++)
    delete AllocatedList[i];
}

/** check if named attribute exists
 *  \param aname String for the attribute name
 *  \return true if the name is found.
 */
template<class PL>
bool ParticleBase<PL>::hasAttrib(const std::string& aname)
{
  return (AttribList.find(aname) != AttribList.end());
}

template<class PL>
void ParticleBase<PL>::createBase(size_t m)
{
  std::map<std::string,OhmmsObject*>::iterator it= AttribList.begin();
  while(it!=AttribList.end())
  {
    (*it).second->create(m);
    ++it;
  }

  //curR.create(m);
  LocalNum += m;
  GlobalNum += m;
}

template<class PL>
void ParticleBase<PL>::resize(size_t m)
{
  std::map<std::string,OhmmsObject*>::iterator it= AttribList.begin();
  while(it!=AttribList.end())
  {
    (*it).second->resize(m);
    ++it;
  }

  //curR.resize(m);
  LocalNum = m;
  GlobalNum = m;
}

template<class PL>
void ParticleBase<PL>::clear()
{
  std::map<std::string,OhmmsObject*>::iterator it= AttribList.begin();
  while(it!=AttribList.end())
  {
    (*it).second->clear();
    ++it;
  }
  //curR.clear();
  LocalNum = 0;
  GlobalNum = 0;
}

/**function to create N-particle system
 *@param agroup an integer array containing the number of particles belonging
 *to a subgroup.
 *@brief allocate the ParticleAttributes, such as R, G, L.
 *The size of n is the number of distinct subgroups. SubPtcl class
 *is used to efficient evaluate the subgroup properties.
 */
template<class PL>
void ParticleBase<PL>::createBase(const std::vector<int>& agroup)
{
  SubPtcl.resize(agroup.size()+1);
  SubPtcl[0] = 0;
  for(int is=0; is<agroup.size(); is++)
    SubPtcl[is+1] = SubPtcl[is]+agroup[is];
  size_t nsum = SubPtcl[agroup.size()];
  resize(nsum);
  int loc=0;
  for(int i=0; i<agroup.size(); i++)
  {
    for(int j=0; j<agroup[i]; j++,loc++)
      GroupID[loc] = i;
  }
}

template<class T, unsigned D>
inline T* get_first_address(ParticleAttrib<TinyVector<T,D> >& a)
{
  return &(a[0][0]);
}

template<class T, unsigned D>
inline T* get_last_address(ParticleAttrib<TinyVector<T,D> >& a)
{
  return &(a[0][0])+D*a.size();
}

}
#endif

