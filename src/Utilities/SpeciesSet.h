//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef OHMMS_SPECIES_SET_H
#define OHMMS_SPECIES_SET_H

#include <string>
#include <vector>

/*! \class SpeciesSet
 *  \brief A class containing a set of attributes for a set of species.
*/
class SpeciesSet
{

public:

  typedef double                     Scalar_t;
  typedef std::vector<Scalar_t>           SpeciesAttrib_t;
  typedef std::vector<SpeciesAttrib_t*>   AttribList_t;

  //! The number of species
  unsigned        TotalNum;

  //! Species name list
  std::vector<std::string> speciesName;

  //! attribute name list
  std::vector<std::string> attribName;

  //! List of species attributes
  AttribList_t    d_attrib;

  //! Constructor
  SpeciesSet();

  SpeciesSet(const SpeciesSet& species);

  SpeciesSet& operator=(const SpeciesSet& species);

  //! Destructor
  virtual ~SpeciesSet();

  ///return the number of species
  inline int size() const
  {
    return TotalNum;
  }
  ///return the number of species
  inline int getTotalNum() const
  {
    return TotalNum;
  }
  ///set the number of species
  inline void setTotalNum(const unsigned n)
  {
    TotalNum = n;
  }
  //! return the number of attributes in our list
  inline int numAttributes() const
  {
    return d_attrib.size();
  }

  /**
   * @param aname Unique name of the species be added.
   * @return the index of the species
   * @brief When a name species does not exist, add a new species
   */
  int addSpecies(const std::string& aname);

  /**
   * @param aname a unique name of an attribute
   * @return the index of a new attribute
   * @brief for a new attribute, allocate the data
   */
  int addAttribute(const std::string& aname);

  /**
   * @param aname Unique name of the species to be looked up.
   * @return the index of the species
   * @brief When a name species does not exist, return attribute.size()
   */
  int getAttribute(const std::string& aname);

  /**
   * @param i attribute index
   * @param j species index
   * @return the value of i-th attribute for the j-th species
   */
  inline double operator()(int i, int j) const
  {
    return d_attrib[i]->operator[](j);
  }

  /**
   * assignment operator
   * @param i attribute index
   * @param j species index
   * @return the value of i-th attribute for the j-th species
   */
  inline double& operator()(int i, int j)
  {
    return d_attrib[i]->operator[](j);
  }

  /**
   * @param m the number of species to be added
   */
  void create(unsigned m);

  /**
   * @param name a name of species
   * @return an ID for the species with name.
   * @brief if the input species is not found, add a new species
   */
  inline int findSpecies(const std::string& name) const
  {
    int i = 0;
    while(i< speciesName.size())
    {
      if(speciesName[i] == name)
        return i;
      i++;
    }
    return i;
  }

  inline int findAttribute(const std::string& name) const
  {
    return findIndex(name,attribName);
  }

  inline int findIndex(const std::string& name, const std::vector<std::string>& alist) const
  {
    int i = 0;
    while(i< alist.size())
    {
      if(alist[i] == name) return i;
      i++;
    }
    return -1;
  }
};
#endif

