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
    
    



#ifndef OHMMS_SPECIESBASE_H
#define OHMMS_SPECIESBASE_H

#include <string>
#include <vector>

/*! \class SpeciesBase
 *  \brief A class containing species attributes.
*/
class SpeciesBase
{
public:

  typedef double                     Scalar_t;
  typedef std::vector<Scalar_t>           SpeciesAttrib_t;
  typedef std::vector<SpeciesAttrib_t*>   AttribList_t;

  //!< The number of species
  unsigned        TotalNum;

  //!< Species Name
  std::vector<std::string>  Name;

  //!< List of species attributes
  AttribList_t    d_attrib;

  //!< Constructor
  SpeciesBase();

  //!< Destructor
  ~SpeciesBase();

  unsigned getTotalNum() const
  {
    return TotalNum;
  }
  void setTotalNum(const unsigned n)
  {
    TotalNum = n;
  }

  //!< returns the number of attributes in our list
  int numAttributes() const
  {
    return d_attrib.size();
  }

  //!< Adding a new attribute to the list
  int addAttrib();

  //!< Returns the value of j-th attribute for the i-th species
  inline double operator()(int i, int j) const
  {
    return d_attrib[i]->operator[](j);
  }

  //!< Assigns the value of i-th attribute for the j-th species
  inline double& operator()(int i, int j)
  {
    return d_attrib[i]->operator[](j);
  }

  //!< Add m species to the list by adding m elements to each attribute.
  void create(unsigned m)
  {
    if(m > 0)
    {
      Name.insert(Name.end(), m, std::string("none"));
      AttribList_t::iterator dit = d_attrib.begin();
      for(; dit != d_attrib.end(); ++dit)
        (*dit)->insert((*dit)->end(), m, 0);
      TotalNum += m;
    }
  }

  //!< Returns an ID for the species with name.
  unsigned find(const std::string& name) const
  {
    unsigned i=0;
    for(; i< TotalNum; i++)
    {
      if(Name[i] == name)
        return i;
    }
    return i;//Not found. Returns TotalNum.
  }

  const std::string& getName(unsigned int i) const
  {
    return Name[i];
  }

  //!< Returns an ID for the species with name.
  unsigned getSpeciesID(const std::string& name)
  {
    unsigned i = find(name); // check if the name is registered
    if(i == TotalNum)
      // if not found, add a new species
    {
      create(1);
      Name[i] = name;
    }
    return i; // return an index for a species
  }
};
#endif

