//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//
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
/**@file SpeciesSet
 * @brief Declaraton of SpeciesSet
 */
#ifndef OHMMS_SPECIESBASE_H
#define OHMMS_SPECIESBASE_H
#ifdef HAVE_CONFIG_H
#include "ohmms-config.h"
#endif
#include <string>
#include <vector>
#include <map>
#include "OhmmsData/OhmmsElementBase.h"

/** A class containing species attributes.
*/
class SpeciesSet: public OhmmsElementBase {

public:

  typedef OHMMS_PRECISION                  Scalar_t;
  typedef std::size_t                      size_type;
  typedef std::vector<Scalar_t>            SpeciesAttrib_t;
  typedef std::vector<SpeciesAttrib_t* >   AttribList_t;
  typedef std::map<std::string,size_type>  AttribMap_t;
  ///the number of species
  size_type        TotalNum;

  /// Species Name list
  std::vector<std::string>  Name;   

  /// Constructor
  SpeciesSet();

  /// Copy constructor, deep copy
  SpeciesSet(const SpeciesSet& species);

  /// Destructor
  ~SpeciesSet();

  SpeciesSet& operator=(const SpeciesSet& species);

  ///return the number of species
  size_type getTotalNum() const { return TotalNum; }

  ///set the number of species
  void setTotalNum(size_type n) { TotalNum = n; }

  ///return the number of attributes in our list
  size_type numAttributes() const { return d_attrib.size();}

  /// add a new attribute to the list
  size_type addAttrib(const char* aname);

  /// return the value of i-th attribute for the j-th species
  inline double operator()(size_type i, size_type j) const { 
    return d_attrib[i]->operator[](j); 
  }

  ///assign the value of i-th attribute for the j-th species
  inline double& operator()(size_type i, size_type j) { 
    return d_attrib[i]->operator[](j); 
  }

  inline const SpeciesAttrib_t& getAttribute(size_type i) const {
    return *d_attrib[i];
  }

  /// add m species to the list by adding m elements to each attribute.
  inline void create(size_type m) {
    if(m > 0) {
      Name.insert(Name.end(), m, std::string("none"));
      AttribList_t::iterator dit(d_attrib.begin());      
      AttribList_t::iterator dit_end(d_attrib.end());      
      while(dit != dit_end) { 
        (*dit)->insert((*dit)->end(), m, 0);
        ++dit;
      }
      TotalNum += m;
    }
  }

  /// return an ID for the species with name. 
  size_type find(const std::string& name) const{
    size_type i=0;
    for(; i< TotalNum; i++) {
      if(Name[i] == name) return i; 
    }
    return i;//Not found. Returns TotalNum.
  }

  /// return the name of the ith species
  const std::string& getName(size_type i) const { return Name[i]; }

  /// return an ID for the species with name. Add a species if not found
  inline size_type getSpeciesID(const std::string& name) { 
    size_type i(find(name)); // check if the name is registered
    if(i == TotalNum) { // if not found, add a new species
      create(1);
      Name[i] = name;
    }
    return i; // return an index for a species
  }

  inline size_type addSpecies(const std::string&  aname) {    
    return getSpeciesID(aname);
  }

  bool get(std::ostream& os) const;
  bool put(std::istream& is);
  void reset();
  bool put(xmlNodePtr cur);

private:
  /// List of species attributes
  AttribList_t d_attrib;

  ///map for the name of an attribute to the attribute index
  AttribMap_t attrib_ind;

};
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
