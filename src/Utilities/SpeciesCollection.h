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
    
    


#ifndef OHMMS_SPECIESCOLLECTION_H
#define OHMMS_SPECIESCOLLECTION_H

#include <iostream>
#include <vector>

#ifndef OHMMS_SPECIESBASE_H
#include "Utilities/OhmmsSpecies.h"
#endif

/*! \class SpeciesCollection
 *  \brief Responsible for instatiating singleton SpeciesBase
(Singleton Pattern). For a simulation, only one set of species should exist.
 *
 *  Any object can request a pointer to the SpeciesBase by calling SpeciesCollection::getSpecies().
*/

class SpeciesCollection
{

public:

  typedef SpeciesBase::SpeciesAttrib_t SpeciesAttrib_t;
  typedef std::vector<std::string>               SpeciesAttribMap_t;

  static SpeciesBase* getSpecies();//!< Returns singleton SpeciesBase*
  static int addAttrib(const char* name);//!< Adds an attribute
  static unsigned int addSpecies(const char* name)
  {
    return mySpecies->getSpeciesID(name);
  }//!< Adds a species

  static void print(std::ostream& );

protected:

  SpeciesCollection() {}
  ~SpeciesCollection();

private:
  static SpeciesAttribMap_t attribMap;//!< Mapping between attribute name and its index
  static SpeciesBase* mySpecies;//!< Singleton SpeciesBase*
};


#endif

