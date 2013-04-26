//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#ifndef OHMMS_SPECIESCOLLECTION_H
#define OHMMS_SPECIESCOLLECTION_H

#include <iostream>
#include <vector>
using namespace std;

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
  typedef vector<string>               SpeciesAttribMap_t;

  static SpeciesBase* getSpecies();//!< Returns singleton SpeciesBase*
  static int addAttrib(const char* name);//!< Adds an attribute
  static unsigned int addSpecies(const char* name)
  {
    return mySpecies->getSpeciesID(name);
  }//!< Adds a species

  static void print(ostream& );

protected:

  SpeciesCollection() {}
  ~SpeciesCollection();

private:
  static SpeciesAttribMap_t attribMap;//!< Mapping between attribute name and its index
  static SpeciesBase* mySpecies;//!< Singleton SpeciesBase*
};


#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
