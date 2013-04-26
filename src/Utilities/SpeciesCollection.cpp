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

#include "Utilities/SpeciesCollection.h"
#include "Utilities/OhmmsInfo.h"

SpeciesBase* SpeciesCollection::mySpecies = NULL;
SpeciesCollection::SpeciesAttribMap_t SpeciesCollection::attribMap;

SpeciesCollection::~SpeciesCollection()
{
  DEBUGMSG("Calling SpeciesCollection::~SpeciesCollection()");
  if(mySpecies)
    delete mySpecies;
}
SpeciesBase* SpeciesCollection::getSpecies()
{
  if(!mySpecies)
  {
    mySpecies = new SpeciesBase;
  }
  return mySpecies;
}

/*! \fn SpeciesCollecton::addAttrib(const char* name)
 *  \param name Unique name of the species to be added.
 */
int SpeciesCollection::addAttrib(const char* name)
{
  for(int i=0; i<attribMap.size(); i++)
  {
    if(attribMap[i] == name)
      return i;
  }
  attribMap.push_back(name);
  mySpecies->addAttrib();
  return attribMap.size() -1;
}

void SpeciesCollection::print(ostream& os)
{
  os << "Total Number of Species Attributes = " << mySpecies->numAttributes() << endl;
  os << "Total Number of Species  = " << mySpecies->getTotalNum() << endl;
  for(int isp=0; isp<mySpecies->getTotalNum(); isp++)
  {
    os  << mySpecies->Name[isp] << " " ;
    for(int id=0; id< mySpecies->numAttributes(); id++)
    {
      os << mySpecies->operator()(id,isp) << " ";
    }
    os << endl;
  }
}
/*
int SpeicesCollection::addSpecies(int argc, char **argv) {
  int i=0;
  int id = -1;
  do { // first check the the species id
    if(!strcmp(argv[i], "species")) {
      id = mySpecies->getSpeciesID(argv[++i].c_str());
    }
  } while(id<0);
  return id;
}

int SpeciesCollection::addSpecies(vector<string>& argv) {

  int i=0;
  int id = -1;
  do { // first check the the species id
    if(argv[i] == "species") {
      id = mySpecies->getSpeciesID(argv[++i].c_str());
    }
  } while(id<0);
  return id;
}
*/

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
