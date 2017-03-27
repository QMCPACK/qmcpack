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

void SpeciesCollection::print(std::ostream& os)
{
  os << "Total Number of Species Attributes = " << mySpecies->numAttributes() << std::endl;
  os << "Total Number of Species  = " << mySpecies->getTotalNum() << std::endl;
  for(int isp=0; isp<mySpecies->getTotalNum(); isp++)
  {
    os  << mySpecies->Name[isp] << " " ;
    for(int id=0; id< mySpecies->numAttributes(); id++)
    {
      os << mySpecies->operator()(id,isp) << " ";
    }
    os << std::endl;
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

int SpeciesCollection::addSpecies(std::vector<std::string>& argv) {

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

