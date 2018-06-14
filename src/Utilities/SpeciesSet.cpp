//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    




#include "Utilities/SpeciesSet.h"

SpeciesSet::SpeciesSet()
{
  TotalNum = 0;
  speciesName.reserve(4);
  attribName.reserve(4);
}

SpeciesSet::~SpeciesSet()
{
  AttribList_t::iterator it(d_attrib.begin());
  AttribList_t::iterator it_end(d_attrib.end());
  while(it != it_end)
  {
    delete *it;
    ++it;
  }
}


void SpeciesSet::create(unsigned m)
{
  if(m > 0)
  {
    speciesName.insert(speciesName.end(), m, std::string("none"));
    AttribList_t::iterator dit = d_attrib.begin();
    for(; dit != d_attrib.end(); ++dit)
      (*dit)->insert((*dit)->end(), m, 0);
    TotalNum += m;
  }
}

int SpeciesSet::addSpecies(const std::string& aname)
{
  int i = findSpecies(aname); // check if the name is registered
  if(i == TotalNum)
    // if not found, add a new species
  {
    create(1);
    speciesName[i] = aname;
  }
  return i; // return an index for a species
}

int SpeciesSet::addAttribute(const std::string& aname)
{
  int i = 0;
  while(i< attribName.size())
  {
    if(attribName[i] == aname)
      return i;
    i++;
  }
  attribName.push_back(aname);
  int n = d_attrib.size();
  d_attrib.push_back(new SpeciesAttrib_t(TotalNum));
  return n;
}

int SpeciesSet::getAttribute(const std::string& aname)
{
  for (int i=0; i<attribName.size(); i++)
  {
	if (attribName[i]== aname)
	  return i;  
  }
  return attribName.size();
}

SpeciesSet::SpeciesSet(const SpeciesSet& species):
  TotalNum(species.TotalNum),
  speciesName(species.speciesName),
  attribName(species.attribName)
{
  AttribList_t::const_iterator dit(species.d_attrib.begin());
  AttribList_t::const_iterator dit_end(species.d_attrib.end());
  while(dit != dit_end)
  {
    d_attrib.push_back(new SpeciesAttrib_t(**dit));
    ++dit;
  }
}

SpeciesSet& SpeciesSet::operator=(const SpeciesSet& species)
{
  if(this != &species)
  {
    TotalNum = species.TotalNum;
    speciesName = species.speciesName;
    attribName = species.attribName;
    AttribList_t::iterator it(d_attrib.begin());
    AttribList_t::iterator it_end(d_attrib.end());
    while(it != it_end)
    {
      delete *it;
      ++it;
    }
    d_attrib.clear();
    AttribList_t::const_iterator dit(species.d_attrib.begin());
    AttribList_t::const_iterator dit_end(species.d_attrib.end());
    while(dit != dit_end)
    {
      d_attrib.push_back(new SpeciesAttrib_t(**dit));
      ++dit;
    }
  }
  return *this;
}


