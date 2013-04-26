//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002,2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-

#include "Utilities/SpeciesSet.h"
using namespace std;

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
    speciesName.insert(speciesName.end(), m, string("none"));
    AttribList_t::iterator dit = d_attrib.begin();
    for(; dit != d_attrib.end(); ++dit)
      (*dit)->insert((*dit)->end(), m, 0);
    TotalNum += m;
  }
}

int SpeciesSet::addSpecies(const string& aname)
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

int SpeciesSet::addAttribute(const string& aname)
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

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

