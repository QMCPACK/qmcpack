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

#include "Utilities/SpeciesSet.h"
using namespace std;

SpeciesSet::SpeciesSet() {
  TotalNum = 0;  
  speciesName.reserve(10); 
  attribName.reserve(10);
}

SpeciesSet::~SpeciesSet() {
  AttribList_t::iterator dit = d_attrib.begin();      
  for(; dit != d_attrib.end(); ++dit) delete (*dit);
}


void SpeciesSet::create(unsigned m) {
  if(m > 0) {
    speciesName.insert(speciesName.end(), m, string("none"));
    AttribList_t::iterator dit = d_attrib.begin();      
    for(; dit != d_attrib.end(); ++dit) (*dit)->insert((*dit)->end(), m, 0);
    TotalNum += m;
  }
}

int SpeciesSet::addSpecies(const string& aname) {    

  int i = findSpecies(aname); // check if the name is registered
  if(i == TotalNum) { // if not found, add a new species
    create(1);
    speciesName[i] = aname;
  }
  return i; // return an index for a species
}

int SpeciesSet::addAttribute(const string& aname) {

  int i = 0;
  while(i< attribName.size()) {
    if(attribName[i] == aname) return i; 
    i++;
  }
  attribName.push_back(aname);
  int n = d_attrib.size();
  d_attrib.push_back(new SpeciesAttrib_t(TotalNum));
  return n;
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

