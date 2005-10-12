//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002,2003- by Jeongnim Kim
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

#include "Utilities/SpeciesSet.h"
using namespace std;

SpeciesSet::SpeciesSet() {
  TotalNum = 0;  
  Name.reserve(10); // expect less than 10 species
}

SpeciesSet::SpeciesSet(const SpeciesSet& species): 
  TotalNum(species.TotalNum),
  Name(species.Name),
  attrib_ind(species.attrib_ind){
  AttribList_t::const_iterator dit(species.d_attrib.begin());
  AttribList_t::const_iterator dit_end(species.d_attrib.end());
  while(dit != dit_end) {
    d_attrib.push_back(new SpeciesAttrib_t(**dit));++dit;
  }
}

SpeciesSet& SpeciesSet::operator=(const SpeciesSet& species) {

  TotalNum = species.TotalNum;
  Name = species.Name;
  attrib_ind = species.attrib_ind;

  AttribList_t::iterator it(d_attrib.begin());
  AttribList_t::iterator it_end(d_attrib.end());
  while(it != it_end) { delete *it; ++it;}
  d_attrib.clear();

  AttribList_t::const_iterator dit(species.d_attrib.begin());
  AttribList_t::const_iterator dit_end(species.d_attrib.end());
  while(dit != dit_end) {
    d_attrib.push_back(new SpeciesAttrib_t(**dit));++dit;
  }
  return *this;
}

SpeciesSet::~SpeciesSet() {
  AttribList_t::iterator dit = d_attrib.begin();      
  for(; dit != d_attrib.end(); ++dit) delete (*dit);
}

// size_type SpeciesSet::addAttrib() {
//   size_type n = d_attrib.size();
//   d_attrib.push_back(new SpeciesAttrib_t(TotalNum));
//   return n;
// }

/** add a new attribute to the list
 *@param aname the name of an attribute
 */
SpeciesSet::size_type
SpeciesSet::addAttrib(const char* aname) {
  AttribMap_t::iterator it = attrib_ind.find(aname);
  if(it == attrib_ind.end()) {
    //size_type a=d_attrib.size();
    size_type a=d_attrib.size();
    attrib_ind[aname] = a;
    d_attrib.push_back(new SpeciesAttrib_t(TotalNum));
    return a;
  } else {
    return (*it).second;
  }
}

bool SpeciesSet::get(ostream& os) const {

  os << "Total Number of Species Attributes = " << d_attrib.size() << endl;
  os << "Total Number of Species  = " << TotalNum << endl;

  for(size_type isp=0; isp<TotalNum; isp++)  {
    os  << Name[isp] << " " ;
    for(size_type id=0; id< d_attrib.size(); id++) {
      os << d_attrib[id]->operator[](isp) << " ";
    }
    os << endl;
  }
  return true;
}

bool SpeciesSet::put(istream& is) {
  return true;
}

void SpeciesSet::reset() {

}

bool SpeciesSet::put(xmlNodePtr cur) {
  Scalar_t val;
  cur = cur->children;
  while(cur != NULL) {
    string cname((const char*)(cur->name));
    if(cname == "group") {
      const xmlChar* s=xmlGetProp(cur,(const xmlChar*)"name");
      if(s) {
	size_type species_id= getSpeciesID((const char*)s);
        xmlNodePtr tcur = cur->children;
        while(tcur != NULL) {	//check parameters
          string tcname((const char*)tcur->name);
          if(tcname == "parameter") {
            const xmlChar* t=xmlGetProp(tcur,(const xmlChar*)"name");
            if(t) {
              size_type attrib_id = addAttrib((const char*)t);
              putContent(val,tcur);
              d_attrib[attrib_id]->operator[](species_id)=val;
            }
          }
          tcur = tcur->next;
        }
      }
    } 
    cur = cur->next;
  }
  return true;
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

