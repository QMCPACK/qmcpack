//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
#ifndef OHMMS_QMC_SIMPLE_UTILITIES_H
#define OHMMS_QMC_SIMPLE_UTILITIES_H

#include "QMCWaveFunctions/OrbitalBuilderBase.h"

namespace ohmmsqmc {

  inline  
  bool 
  determineNumOfElectrons(ParticleSet& el, xmlXPathContextPtr acontext){
    
    //initialize el with the wave function information
    //This is to determine the number of up and down electrons
    vector<int> N;
    string sdname("//");
    sdname.append(OrbitalBuilderBase::sd_tag);
    cout << "***** The xpath for to find SlaterDeterminant " << sdname << endl;
    xmlXPathObjectPtr result
      = xmlXPathEvalExpression((const xmlChar*)(sdname.c_str()),acontext);
    
    bool found_el=false;
    int nsd= result->nodesetval->nodeNr;
    XMLReport("Found " << nsd << " SlaterDeterminant")
    if(nsd) {
      vector<xmlNodePtr> dset;
      xmlNodePtr cur=result->nodesetval->nodeTab[0]->children;
      while(cur != NULL) {
	string cname((const char*)(cur->name));
	if(cname == OrbitalBuilderBase::det_tag) dset.push_back(cur);
	cur=cur->next;
      }
      if(dset.size()) {
	XMLReport("Found " << dset.size() << OrbitalBuilderBase::det_tag)
	for(int i=0; i<dset.size(); i++) {
	  int norb = 0;
	  xmlChar* orb=xmlGetProp(dset[i],(const xmlChar*)"orbitals");
	  if(orb){
	    norb = atoi((const char*)orb);
	    XMLReport("Using attribute orbitals " << norb)
	  } else {
	    cur = dset[i]->children;
	    while(cur != NULL) {
	      string cname((const char*)(cur->name));
	      if(cname == OrbitalBuilderBase::spo_tag) norb++;
	      cur=cur->next;
	    }
	    XMLReport("Counting the number or ortbials " << norb)
	  }
	  N.push_back(norb);
	}
	el.create(N);
	found_el=true;
      }
    }
    xmlXPathFreeObject(result);
    
    if(!found_el) {
      //if reached here, need to create the 2-electron system with a error message
      ERRORMSG("Cannot determine the number of electrons: assume 2-electron system")
      N.resize(2);
      N[0] = 1; N[1] = 1;
      el.create(N);    
    }
    return true;
  }
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

  
