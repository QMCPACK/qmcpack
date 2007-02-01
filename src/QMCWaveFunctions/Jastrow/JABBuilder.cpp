//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCWaveFunctions/Jastrow/JABBuilder.h"
#include "QMCWaveFunctions/Jastrow/PadeJastrow.h"
#include "QMCWaveFunctions/Jastrow/NoCuspJastrow.h"
#include "QMCWaveFunctions/Jastrow/ModPadeJastrow.h"
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowFunction.h"

namespace qmcplusplus {

  template<class FN>
  bool JABBuilder::createJAB(xmlNodePtr cur, FN* dummy) {

    int cur_var = targetPsi.VarList.size();
    string corr_tag("correlation");

    vector<FN*> jastrow;
    int ng = 0;
    ParticleSet* sourcePtcl=0;
    const xmlChar* s=xmlGetProp(cur,(const xmlChar*)"source");
    if(s != NULL) {
      map<string,ParticleSet*>::iterator pa_it(ptclPool.find((const char*)s));
      if(pa_it == ptclPool.end()) return false;
      sourcePtcl = (*pa_it).second;
      ng=sourcePtcl->getSpeciesSet().getTotalNum();
      for(int i=0; i<ng; i++) jastrow.push_back(0);
    }

    cur = cur->xmlChildrenNode;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == dtable_tag) {
      	string source_name((const char*)(xmlGetProp(cur,(const xmlChar *)"source")));
        map<string,ParticleSet*>::iterator pa_it(ptclPool.find(source_name));
        if(pa_it == ptclPool.end()) return false;
        sourcePtcl = (*pa_it).second;
        ng=sourcePtcl->getSpeciesSet().getTotalNum();
	XMLReport("Number of sources " << ng)
	for(int i=0; i<ng; i++) jastrow.push_back(0);
      } else if(cname == corr_tag) {
        if(sourcePtcl == 0) return false;
	string speciesA((const char*)(xmlGetProp(cur,(const xmlChar *)"speciesA")));
	//string speciesB((const char*)(xmlGetProp(cur,(const xmlChar *)"speciesB")));
	int ia = sourcePtcl->getSpeciesSet().findSpecies(speciesA);
	if(!(jastrow[ia])) {
	  jastrow[ia]= new FN;
	  jastrow[ia]->put(cur);
          jastrow[ia]->addOptimizables(targetPsi.VarList);
	  LOGMSG("  Added Jastrow Correlation between " << speciesA)
	}
      }
      cur = cur->next;
    } // while cur
    
    if(sourcePtcl == 0) {//invalid input: cleanup
      for(int ig=0; ig<ng; ig++)  if(jastrow[ig]) delete jastrow[ig];
      return false;
    }

    typedef OneBodyJastrow<FN> JneType;
    JneType* J1 = new JneType(*sourcePtcl,targetPtcl);
    for(int ig=0; ig<ng; ig++) {
      J1->addFunc(ig,jastrow[ig]);
    }

    //set this jastrow function to be not optimizable
    if(targetPsi.VarList.size() == cur_var) {
      J1->setOptimizable(false);
    }

    targetPsi.addOrbital(J1);
    XMLReport("Added a One-Body Jastrow Function")
    return true;
  }

  bool JABBuilder::put(xmlNodePtr cur) {

    string jastfunction("pade");
    const xmlChar *ftype = xmlGetProp(cur, (const xmlChar *)"function");
    if(ftype) jastfunction = (const char*) ftype;

    bool success=false;
    app_log() << "  One-Body Jastrow Function = " << jastfunction << endl;
    if(jastfunction == "nocusp") {
      NoCuspJastrow<RealType> *dummy = 0;
      success = createJAB(cur,dummy);
    } else if(jastfunction == "pade") {
      PadeJastrow<RealType> *dummy = 0;
      success = createJAB(cur,dummy);
    } else if(jastfunction == "pade2") {
      PadeJastrow2<RealType> *dummy = 0;
      success = createJAB(cur,dummy);
    } else if(jastfunction == "short") {
      ModPadeJastrow<RealType> *dummy = 0;
      success = createJAB(cur,dummy);
    }
    return success;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
