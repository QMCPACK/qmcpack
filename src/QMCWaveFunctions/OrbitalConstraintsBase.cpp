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
#include "QMCWaveFunctions/OrbitalConstraintsBase.h"

namespace qmcplusplus {

  void OrbitalConstraintsBase::getParam(xmlNodePtr cur) {
    const xmlChar* a=xmlGetProp(cur,(const xmlChar*)"id");
    const xmlChar* b=xmlGetProp(cur,(const xmlChar*)"name");
    if(a == NULL || b == NULL) return;
    RealType val;
    string vname((const char*)b);
    putContent(val,cur);
    map<string,pair<string,RealType> >::iterator vit(inVars.find(vname));
    if(vit == inVars.end()) {
      inVars[vname]=pair<string,RealType>((const char*)a,val);
    }
  }

  bool OrbitalConstraintsBase::getVariables(xmlNodePtr cur) {
    cur = cur->children;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == "correlation") {
        xmlNodePtr cur1=cur->children;
        while(cur1 != NULL) {
          string cname1((const char*)(cur1->name));
          if(cname1 == "parameter") {
            getParam(cur1);
          }
          cur1=cur1->next;
        }
      } else if(cname == "parameter") {
        getParam(cur);
      }
      cur = cur->next;
    } // while cur
    return true;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
