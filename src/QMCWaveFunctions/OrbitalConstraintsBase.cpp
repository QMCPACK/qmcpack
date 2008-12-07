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

  OrbitalConstraintsBase::OrbitalConstraintsBase(ParticleSet& p, TrialWaveFunction& psi)
    : OrbitalBuilderBase(p,psi), myGrid(0), PrintTables(false) 
  {
  }

  void OrbitalConstraintsBase::getParam(xmlNodePtr cur) 
  {
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

  bool OrbitalConstraintsBase::getVariables(xmlNodePtr cur) 
  {

    xmlChar* prnode=xmlGetProp(cur,(const xmlChar*)"print");
    if(prnode != NULL)
    {
      PrintTables = xmlStrEqual(prnode,(const xmlChar*)"yes");
    }

    //save the xml node
    myNode=cur;
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

  void OrbitalConstraintsBase::setRadialGrid() 
  {
    if(myGrid) return;
    xmlNodePtr gridPtr=NULL;
    xmlNodePtr cur=myNode->children;
    while(cur != NULL) {
      if(xmlStrEqual(cur->name,(const xmlChar*)"grid")) {
        gridPtr=cur;
      }
      cur=cur->next;
    }

    if(gridPtr == NULL) {
      gridPtr = xmlNewNode(NULL,(const xmlChar*)"grid");
      xmlNewProp(gridPtr,(const xmlChar*)"type",(const xmlChar*)"log");
      xmlNewProp(gridPtr,(const xmlChar*)"ri",(const xmlChar*)"1.0e-6");
      std::ostringstream rf;
      if(targetPtcl.Lattice.SuperCellEnum)
        rf << targetPtcl.Lattice.LR_rc;
      else
        rf << 100.;
      xmlNewProp(gridPtr,(const xmlChar*)"rf",(const xmlChar*)rf.str().c_str());
      xmlNewProp(gridPtr,(const xmlChar*)"npts",(const xmlChar*)"101");
      xmlAddChild(cur,gridPtr);
    }
    //create grid and initialize CubicSplineFunctions
    myGrid = OneDimGridFactory::createGrid(gridPtr);
    //get the cutoff radius
    Rcut = OneDimGridFactory::setSmoothCutoff(myGrid,gridPtr);
    app_log() << "  Radial Grid Rcut=" << Rcut << " Rmax=" << myGrid->rmax() <<endl;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
