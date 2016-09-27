//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//                    Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCWaveFunctions/OrbitalConstraintsBase.h"

namespace qmcplusplus
{

OrbitalConstraintsBase::OrbitalConstraintsBase(ParticleSet& p, TrialWaveFunction& psi)
  : OrbitalBuilderBase(p,psi), myGrid(0), PrintTables(false)
{
}

void OrbitalConstraintsBase::getParam(xmlNodePtr cur)
{
  const xmlChar* a=xmlGetProp(cur,(const xmlChar*)"id");
  const xmlChar* b=xmlGetProp(cur,(const xmlChar*)"name");
  if(a == NULL || b == NULL)
    return;
  RealType val;
  std::string vname((const char*)b);
  putContent(val,cur);
  std::map<std::string,std::pair<std::string,RealType> >::iterator vit(inVars.find(vname));
  if(vit == inVars.end())
  {
    inVars[vname]= std::pair<std::string,RealType>((const char*)a,val);
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
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "correlation")
    {
      xmlNodePtr cur1=cur->children;
      while(cur1 != NULL)
      {
        std::string cname1((const char*)(cur1->name));
        if(cname1 == "parameter")
        {
          getParam(cur1);
        }
        cur1=cur1->next;
      }
    }
    else
      if(cname == "parameter")
      {
        getParam(cur);
      }
    cur = cur->next;
  } // while cur
  return true;
}

void OrbitalConstraintsBase::setRadialGrid()
{
  if(myGrid)
    return;
  xmlNodePtr gridPtr=NULL;
  xmlNodePtr cur=myNode->children;
  while(cur != NULL)
  {
    if(xmlStrEqual(cur->name,(const xmlChar*)"grid"))
    {
      gridPtr=cur;
    }
    cur=cur->next;
  }
  if(gridPtr == NULL)
  {
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
  app_log() << "  Radial Grid Rcut=" << Rcut << " Rmax=" << myGrid->rmax() << std::endl;
}
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
