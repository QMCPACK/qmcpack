//////////////////////////////////////////////////////////////////
// (c) Copyright 2005-  by Jeongnim Kim
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
/**@file ThreeBodyGeminalBuilder.cpp
 *@brief definition of three-body jastrow of Geminal functions
 */
#include "QMCWaveFunctions/AGPDeterminant.h"
#include "QMCWaveFunctions/AGPDeterminantBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/MolecularBasisBuilder.h"
//#include "QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.h"
//#include "QMCWaveFunctions/MolecularOrbitals/STOMolecularOrbitals.h"
//#include "QMCWaveFunctions/MolecularOrbitals/GTOMolecularOrbitals.h"
namespace qmcplusplus
{

AGPDeterminantBuilder::AGPDeterminantBuilder(ParticleSet& els, TrialWaveFunction& wfs,
    PtclPoolType& pset):
  OrbitalBuilderBase(els,wfs), ptclPool(pset), myBasisSetFactory(0), agpDet(0)
{
}

template<class BasisBuilderT>
bool AGPDeterminantBuilder::createAGP(BasisBuilderT *abuilder, xmlNodePtr cur)
{
  bool spinpolarized=false;
  typename BasisBuilderT::BasisSetType *basisSet=0;
  cur = cur->xmlChildrenNode;
  while(cur != NULL)
  {
    string cname((const char*)(cur->name));
    if(cname == basisset_tag)
    {
      basisSet = abuilder->addBasisSet(cur);
      if(!basisSet)
        return false;
    }
    else if(cname == "coefficient" || cname == "coefficients")
    {
      if(agpDet == 0)
      {
        int nup=targetPtcl.first(1), ndown=0;
        if(targetPtcl.groups()>1)
          ndown = targetPtcl.first(2)-nup;
        basisSet->resize(nup+ndown);
        agpDet = new AGPDeterminant(basisSet);
        agpDet->resize(nup,ndown);
      }
      int offset=1;
      xmlNodePtr tcur=cur->xmlChildrenNode;
      while(tcur != NULL)
      {
        if(xmlStrEqual(tcur->name,(const xmlChar*)"lambda"))
        {
          int i=atoi((const char*)(xmlGetProp(tcur,(const xmlChar*)"i")));
          int j=atoi((const char*)(xmlGetProp(tcur,(const xmlChar*)"j")));
          double c=atof((const char*)(xmlGetProp(tcur,(const xmlChar*)"c")));
          agpDet->Lambda(i-offset,j-offset)=c;
          if(i != j)
          {
            agpDet->Lambda(j-offset,i-offset)=c;
          }
        }
        tcur=tcur->next;
      }
    }
    else if(cname == "unpaired")
    {
      spinpolarized=true;
      int offset=1;
      xmlNodePtr tcur=cur->xmlChildrenNode;
      while(tcur != NULL)
      {
        if(xmlStrEqual(tcur->name,(const xmlChar*)"lambda"))
        {
          int i=atoi((const char*)(xmlGetProp(tcur,(const xmlChar*)"i")));
          int j=atoi((const char*)(xmlGetProp(tcur,(const xmlChar*)"j")));
          double c=atof((const char*)(xmlGetProp(tcur,(const xmlChar*)"c")));
          agpDet->LambdaUP(j-offset,i-offset)=c;
        }
        tcur=tcur->next;
      }
    }
    cur=cur->next;
  }
  //app_log() << agpDet->Lambda << endl;
  if(spinpolarized)
  {
    app_log() << "  Coefficients for the unpaired electrons " << endl;
    app_log() << agpDet->LambdaUP << endl;
  }
  return true;
}

bool AGPDeterminantBuilder::put(xmlNodePtr cur)
{
  if(agpDet)
  {
    app_error() << "  AGPDeterminantBuilder::put exits. AGPDeterminant has been already created."
                << endl;
    return false;
  }
  app_log() << "  AGPDeterminantBuilder Creating AGPDeterminant." << endl;
  xmlNodePtr curRoot=cur;
  bool success=true;
  string cname, tname;
  xmlNodePtr bPtr=NULL;
  xmlNodePtr cPtr=NULL;
  xmlNodePtr uPtr=NULL;
  OhmmsAttributeSet oAttrib;
  oAttrib.add(funcOpt,"function");
  oAttrib.add(transformOpt,"transform");
  oAttrib.put(cur);
  cur = cur->children;
  while(cur != NULL)
  {
    getNodeName(cname,cur);
    if(cname == basisset_tag)
    {
      bPtr=cur;
    }
    else if(cname.find("coeff")<cname.size())
    {
      cPtr=cur;
    }
    else if(cname.find("un")<cname.size())
    {
      uPtr=cur;
    }
    cur=cur->next;
  }
  if(bPtr == NULL || cPtr == NULL)
  {
    app_error() << "  AGPDeterminantBuilder::put Cannot create AGPDeterminant. " << endl;
    app_error() << "    Missing <basisset/> or <coefficients/>" << endl;
    return false;
  }
  if(myBasisSetFactory == 0)
  {
    myBasisSetFactory = new BasisSetFactory(targetPtcl,targetPsi,ptclPool);
    myBasisSetFactory->createBasisSet(bPtr,curRoot);
  }
// mmorales: this needs to be fixed after changes to BasisSetfactory
//    BasisSetBase<RealType>* myBasisSet=myBasisSetFactory->getBasisSet();
  BasisSetBase<RealType>* myBasisSet=0; //=myBasisSetFactory->getBasisSet();
  int nup=targetPtcl.first(1), ndown=0;
  if(targetPtcl.groups()>1)
    ndown = targetPtcl.first(2)-nup;
  myBasisSet->resize(nup+ndown);
  agpDet = new AGPDeterminant(myBasisSet);
  agpDet->resize(nup,ndown);
  //set Lambda: possible to move to AGPDeterminant
  int offset=1;
  xmlNodePtr tcur=cPtr->xmlChildrenNode;
  while(tcur != NULL)
  {
    if(xmlStrEqual(tcur->name,(const xmlChar*)"lambda"))
    {
      int i=atoi((const char*)(xmlGetProp(tcur,(const xmlChar*)"i")));
      int j=atoi((const char*)(xmlGetProp(tcur,(const xmlChar*)"j")));
      double c=atof((const char*)(xmlGetProp(tcur,(const xmlChar*)"c")));
      agpDet->Lambda(i-offset,j-offset)=c;
      if(i != j)
      {
        agpDet->Lambda(j-offset,i-offset)=c;
      }
    }
    tcur=tcur->next;
  }
  if(uPtr != NULL)
  {
    tcur=uPtr->xmlChildrenNode;
    while(tcur != NULL)
    {
      if(xmlStrEqual(tcur->name,(const xmlChar*)"lambda"))
      {
        int i=atoi((const char*)(xmlGetProp(tcur,(const xmlChar*)"i")));
        int j=atoi((const char*)(xmlGetProp(tcur,(const xmlChar*)"j")));
        double c=atof((const char*)(xmlGetProp(tcur,(const xmlChar*)"c")));
        agpDet->LambdaUP(j-offset,i-offset)=c;
      }
      tcur=tcur->next;
    }
    app_log() << "  AGPDeterminantBuilder::put Coefficients for the unpaired electrons " << endl;
    app_log() << agpDet->LambdaUP << endl;
  }
  if(agpDet)
    targetPsi.addOrbital(agpDet,"AGP");
  return success;
}
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
