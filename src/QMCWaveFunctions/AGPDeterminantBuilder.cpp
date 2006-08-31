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
#include "QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.h"
//#include "QMCWaveFunctions/MolecularOrbitals/STOMolecularOrbitals.h"
//#include "QMCWaveFunctions/MolecularOrbitals/GTOMolecularOrbitals.h"
#include "QMCWaveFunctions/AGPDeterminantBuilder.h"
namespace qmcplusplus {

  AGPDeterminantBuilder::AGPDeterminantBuilder(ParticleSet& els, 
      TrialWaveFunction& wfs, ParticleSet& ions): 
    OrbitalBuilderBase(els,wfs), ionRef(ions), basisBuilder(0), agpDet(0) {
  }

  template<class BasisBuilderT>
  bool AGPDeterminantBuilder::createAGP(BasisBuilderT *abuilder, xmlNodePtr cur) {
    bool spinpolarized=false;
    typename BasisBuilderT::BasisSetType *basisSet=0;
    cur = cur->xmlChildrenNode;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == basisset_tag) {
        basisSet = abuilder->addBasisSet(cur);
        if(!basisSet) return false;
      } else if(cname == "coefficient" || cname == "coefficients") {
        if(agpDet == 0) {
          int nup=targetPtcl.first(1), ndown=0;
          if(targetPtcl.groups()>1) ndown = targetPtcl.first(2)-nup;
          basisSet->resize(nup+ndown);
          agpDet = new AGPDeterminant(basisSet);
          agpDet->resize(nup,ndown);
        }
        int offset=1;
        xmlNodePtr tcur=cur->xmlChildrenNode;
        while(tcur != NULL) {
          if(xmlStrEqual(tcur->name,(const xmlChar*)"lambda")) {
            int i=atoi((const char*)(xmlGetProp(tcur,(const xmlChar*)"i")));
            int j=atoi((const char*)(xmlGetProp(tcur,(const xmlChar*)"j")));
            double c=atof((const char*)(xmlGetProp(tcur,(const xmlChar*)"c")));
            agpDet->Lambda(i-offset,j-offset)=c;
            if(i != j) {
              agpDet->Lambda(j-offset,i-offset)=c;
            }
          }
          tcur=tcur->next;
        }
      } else if(cname == "unpaired") {
        spinpolarized=true;
        int offset=1;
        xmlNodePtr tcur=cur->xmlChildrenNode;
        while(tcur != NULL) {
          if(xmlStrEqual(tcur->name,(const xmlChar*)"lambda")) {
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
    if(spinpolarized) {
      app_log() << "  Coefficients for the unpaired electrons " << endl;
      app_log() << agpDet->LambdaUP << endl;
    }
    return true;
  }

  bool AGPDeterminantBuilder::put(xmlNodePtr cur) {

    bool success = false;
    if(basisBuilder == 0) {
      GridMolecularOrbitals* gtoBuilder = new GridMolecularOrbitals(targetPtcl,targetPsi,ionRef);
      //STOMolecularOrbitals* gtoBuilder = new STOMolecularOrbitals(targetPtcl,targetPsi,ionRef);
      //GTOMolecularOrbitals* gtoBuilder = new GTOMolecularOrbitals(targetPtcl,targetPsi,ionRef);
      success = createAGP(gtoBuilder,cur);
      basisBuilder=gtoBuilder;
    }

    if(agpDet) targetPsi.addOrbital(agpDet);
    return success;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
