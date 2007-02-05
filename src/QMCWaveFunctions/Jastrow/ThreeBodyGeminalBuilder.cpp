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
#include "QMCWaveFunctions/Jastrow/ThreeBodyGeminalBuilder.h"
#include "QMCWaveFunctions/Jastrow/JastrowBasisBuilder.h"
//#include "QMCWaveFunctions/MolecularOrbitals/GTOMolecularOrbitals.h"
//#include "QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.h"
#include "QMCWaveFunctions/Jastrow/ThreeBodyGeminal.h"
namespace qmcplusplus {

  ThreeBodyGeminalBuilder::ThreeBodyGeminalBuilder(ParticleSet& els, 
      TrialWaveFunction& wfs, 
      ParticleSet& ions):
    OrbitalBuilderBase(els,wfs) {
    //gtoBuilder = new GTOMolecularOrbitals(els,wfs,ions);
    //gtoBuilder = new GridMolecularOrbitals(els,wfs,ions);
    basisBuilder = new JastrowBasisBuilder(els,ions);
    J3 = new ThreeBodyGeminal(ions, els);
  }

  bool ThreeBodyGeminalBuilder::put(xmlNodePtr cur) {

    ThreeBodyGeminal::BasisSetType *basisSet=0;
    bool foundBasisSet=false;

    xmlNodePtr basisPtr=NULL;
    xmlNodePtr coeffPtr=NULL;
    cur = cur->xmlChildrenNode;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == basisset_tag) {
        basisPtr=cur;
        //call the BasisSet builder
        //basisSet = gtoBuilder->addBasisSet(cur);
      } else if(cname == "coefficient" || cname == "coefficients") {
        coeffPtr=cur;
      }
      cur=cur->next;
    }

    if(basisPtr != NULL)
    {
      foundBasisSet = basisBuilder->put(basisPtr);
      J3->setBasisSet(basisBuilder->myBasisSet);

      if(coeffPtr != NULL)
      {
        J3->put(coeffPtr,targetPsi.VarList);
      }
      else
      {
        cout << "Coefficients are not given." << endl;
      }
    }

    //add three-body jastrow
    targetPsi.addOrbital(J3);
    return true;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
