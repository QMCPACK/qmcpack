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
#include "QMCWaveFunctions/ThreeBodyGeminalBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/GTOMolecularOrbitals.h"
//#include "QMCWaveFunctions/MolecularOrbitals/GridMolecularOrbitals.h"
#include "QMCWaveFunctions/ThreeBodyGeminal.h"
namespace qmcplusplus {

  ThreeBodyGeminalBuilder::ThreeBodyGeminalBuilder(ParticleSet& els, 
      TrialWaveFunction& wfs, 
      ParticleSet& ions):
    OrbitalBuilderBase(els,wfs) {
    gtoBuilder = new GTOMolecularOrbitals(els,wfs,ions);
    //gtoBuilder = new GridMolecularOrbitals(els,wfs,ions);
    J3 = new ThreeBodyGeminal(ions, els);
  }

  bool ThreeBodyGeminalBuilder::put(xmlNodePtr cur) {

    ThreeBodyGeminal::BasisSetType *basisSet=0;
    cur = cur->xmlChildrenNode;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == basisset_tag) {
        //call the BasisSet builder
        basisSet = gtoBuilder->addBasisSet(cur);
        if(!basisSet) return 0;
      } else if(cname == "coefficient" || cname == "coefficients") {
        J3->setBasisSet(basisSet);
        J3->put(cur,targetPsi.VarList);
      }
      cur=cur->next;
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
