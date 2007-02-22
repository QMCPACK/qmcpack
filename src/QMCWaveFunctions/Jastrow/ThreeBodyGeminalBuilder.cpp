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
#include "QMCWaveFunctions/MolecularOrbitals/NGOBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/GTOBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/STOBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/MolecularBasisBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/GTOMolecularOrbitals.h"
#include "QMCWaveFunctions/Jastrow/ThreeBodyGeminal.h"
namespace qmcplusplus {

  ThreeBodyGeminalBuilder::ThreeBodyGeminalBuilder(ParticleSet& els, 
      TrialWaveFunction& wfs, 
      ParticleSet& ions):
    OrbitalBuilderBase(els,wfs), sourcePtcl(ions), basisBuilder(0)
    {
    }

  bool ThreeBodyGeminalBuilder::put(xmlNodePtr cur) {

    ThreeBodyGeminal::BasisSetType *basisSet=0;
    bool foundBasisSet=false;

    xmlNodePtr basisPtr=NULL;
    xmlNodePtr coeffPtr=NULL;
    cur = cur->xmlChildrenNode;
    while(cur != NULL) 
    {
      string cname((const char*)(cur->name));
      if(cname == basisset_tag) 
      {
        basisPtr=cur;
      } 
      else if(cname == "coefficient" || cname == "coefficients") 
      {
        coeffPtr=cur;
      }
      cur=cur->next;
    }

    if(basisPtr != NULL)
    {
      ThreeBodyGeminal* J3 = new ThreeBodyGeminal(sourcePtcl, targetPtcl);
      basisBuilder = new JastrowBasisBuilder(targetPtcl,sourcePtcl,"gto",false);
      basisBuilder->put(basisPtr);
      J3->setBasisSet(basisBuilder->myBasisSet);
      J3->put(coeffPtr,targetPsi.VarList);
      targetPsi.addOrbital(J3);
      return true;
    }
    return false;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
