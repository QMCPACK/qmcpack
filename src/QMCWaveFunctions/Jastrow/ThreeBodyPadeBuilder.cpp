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
/**@file ThreeBodyPadeBuilder.cpp
 *@brief definition of three-body jastrow of Pade functions
 */
#include "QMCWaveFunctions/Jastrow/ThreeBodyPadeBuilder.h"
#include "QMCWaveFunctions/Jastrow/ThreeBodyPade.h"
namespace qmcplusplus
{

ThreeBodyPadeBuilder::ThreeBodyPadeBuilder(ParticleSet& els, TrialWaveFunction& wfs, ParticleSet& ions):
  OrbitalBuilderBase(els,wfs)
{
  J3 = new ThreeBodyPade(ions, els);
}

bool ThreeBodyPadeBuilder::put(xmlNodePtr cur)
{
  // ThreeBodyPade::BasisSetType *basisSet=0;
  cur = cur->xmlChildrenNode;
  while(cur != NULL)
  {
    string cname((const char*)(cur->name));
    /* if(cname == basisset_tag) {
       //call the BasisSet builder
       basisSet = gtoBuilder->addBasisSet(cur);
       if(!basisSet) return 0;
     } else*/
    if(cname == "coefficient" || cname == "coefficients")
    {
      //J3->setBasisSet(basisSet);
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
