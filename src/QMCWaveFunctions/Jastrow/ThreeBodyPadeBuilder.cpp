//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
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
    std::string cname((const char*)(cur->name));
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
