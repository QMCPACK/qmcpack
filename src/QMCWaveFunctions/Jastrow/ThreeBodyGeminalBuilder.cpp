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
namespace qmcplusplus
{

ThreeBodyGeminalBuilder::ThreeBodyGeminalBuilder(ParticleSet& els,
    TrialWaveFunction& wfs,
    ParticleSet& ions):
  OrbitalBuilderBase(els,wfs), sourcePtcl(ions), basisBuilder(0)
{
}

bool ThreeBodyGeminalBuilder::put(xmlNodePtr cur)
{
  ThreeBodyGeminal::BasisSetType *basisSet=0;
  bool foundBasisSet=false;
  xmlNodePtr basisPtr=NULL;
  xmlNodePtr coeffPtr=NULL;
  cur = cur->xmlChildrenNode;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == basisset_tag)
    {
      basisPtr=cur;
    }
    else
      if(cname == "coefficient" || cname == "coefficients")
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
