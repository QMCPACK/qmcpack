//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCWaveFunctions/MolecularOrbitals/GTOBuilder.h"
namespace qmcplusplus
{

GTOBuilder::GTOBuilder(xmlNodePtr cur): Normalized(true), m_orbitals(0)
{
  if(cur != NULL)
  {
    putCommon(cur);
  }
}

bool GTOBuilder::putCommon(xmlNodePtr cur)
{
  const xmlChar* a=xmlGetProp(cur,(const xmlChar*)"normalized");
  if(a)
  {
    if(xmlStrEqual(a,(const xmlChar*)"no"))
      Normalized=false;
  }
  return true;
}

bool
GTOBuilder::addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms)
{
  if(!m_orbitals)
  {
    ERRORMSG("m_orbitals, SphericalOrbitals<ROT,GT>*, is not initialized")
    return false;
  }
  RadialOrbitalType* radorb= new RadialOrbitalType(nlms[q_l],Normalized);
  radorb->putBasisGroup(cur);
  m_orbitals->Rnl.push_back(radorb);
  m_orbitals->RnlID.push_back(nlms);
  return true;
}

}
