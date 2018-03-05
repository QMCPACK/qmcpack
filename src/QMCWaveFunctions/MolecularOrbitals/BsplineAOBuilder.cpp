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


#include "QMCWaveFunctions/MolecularOrbitals/BsplineAOBuilder.h"
namespace qmcplusplus
{

BsplineAOBuilder::BsplineAOBuilder(xmlNodePtr cur): m_orbitals(0)
{
  if(cur != NULL)
    putCommon(cur);
}

bool BsplineAOBuilder::putCommon(xmlNodePtr cur)
{
  return true;
}

bool
BsplineAOBuilder::addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms, bool useSphericalHarmonicsNormalization)
{
  if(!m_orbitals)
  {
    ERRORMSG("m_orbitals, SphericalOrbitals<ROT,GT>*, is not initialized")
    return false;
  }
  std::cout << "#### BsplineAOBuilder::addRadialOrbital " << std::endl;
  RadialOrbitalType* radorb= new RadialOrbitalType(0.0);
  radorb->put(cur);
  m_orbitals->Rnl.push_back(radorb);
  m_orbitals->RnlID.push_back(nlms);
  return true;
}

}
