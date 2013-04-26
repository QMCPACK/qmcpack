//////////////////////////////////////////////////////////////////
// (c) Copyright 2008- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
  cout << "#### BsplineAOBuilder::addRadialOrbital " << endl;
  RadialOrbitalType* radorb= new RadialOrbitalType(0.0);
  radorb->put(cur);
  m_orbitals->Rnl.push_back(radorb);
  m_orbitals->RnlID.push_back(nlms);
  return true;
}

}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1237 $   $Date: 2006-08-08 11:52:16 -0500 (Tue, 08 Aug 2006) $
 * $Id: BsplineAOBuilder.cpp 1237 2006-08-08 16:52:16Z jnkim $
 ***************************************************************************/
