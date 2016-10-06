//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "CoulombPot.h"

double CoulombPot::V (double r)
{
  return (Z1Z2 / r);
}

double CoulombPot::dVdr (double r)
{
  return (-Z1Z2/(r*r));
}

double CoulombPot::d2Vdr2 (double r)
{
  return (2.0*Z1Z2/(r*r*r));
}

void CoulombPot::Write(IOSectionClass &out)
{
  out.WriteVar("Type", "Coulomb");
  out.WriteVar("Z1Z2", Z1Z2);
}

void CoulombPot::Read(IOSectionClass &in)
{
  assert(in.ReadVar("Z1Z2", Z1Z2));
}

bool CoulombPot::NeedsRel()
{
  return true;
}

double
CoulombPot::X_k(double rcut, double k)
{
  return -4.0*M_PI*Z1Z2/(k*k)*cos(k*rcut);
}
