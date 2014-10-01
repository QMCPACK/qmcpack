/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
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
