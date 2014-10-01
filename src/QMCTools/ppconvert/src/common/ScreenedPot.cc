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

#include "ScreenedPot.h"

bool ScreenedPot::IsPH()
{
  return (BarePot->IsPH());
}

bool ScreenedPot::NeedsRel()
{
  return (BarePot->NeedsRel());
}

double ScreenedPot::GetCoreRadius()
{
  return (BarePot->GetCoreRadius());
}

double ScreenedPot::A (double r)
{
  return (BarePot->A(r));
}

double ScreenedPot::B (double r)
{
  return (BarePot->B(r));
}

double ScreenedPot::dAdr (double r)
{
  return (BarePot->dAdr(r));
}

double ScreenedPot::d2Adr2 (double r)
{
  return (BarePot->d2Adr2(r));
}


double ScreenedPot::V (double r)
{
  double VHXC;
  VHXC = (r <= HXC.grid->End) ? HXC(r) : Charge/r;
  return (VHXC + BarePot->V(r));
}


double ScreenedPot::dVdr (double r)
{
  double dVHXC;
  dVHXC = (r <= HXC.grid->End) ? HXC.Deriv(r) : -Charge/(r*r);
  return (dVHXC + BarePot->dVdr(r));
}

double ScreenedPot::d2Vdr2(double r)
{
  double d2VHXC;
  d2VHXC = (r <= HXC.grid->End) ? HXC.Deriv2(r) : 2.0*Charge/(r*r*r);
  return (d2VHXC + BarePot->d2Vdr2(r));
}


double ScreenedPot::V (int l, double r)
{
  double VHXC;
  VHXC = (r <= HXC.grid->End) ? HXC(r) : Charge/r;
  return (VHXC + BarePot->V(l, r));
}


double ScreenedPot::dVdr (int l, double r)
{
  double dVHXC;
  dVHXC = (r <= HXC.grid->End) ? HXC.Deriv(r) : -Charge/(r*r);
  return (dVHXC + BarePot->dVdr(l, r));
}

double ScreenedPot::d2Vdr2(int l, double r)
{
  double d2VHXC;
  d2VHXC = (r <= HXC.grid->End) ? HXC.Deriv2(r) : 2.0*Charge/(r*r*r);
  return (d2VHXC + BarePot->d2Vdr2(l, r));
}


void ScreenedPot::Write (IOSectionClass &out)
{
  out.WriteVar ("Type", "Screened");
  out.WriteVar ("Charge", Charge);

  out.NewSection("HXC");
  HXC.Write (out);
  out.CloseSection();

  out.NewSection ("BarePot");
  BarePot->Write(out);
  out.CloseSection();
}


void ScreenedPot::Read (IOSectionClass &in)
{
  assert (in.ReadVar ("Charge", Charge));

  assert (in.OpenSection("HXC"));
  HXC.Read(in);
  in.CloseSection();

  assert (in.OpenSection("BarePot"));
  BarePot=ReadPotential(in);
  in.CloseSection();
}

