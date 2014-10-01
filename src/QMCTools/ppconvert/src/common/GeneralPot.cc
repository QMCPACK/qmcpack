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

#include "GeneralPot.h"

void
GeneralPot::Read(IOSectionClass &in)
{
  assert (in.OpenSection("Grid"));
  PotGrid = ReadGrid (in);
  in.CloseSection();
  Array<double,1> vdata;
  assert (in.ReadVar("V", vdata));
  assert (in.ReadVar("Z", Z));
  PotSpline.Init (PotGrid, vdata);
}

void
GeneralPot::Write(IOSectionClass &out)
{
  out.NewSection("Grid");
  PotGrid->Write(out);
  out.CloseSection();
  out.WriteVar("V", PotSpline.Data());
  out.WriteVar("Z", Z);
  out.WriteVar("Type", "General");
}

double
GeneralPot::V(double r)
{
  if (r < PotGrid->End)
    return PotSpline(r);
  else
    return (-Z/r);
}

double
GeneralPot::dVdr (double r)
{
  if (r < PotGrid->End)
    return PotSpline.Deriv(r);
  else
    return Z/(r*r);
}

double
GeneralPot::d2Vdr2 (double r)
{
  if (r < PotGrid->End)
    return PotSpline.Deriv2(r);
  else
    return -2.0*Z/(r*r*r);
}


GeneralPot::GeneralPot() : PotGrid(NULL)
{

}

GeneralPot::~GeneralPot()
{
  if (PotGrid != NULL)
    delete PotGrid;
}
