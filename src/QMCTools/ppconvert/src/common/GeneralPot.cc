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
