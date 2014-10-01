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

#include "SplinePot.h"

double SplinePot::V(double r)
{
  if (r <= Spline.grid->End)
    return Spline(r);
  else if (Vouter != NULL) 
    return (Vouter->V(r));
  else {
#ifdef BZ_DEBUG
    cerr << "r outside grid in SplinePot:  " << r << endl;
#endif
    return 0.0;
  }
}

double SplinePot::dVdr(double r)
{
  if (r <= Spline.grid->End)
    return Spline.Deriv(r);
  else if (Vouter != NULL) 
    return (Vouter->dVdr(r));
  else {
#ifdef BZ_DEBUG
    cerr << "r outside grid in SplinePot:  " << r << endl;
#endif
    return 0.0;
  }
}

double SplinePot::d2Vdr2(double r)
{
  if (r <= Spline.grid->End)
    return Spline.Deriv2(r);
  else if (Vouter != NULL) 
    return (Vouter->d2Vdr2(r));
  else {
#ifdef BZ_DEBUG
    cerr << "r outside grid in SplinePot:  " << r << endl;
#endif
    return 0.0;
  }
}

void SplinePot::Read(IOSectionClass &in)
{
  assert(in.OpenSection("Grid"));
  Grid *grid = ReadGrid(in);
  in.CloseSection(); // "Grid" 
  Array<double,1> data;
  assert(in.ReadVar("SplineData", data));
  Spline.Init (grid, data);
  if(in.OpenSection("Vouter")) {
    Vouter = ReadPotential(in);
    in.CloseSection();
  }
  else
    Vouter = NULL;
 }


void SplinePot::Write(IOSectionClass &out)
{
  out.WriteVar ("Type", "Spline");
  out.NewSection("Grid");
  Spline.grid->Write(out);
  out.CloseSection();
  out.WriteVar ("SplineData", Spline.Data());
  if (Vouter != NULL) {
    out.NewSection ("Vouter");
    Vouter->Write(out);
    out.CloseSection();
  }
}
