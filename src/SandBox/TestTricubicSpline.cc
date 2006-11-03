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

#include "TricubicSpline.h"

main()
{
  LinearGrid xgrid(0.0, 5.0, 20);
  LinearGrid ygrid(0.0, 5.0, 30);
  LinearGrid zgrid(0.0, 5.0, 40);
  Array<double,1> f(20,30,40);

  for (int ix=0; ix<xgrid.NumPoints; ix++)
    for (int iy=0; iy<xgrid.NumPoints; iy++)
      for (int iz=0; iz<xgrid.NumPoints; iz++) {
	x = xgrid(ix);	y = ygrid(iy);	z = zgrid(iz);
	f(ix,iy,iz) = cos(2.0*x)*cos(2.0*y)*cos(2.0*z);
      }
  
  TricubicSpline(&xgrid, &ygrid, &zgrid, f);
}

