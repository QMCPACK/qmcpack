//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
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

