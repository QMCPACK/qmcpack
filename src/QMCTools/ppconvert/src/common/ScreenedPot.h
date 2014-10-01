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

#ifndef SCREENED_POT_H
#define SCREENED_POT_H

#include "PotentialBase.h"
#include "CubicSplineCommon.h"

class ScreenedPot : public Potential
{
public:
  double Charge;
  CubicSplineCommon HXC;
  Potential *BarePot;

  bool IsPH();
  bool NeedsRel();
  double GetCoreRadius();
  double A      (double r);
  double B      (double r);
  double dAdr   (double r);
  double d2Adr2 (double r);

  double V      (int l, double r);
  double dVdr   (int l, double r);
  double d2Vdr2 (int l, double r);

  double V      (double r);
  double dVdr   (double r);
  double d2Vdr2 (double r);
  void   Write(IOSectionClass &out);
  void   Read (IOSectionClass &in);
};



#endif
