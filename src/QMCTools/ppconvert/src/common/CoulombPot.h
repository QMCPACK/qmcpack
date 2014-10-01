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

#ifndef COULOMB_POT_H
#define COULOMB_POT_H

#include "PotentialBase.h"

class CoulombPot : public Potential
{
public:
  double Z1Z2;
  
  double V      (double r);
  double dVdr   (double r);
  double d2Vdr2 (double r);
  double X_k (double rcut, double k);
  bool NeedsRel();
  void Write(IOSectionClass &out);
  void Read (IOSectionClass &in);
};

#endif
