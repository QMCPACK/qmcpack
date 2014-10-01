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

#include "Potential.h"

Potential* ReadPotential (IOSectionClass &in)
{
  string type;
  Potential *pot;
  assert (in.ReadVar ("Type", type));
  if (type == "Coulomb")
    pot = new CoulombPot;
  else if (type == "General")
    pot = new GeneralPot;
  else if (type == "Screened")
    pot = new ScreenedPot;
  else if (type == "Spline")
    pot = new SplinePot;
  else if (type == "NLPP")
    pot = new NLPPClass;
  else {
    cerr << "Unrecognize potential type \"" << type << "\".  Exitting.\n";
    exit(1);
  }
  pot->Read(in);
  return pot;
}
