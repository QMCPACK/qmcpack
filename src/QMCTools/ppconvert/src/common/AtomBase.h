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

#ifndef ATOM_BASE_H
#define ATOM_BASE_H

#include "RadialWF.h"

typedef enum {DFTType, HFType, OEPType} AtomType;

class Atom
{
protected:
  Grid *grid;
public:
  Array<RadialWF,1> RadialWFs;
  virtual AtomType Type() = 0;
  virtual void UpdateVHXC() = 0;
  virtual void CalcEnergies (double &kinetic, double &potential,
			     double &hartree, double &HXC) = 0;
  virtual void Solve() = 0;
  virtual void Write(IOSectionClass &out) = 0;
  virtual void Read(IOSectionClass &in) = 0;
  virtual void SetGrid(Grid *newGrid) = 0;
  inline Grid* GetGrid(){ return grid; }
  inline double NumElecs();
  virtual void SetBarePot (Potential *pot) = 0;
};

inline double
Atom::NumElecs()
{
  double num = 0.0;
  for (int i=0; i<RadialWFs.size(); i++) 
    num += RadialWFs(i).Occupancy;
  return num;
}


#endif
