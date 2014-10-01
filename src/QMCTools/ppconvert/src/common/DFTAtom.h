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

#ifndef DFT_ATOM_H
#define DFT_ATOM_H

#include "AtomBase.h"

class DFTAtom : public Atom
{
private:
  CubicSplineCommon ChargeDensity;
  CubicSplineCommon Hartree;
  CubicSplineCommon ExCorr;
  void UpdateChargeDensity();
  void UpdateHartree();
  void UpdateExCorr();
  Array<double,1> temp, temp2;
  Potential *BarePot;
  Array<double,1> OldEnergies;
public:
  ScreenedPot V;
  double NewMix;
  
  AtomType Type();
  /// This function calculates the charge density, hartree and exchange
  /// potentials and places them in pot.
  void UpdateVHXC();
  void CalcEnergies (double &kinetic, double &potential, 
		     double &hartree, double &XC);
  void Solve();
  void SolveInit();
  double SolveIter();
  void Write (IOSectionClass &out);
  void Read  (IOSectionClass &in);
  void SetGrid (Grid *newGrid);
  void SetBarePot (Potential *pot);

  inline double rho (double r) { return ChargeDensity(r); }
  
  inline double Hartree1 (double r, double sum);
  inline double Hartree2 (double r, double sum);
};

inline double DFTAtom::Hartree1 (double r, double sum)
{
  double rho = ChargeDensity(r);
  return (rho*r*r);
}

inline double DFTAtom::Hartree2 (double r, double sum)
{
  double rho = ChargeDensity(r);
  return (rho*r);
}

#endif
