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
