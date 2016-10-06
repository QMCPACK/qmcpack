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

#include "DFTAtom.h"
#include "RungeKutta.h"
#include "Functionals.h"

AtomType DFTAtom::Type()
{
  return (DFTType);
}

void DFTAtom::UpdateVHXC()
{
  UpdateChargeDensity();
  UpdateHartree();
  UpdateExCorr();

  for (int i=0; i < grid->NumPoints; i++)
    V.HXC(i) = Hartree(i) + ExCorr(i);
}

/// Radial WFs must be normalized before calling this function
void DFTAtom::UpdateChargeDensity()
{
  for (int i=0; i<grid->NumPoints; i++) {
    ChargeDensity(i) = 0.0;
    double r = (*grid)(i);
    double rinv2 = 1.0/(r*r);
    for (int j=0; j<RadialWFs.size(); j++) {
      double u = RadialWFs(j).u(i);
      ChargeDensity(i) += RadialWFs(j).Occupancy * u*u * rinv2;
    }
  }
}


class HartreeDeriv1
{
  DFTAtom &atom;
public:
  inline double operator()(double r, double sum)
  { return atom.Hartree1(r,sum); }

  HartreeDeriv1(DFTAtom &newAtom) : atom(newAtom) 
  { /* Do nothing */ }
};

class HartreeDeriv2
{
  DFTAtom &atom;
public:
  inline double operator()(double r, double sum)
  { return atom.Hartree2(r,sum); }

  HartreeDeriv2(DFTAtom &newAtom) : atom(newAtom) 
  { /* Do nothing */ }
};


void DFTAtom::UpdateHartree()
{
  int N = grid->NumPoints;
  
  HartreeDeriv1 H1(*this);
  HartreeDeriv2 H2(*this);
  RungeKutta<HartreeDeriv1,double> integrator1(H1);
  RungeKutta<HartreeDeriv2,double> integrator2(H2);
  temp(0) = 0.0;  temp2(0) = 0.0;
  integrator1.Integrate(*grid, 0, N-1, temp);
  integrator2.Integrate(*grid, 0, N-1, temp2);

  for (int i=0; i<N; i++) {
    double r = (*grid)(i);
    Hartree(i) = 4.0*M_PI*(temp(i)/r + (temp2(N-1)-temp2(i)));
  }
}

void DFTAtom::UpdateExCorr()
{
  int N = grid->NumPoints;
  for (int i=0; i<N; i++) {
    double Vxcup, Vxcdown;
    double upCharge = 0.5*ChargeDensity(i);
    double downCharge = upCharge;
    CPPExCorr(upCharge, downCharge, Vxcup, Vxcdown);
    //    FortranExCorr(upCharge, downCharge, Vxcup, Vxcdown);
    ExCorr(i) = Vxcup;
  }
}


void DFTAtom::SetGrid(Grid *newgrid)
{
  grid = newgrid;
  int N = grid->NumPoints;
  temp.resize(N);
  temp2.resize(N);
  temp = 0.0;
  V.HXC.Init(grid, temp);
  ChargeDensity.Init(grid,temp);
  Hartree.Init(grid, temp);
  ExCorr.Init(grid,temp);
  for (int i=0; i<RadialWFs.size(); i++)
    RadialWFs(i).SetGrid(grid);
}

void DFTAtom::SetBarePot(Potential *newPot)
{
  BarePot = newPot;
  V.BarePot = BarePot;
  for (int i=0; i<RadialWFs.size(); i++)
    RadialWFs(i).SetPotential(&V);
}


void
DFTAtom::SolveInit()
{
  int N = grid->NumPoints; 
  // Just rename temp and temp2 for clarity
  Array<double,1> &oldCharge = temp;
  Array<double,1> &newCharge = temp2;

  // First, zero out screening
  for (int i=0; i<N; i++)
    V.HXC(i) = 0.0;

  OldEnergies.resize(RadialWFs.size());  
  // Now solve radial equations
  for (int i=0; i<RadialWFs.size(); i++) {
    RadialWFs(i).Solve();
    //fprintf (stderr, "Energy(%d) = %1.16f\n", i, RadialWFs(i).Energy);
    OldEnergies(i) = RadialWFs(i).Energy;
  }

  if (NewMix < 1.0e-8) // We don't want to do self-consistent
    return;

  UpdateChargeDensity();
  oldCharge = 0.0;
  newCharge = ChargeDensity.Data();
}

double
DFTAtom::SolveIter()
{
  int N = grid->NumPoints; 
  // Just rename temp and temp2 for clarity
  Array<double,1> &oldCharge = temp;
  Array<double,1> &newCharge = temp2;

  for (int i=0; i<N; i++)
    ChargeDensity(i) = NewMix*(newCharge(i)) + (1.0-NewMix)*oldCharge(i);
  UpdateHartree();
  UpdateExCorr();
  for (int i=0; i < grid->NumPoints; i++)
    V.HXC(i) = Hartree(i) + ExCorr(i);
  
  double maxDiff = 0.0;
  for (int i=0; i<RadialWFs.size(); i++) {
    RadialWFs(i).Solve();
    maxDiff = max(fabs(OldEnergies(i)-RadialWFs(i).Energy), maxDiff);
    OldEnergies(i) = RadialWFs(i).Energy;
    // fprintf (stderr, "Energy(%d) = %1.16f\n", i, RadialWFs(i).Energy);
  }
  oldCharge = ChargeDensity.Data();
  UpdateChargeDensity();
  newCharge = ChargeDensity.Data();
  return (maxDiff);
}


void DFTAtom::Solve()
{
  int N = grid->NumPoints; 
  // Just rename temp and temp2 for clarity
  Array<double,1> &oldCharge = temp;
  Array<double,1> &newCharge = temp2;

  // First, zero out screening
  for (int i=0; i<N; i++)
    V.HXC(i) = 0.0;

  Array<double,1> oldEnergies(RadialWFs.size());  
  // Now solve radial equations
  for (int i=0; i<RadialWFs.size(); i++) {
    RadialWFs(i).Solve();
    // fprintf (stderr, "Energy(%d) = %1.16f\n", i, RadialWFs(i).Energy);
    oldEnergies(i) = RadialWFs(i).Energy;
  }

  if (NewMix < 1.0e-8) // We don't want to do self-consistent
    return;

  UpdateChargeDensity();
  oldCharge = 0.0;
  newCharge = ChargeDensity.Data();

  bool done = false;

  while (!done) {
    for (int i=0; i<N; i++)
      ChargeDensity(i) = NewMix*(newCharge(i)) + (1.0-NewMix)*oldCharge(i);
    UpdateHartree();
    UpdateExCorr();
    for (int i=0; i < grid->NumPoints; i++)
      V.HXC(i) = Hartree(i) + ExCorr(i);

    done = true;
    for (int i=0; i<RadialWFs.size(); i++) {
      RadialWFs(i).Solve();
      if (fabs(oldEnergies(i)-RadialWFs(i).Energy) > 1.0e-8)
	done = false;
      oldEnergies(i) = RadialWFs(i).Energy;
    }
    for (int i=0; i<RadialWFs.size(); i++)
      fprintf (stderr, "Energy(%d) = %1.16f\n", i, RadialWFs(i).Energy);
    oldCharge = ChargeDensity.Data();
    UpdateChargeDensity();
    newCharge = ChargeDensity.Data();
  }
  

}


void DFTAtom::Write(IOSectionClass &out)
{
  out.WriteVar ("Type", "DFT");
  out.NewSection("Grid");
  grid->Write(out);
  out.CloseSection();

  for (int i=0; i<RadialWFs.size(); i++) {
    out.NewSection("RadialWF");
    RadialWFs(i).Write(out);
    out.CloseSection();
  }

  out.WriteVar ("ChargeDensity", ChargeDensity.Data());
  out.WriteVar ("Hartree", Hartree.Data());
  out.WriteVar ("ExCorr", ExCorr.Data());
  out.NewSection("Potential");
  BarePot->Write(out);
  out.CloseSection();
  out.WriteVar ("NewMix", NewMix);
}


void DFTAtom::Read(IOSectionClass &in) 
{
  assert(in.OpenSection("Grid"));
  grid = ReadGrid(in);
  in.CloseSection();

  assert(in.OpenSection("Potential"));
  BarePot = ReadPotential(in);
  in.CloseSection();
  V.BarePot = BarePot;

  int numRadialWFs = in.CountSections("RadialWF");
  RadialWFs.resize(numRadialWFs);
  SetGrid(grid);
  in.ReadVar ("ChargeDensity", ChargeDensity.Data());
  in.ReadVar ("Hartree", Hartree.Data());
  in.ReadVar ("ExCorr", ExCorr.Data());
  for (int i=0; i<V.HXC.size(); i++)
    V.HXC(i) = Hartree(i)+ExCorr(i);

  double charge = 0.0;
  for (int i=0; i<numRadialWFs; i++) {
    RadialWFs(i).SetGrid (grid);
    RadialWFs(i).SetPotential (&V);
    assert (in.OpenSection("RadialWF", i));
    RadialWFs(i).Read(in);
    charge += RadialWFs(i).Occupancy;
    in.CloseSection();
  }
  V.Charge = charge;
  assert (in.ReadVar("NewMix", NewMix));
}


void DFTAtom::CalcEnergies(double &kinetic, double &potential,
			   double &hartree, double &XC)
{

}
