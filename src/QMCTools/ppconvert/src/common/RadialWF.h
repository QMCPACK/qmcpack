//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#ifndef RADIAL_WF_H
#define RADIAL_WF_H

#include "Potential.h"
#include "CubicSplineCommon.h"
#include "IO.h"

class RadialWF
{
private:
  int TurningIndex();
  double IntegrateInOut(int &tindex);
  void OriginBC(double r0, double &u0, double &du0);
  void InfinityBC(double rend, double &uend, double &duend);
  Array<Vec2,1> uduVec;
  Array<double,1> normVec;
  Potential *pot;  
public:
  Grid *grid;
  CubicSplineCommon u, dudr;
  int n, l, CoreNodes;
  double Energy, Occupancy, Weight;
  std::string Label;

  inline double NormDeriv(double r, double u);
  inline Vec2 PseudoDerivs (double r, Vec2 &u_and_du);
  inline Vec2 NormalDerivs (double r, Vec2 &u_and_du);
  inline Vec2 ScalarRelDerivs (double r, Vec2 &u_and_du);
  int CountNodes();
  void IntegrateOut();
  double PartialNorm();
  double LogDerivative();
  void Solve (double tolerance=1.0e-8);
  void Normalize();
  void SetGrid(Grid *newgrid);
  void SetPotential (Potential *newPot);
  Potential *GetPotential ();
  void Write (IOSectionClass &out);
  void Read  (IOSectionClass &in);
};



// These are wrappers so we can inline integration routines

/// Derivative used for normalization of the radial function
class NormalizeDeriv
{
private:
  RadialWF &WF;
public:
  inline double operator()(double r, double u)
  { return WF.NormDeriv(r, u); }
  NormalizeDeriv (RadialWF &wf) : WF(wf)
  { /* Do nothing */ }
};

/// Derivatives for the radial equation given a pseudoHamiltonian  
class PHDerivs 
{
private:
  RadialWF &WF;
public:
  inline Vec2 operator() (double r, Vec2& u_and_du)
  { return WF.PseudoDerivs(r, u_and_du); }
  PHDerivs (RadialWF &wf) : WF(wf)
  { /* Do nothing */ }
};

/// Derivatives for the radial equation given a pseudoHamiltonian  
class RegularDerivs 
{
private:
  RadialWF &WF;
public:
  inline Vec2 operator() (double r, Vec2& u_and_du)
  { return WF.NormalDerivs(r, u_and_du); }
  RegularDerivs (RadialWF &wf) : WF(wf)
  { /* Do nothing */ }
};



/// Derivatives for the scalar relativistic radial equation
class NonPHDerivs
{
private:
  RadialWF &WF;
public:
  inline Vec2 operator() (double r, Vec2& u_and_du)
  { return WF.ScalarRelDerivs(r, u_and_du); }
  NonPHDerivs (RadialWF &wf) : WF(wf)
  { /* Do nothing */ }
};


/////////////////////////////////////////////////////////////////
//                     Inline functions                        //
/////////////////////////////////////////////////////////////////

inline double RadialWF::NormDeriv (double r, double cumulative)
{
  double uval = u(r);
  return (uval*uval);
}

inline Vec2 RadialWF::PseudoDerivs (double r, Vec2 &u_and_du)
{
  Vec2 derivs;  
  derivs[0] = u_and_du[1];
  double A = pot->A(r);
  double B = pot->B(r);
  double V = pot->V(l,r);
  double dAdr = pot->dAdr(r);
  derivs[1] =  1.0/A*
    (-dAdr*u_and_du[1] + (dAdr/r + (double)(l*(l+1))*B/(r*r) 
			  + 2.0*(V-Energy))*u_and_du[0]);
  return derivs;
}


inline Vec2 RadialWF::NormalDerivs(double r, Vec2 &u_and_du)
{
  Vec2 derivs;
  double V = pot->V(l,r);
  derivs[0] = u_and_du[1];
  derivs[1] = ((double)(l*(l+1))/(r*r) + 2.0*(V-Energy))*u_and_du[0];
  return derivs;
}


inline Vec2 RadialWF::ScalarRelDerivs (double r, Vec2 &u_and_du)
{
  const double alpha = 1.0/137.036;
  const double kappa = -1.0;

  Vec2 derivs;
  derivs[0] = u_and_du[1];
  double V = pot->V(l,r);
  double dVdr = pot->dVdr(r);
  double M = 1.0 - alpha*alpha*0.5*(V-Energy);
  
  derivs[1] = ((double)(l*(l+1))/(r*r) + 2.0*M*(V-Energy))*u_and_du[0] 
    - 0.5*alpha*alpha/M*dVdr*(u_and_du[1] + u_and_du[0]*kappa/r);
  
  return derivs;
}

#endif
