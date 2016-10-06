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

#ifndef NLPP_CLASS_H
#define NLPP_CLASS_H

#include "PotentialBase.h"
#include "IO.h"
#include "CubicSplineCommon.h"
#include <vector>
using namespace IO;
using namespace blitz;

class ChannelPotential
{
protected:
  // The projector is given by norm*deltaVl(r)*ul(r)
  double ProjectorNorm;
  inline double jl(int l, double x);
  LinearGrid qGrid;
  double qCurrent, rCurrent;
  typedef enum { NORM, EKB, ZETA_Q, CHI_R, CHECK_CHI_R} IntegrandType;
  IntegrandType Job;
  double A(double q, double qp);
public:
  int l;
  // V stores the potential
  // DeltaV stores the potential less to the local potential
  int n_principal;
  CubicSplineCommon V, DeltaV, u;
  double rc, R0;
  double Occupation, Eigenvalue;
  // The Kleinman-Bylander projection energy
  double E_KB;
  // The unfiltered projection operators in real-space and reciprocal
  // space.  
  CubicSplineCommon zeta_r, zeta_q;
  // These store the real and reciprocal-space representation of the
  // filtered projector
  CubicSplineCommon chi_r, chi_q;

  // This is used as the integrand for the various integrals that are required.
  inline double operator()(double x);

  void SetupProjector(double G_max, double G_FFT);
  void Read  (IOSectionClass &in, Grid* grid);
  void Write (IOSectionClass &out);
};

class NLPPClass : public Potential
{
protected:
  int lLocal;
  std::vector<ChannelPotential> Vl;  
  // The charge of the ion.  The potential should have a tail of
  // V(r) = -Zion/r for large r.
  double Zion;
  int AtomicNumber;
  std::string Symbol;
  Grid *PotentialGrid;
public:
  // General accessor functions
  bool IsNonlocal();
  inline int LocalChannel()                { return lLocal;           } 
  inline int NumChannels()                 { return Vl.size();        }
  inline CubicSplineCommon& GetLocalSpline ()    { return Vl[lLocal].V;     }
  inline double GetValenceCharge()         { return Zion;             }

  // Nonlocal part accessor functions:
  inline double GetChi_r (int l, double r)  { return Vl[l].chi_r(r);  }
  inline double GetZeta_r (int l, double r) { return Vl[l].zeta_r(r); }
  inline double GetChi_q (int l, double q)  { return Vl[l].chi_q(q);  }
  inline double GetZeta_q (int l, double q) { return Vl[l].zeta_q(q); }
  // HACK HACK HACK
  //inline double GetChi_r (int l, double r) { return Vl[l].zeta_r(r); }
  inline double GetE_KB (int l)            { return Vl[l].E_KB;       }
  inline double GetR0 (int l)              { return Vl[l].R0;         }
  inline double Getrc(int l)               { return Vl[l].rc;         }
  inline double GetDeltaV(int l, double r) { return Vl[l].DeltaV(r);  }
  // Override default for local potentials
  double V     (int l, double r);
  double dVdr  (int l, double r);
  double d2Vdr2(int l, double r);


  // Required member functions:  These give information about the
  // local part of the pseudopotential only
  double V(double r);
  double dVdr(double r);
  double d2Vdr2(double r);

  // IO routines
  void Write(IOSectionClass &out);
  void Read(IOSectionClass &in);
  void SetupProjectors(double G_max, double G_FFT);
};


inline double
ChannelPotential::jl(int l, double x)
{
  if (std::abs(x) > 0.0) {
    if (l == 0)
      return sin(x)/x;
    else if (l == 1) {
      if (x < 1.0e-4)
	return x/3.0 - x*x*x/30.0 + x*x*x*x*x/840.0 - x*x*x*x*x*x*x/45360.0;
      else
	return sin(x)/(x*x) - cos(x)/x;
    }
    else if (l == 2) {
      if (x < 1.0e-2)
	return x*x/15.0 - x*x*x*x/210.0 + x*x*x*x*x*x/7560.0 - x*x*x*x*x*x*x*x/498960.0;
      else
	return ((3.0/(x*x*x) - 1.0/x)*sin(x) - 3.0/(x*x)*cos(x));
    }
    else {
      std::cerr << "j(l,x) not implemented for l > 2.\n";
      abort();
    }
  }
  else { // x -> 0 limit
    if (l == 0)
      return 1.0;
    else
      return 0.0;
  }
}


inline double
ChannelPotential::operator()(double x)
{
  switch (Job) {
  case NORM:
    return u(x)*DeltaV(x)*DeltaV(x)*u(x);
  case EKB:
    return u(x)*DeltaV(x)*u(x);
  case ZETA_Q:
    return jl(l,qCurrent*x)*x*x*zeta_r(x);
  case CHI_R:
    return 2.0/M_PI * x*x*chi_q(x)*jl(l,x*rCurrent);
  case CHECK_CHI_R:
    return chi_r(x)*chi_r(x)*x*x;
  default:
    return 0.0;
  }
}



#endif
