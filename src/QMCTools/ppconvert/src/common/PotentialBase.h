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

#ifndef POTENTIAL_BASE_H
#define POTENTIAL_BASE_H

#include "IO.h"
using namespace IO;

class Potential
{
public:
  // Optional member functions -- if you're not a pseudoHamiltonian,
  // you do not need to define these
  virtual bool IsPH();
  virtual bool IsNonlocal()        { return false; }
  // Nonlocal version of functions
  virtual double V     (int l, double r);
  virtual double dVdr  (int l, double r);
  virtual double d2Vdr2(int l, double r);
  virtual bool NeedsRel(); 
  virtual double GetCoreRadius()   { return 0.0; }
  virtual double A      (double r) { return 1.0; }
  virtual double B      (double r) { return 1.0; }
  virtual double dAdr   (double r) { return 0.0; }
  virtual double d2Adr2 (double r) { return 0.0; }

  // Required member functions
  virtual double V(double r) = 0;
  virtual double dVdr(double r) = 0;
  virtual double d2Vdr2(double r) = 0;
  virtual void Write(IOSectionClass &out) = 0;
  virtual void Read(IOSectionClass &in) = 0;
  virtual double X_k (double rcut, double k) { return 0.0; }
};

Potential* ReadPotential (IOSectionClass &in);

#endif
