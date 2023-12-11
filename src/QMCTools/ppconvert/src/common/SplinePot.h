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

#ifndef SPLINE_POT_H
#define SPLINE_POT_H

#include "PotentialBase.h"
#include "CubicSplineCommon.h"

/// This class stores a tabulated potential and interpolates the data
/// with a cubic spline.  In case r is outside the tabulated grid, it
/// optionally calls Vouter.
class SplinePot : public Potential
{
protected:
public:
  /// This stores
  CubicSplineCommon Spline;
  /// This is an optionally set potential that kicks in outside the
  /// maximum value of the grid point.
  Potential* Vouter;

  double V(double r) override;
  double dVdr(double r) override;
  double d2Vdr2(double r) override;
  void Write(IOSectionClass& out) override;
  void Read(IOSectionClass& in) override;
  SplinePot() : Vouter(NULL)
  { /* No nothing else for now */
  }
};

#endif
