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

#ifndef SCREENED_POT_H
#define SCREENED_POT_H

#include "PotentialBase.h"
#include "CubicSplineCommon.h"

class ScreenedPot : public Potential
{
public:
  double Charge;
  CubicSplineCommon HXC;
  Potential* BarePot;

  bool IsPH() override;
  bool NeedsRel() override;
  double GetCoreRadius() override;
  double A(double r) override;
  double B(double r) override;
  double dAdr(double r) override;
  double d2Adr2(double r) override;

  double V(int l, double r) override;
  double dVdr(int l, double r) override;
  double d2Vdr2(int l, double r) override;

  double V(double r) override;
  double dVdr(double r) override;
  double d2Vdr2(double r) override;
  void Write(IOSectionClass& out) override;
  void Read(IOSectionClass& in) override;
};


#endif
