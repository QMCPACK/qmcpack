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
  Potential *BarePot;

  bool IsPH();
  bool NeedsRel();
  double GetCoreRadius();
  double A      (double r);
  double B      (double r);
  double dAdr   (double r);
  double d2Adr2 (double r);

  double V      (int l, double r);
  double dVdr   (int l, double r);
  double d2Vdr2 (int l, double r);

  double V      (double r);
  double dVdr   (double r);
  double d2Vdr2 (double r);
  void   Write(IOSectionClass &out);
  void   Read (IOSectionClass &in);
};



#endif
