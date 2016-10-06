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

#ifndef GENERAL_POT_H
#define GENERAL_POT_H

#include "PotentialBase.h"
#include "CubicSplineCommon.h"

class GeneralPot : public Potential
{
protected:
  Grid *PotGrid;
  CubicSplineCommon PotSpline;
  double Z;
public:
  double V      (double r);
  double dVdr   (double r);
  double d2Vdr2 (double r);

  void Write (IOSectionClass &out);
  void Read (IOSectionClass &in);
  GeneralPot();
  ~GeneralPot();
};

#endif
