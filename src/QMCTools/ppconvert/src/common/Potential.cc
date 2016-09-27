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

#include "Potential.h"

Potential* ReadPotential (IOSectionClass &in)
{
  string type;
  Potential *pot;
  assert (in.ReadVar ("Type", type));
  if (type == "Coulomb")
    pot = new CoulombPot;
  else if (type == "General")
    pot = new GeneralPot;
  else if (type == "Screened")
    pot = new ScreenedPot;
  else if (type == "Spline")
    pot = new SplinePot;
  else if (type == "NLPP")
    pot = new NLPPClass;
  else {
    cerr << "Unrecognize potential type \"" << type << "\".  Exitting.\n";
    exit(1);
  }
  pot->Read(in);
  return pot;
}
