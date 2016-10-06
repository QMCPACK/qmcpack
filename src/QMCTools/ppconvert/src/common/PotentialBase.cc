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

#include "PotentialBase.h"

bool 
Potential::IsPH()     
{ 
  return false; 
}

bool
Potential::NeedsRel()
{
  return false;
}

double
Potential::V(int l, double r)
{  return V(r); }

double
Potential::dVdr(int l, double r)
{  return dVdr(r); }

double
Potential::d2Vdr2(int l, double r)
{  return d2Vdr2(r); }
