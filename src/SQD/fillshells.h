//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef ATOMICHF_FILLSHELLS_H
#define ATOMICHF_FILLSHELLS_H
#include "SQD/HFConfiguration.h"

/**@file fillshells.h
  @brief Fill the closed shells for the nuclear and harmonic potentials
  @authors Jeongnim Kim, Jordan Vincent
  @note  The original Prim was written in F90 by Tim Wilkens.
 */

namespace ohmmshf
{

/**
 *@param mo the
 *@param nmax the number of closed shells
 *@brief Fill the closed shells for the \f$ -Z/r \f$ potential.
 *
 *It is convenient to automatically fill the closed shell system
 for the closest noble-gas configuration, all that is required at
 the input level is to add the valence electrons.  For the post-
 transition elements it is necessary to add the d-states.
 *
 The noble-gas corresponding to the closed shell:
 *
 shell   atom    #electrons
 *
 1       He             2
 *
 2       Ne            10
 *
 3       Ar            18
 *
 4       Kr            32
 *
 5       Xe            54
 *
 6       Rn            86
*/
inline
void FillShellsNucPot(SphericalOrbitalTraits::BasisSetType& mo,
                      int nmax)
{
  int up=1;
  int dn=-1;
  switch(nmax)
  {
  case(1):
    //add 2 orbitals
    mo.add(1,0,0,up,1.0);
    mo.add(1,0,0,dn,1.0);
    mo.Nup += 1;
    mo.Ndown += 1;
    break;
  case(2):
    //add 10 orbitals (Ne)
    mo.add(1,0, 0,up,1.0);
    mo.add(1,0, 0,dn,1.0);
    mo.add(2,0, 0,up,1.0);
    mo.add(2,0, 0,dn,1.0);
    mo.add(2,1,-1,up,1.0);
    mo.add(2,1,-1,dn,1.0);
    mo.add(2,1, 0,up,1.0);
    mo.add(2,1, 0,dn,1.0);
    mo.add(2,1, 1,up,1.0);
    mo.add(2,1, 1,dn,1.0);
    mo.Nup += 5;
    mo.Ndown += 5;
    break;
  case(3):
    //add 18 orbitals (Ar)
    mo.add(1,0, 0,up,1.0);
    mo.add(1,0, 0,dn,1.0);
    mo.add(2,0, 0,up,1.0);
    mo.add(2,0, 0,dn,1.0);
    mo.add(2,1,-1,up,1.0);
    mo.add(2,1,-1,dn,1.0);
    mo.add(2,1, 0,up,1.0);
    mo.add(2,1, 0,dn,1.0);
    mo.add(2,1, 1,up,1.0);
    mo.add(2,1, 1,dn,1.0);
    mo.add(3,0, 0,up,1.0);
    mo.add(3,0, 0,dn,1.0);
    mo.add(3,1,-1,up,1.0);
    mo.add(3,1,-1,dn,1.0);
    mo.add(3,1, 0,up,1.0);
    mo.add(3,1, 0,dn,1.0);
    mo.add(3,1, 1,up,1.0);
    mo.add(3,1, 1,dn,1.0);
    mo.Nup += 9;
    mo.Ndown += 9;
    break;
  case(4):
    //add 36 orbitals (Kr)
    mo.add(1,0, 0,up,1.0);
    mo.add(1,0, 0,dn,1.0);
    mo.add(2,0, 0,up,1.0);
    mo.add(2,0, 0,dn,1.0);
    mo.add(2,1,-1,up,1.0);
    mo.add(2,1,-1,dn,1.0);
    mo.add(2,1, 0,up,1.0);
    mo.add(2,1, 0,dn,1.0);
    mo.add(2,1, 1,up,1.0);
    mo.add(2,1, 1,dn,1.0);
    mo.add(3,0, 0,up,1.0);
    mo.add(3,0, 0,dn,1.0);
    mo.add(3,1,-1,up,1.0);
    mo.add(3,1,-1,dn,1.0);
    mo.add(3,1, 0,up,1.0);
    mo.add(3,1, 0,dn,1.0);
    mo.add(3,1, 1,up,1.0);
    mo.add(3,1, 1,dn,1.0);
    mo.add(3,2,-2,up,1.0);
    mo.add(3,2,-2,dn,1.0);
    mo.add(3,2,-1,up,1.0);
    mo.add(3,2,-1,dn,1.0);
    mo.add(3,2, 0,up,1.0);
    mo.add(3,2, 0,dn,1.0);
    mo.add(3,2, 1,up,1.0);
    mo.add(3,2, 1,dn,1.0);
    mo.add(3,2, 2,up,1.0);
    mo.add(3,2, 2,dn,1.0);
    mo.add(4,0, 0,up,1.0);
    mo.add(4,0, 0,dn,1.0);
    mo.add(4,1,-1,up,1.0);
    mo.add(4,1,-1,dn,1.0);
    mo.add(4,1, 0,up,1.0);
    mo.add(4,1, 0,dn,1.0);
    mo.add(4,1, 1,up,1.0);
    mo.add(4,1, 1,dn,1.0);
    mo.Nup += 18;
    mo.Ndown += 18;
    break;
  case(5):
    //add 54 orbitals (Xe)
    mo.add(1,0, 0,up,1.0);
    mo.add(1,0, 0,dn,1.0);
    mo.add(2,0, 0,up,1.0);
    mo.add(2,0, 0,dn,1.0);
    mo.add(2,1,-1,up,1.0);
    mo.add(2,1,-1,dn,1.0);
    mo.add(2,1, 0,up,1.0);
    mo.add(2,1, 0,dn,1.0);
    mo.add(2,1, 1,up,1.0);
    mo.add(2,1, 1,dn,1.0);
    mo.add(3,0, 0,up,1.0);
    mo.add(3,0, 0,dn,1.0);
    mo.add(3,1,-1,up,1.0);
    mo.add(3,1,-1,dn,1.0);
    mo.add(3,1, 0,up,1.0);
    mo.add(3,1, 0,dn,1.0);
    mo.add(3,1, 1,up,1.0);
    mo.add(3,1, 1,dn,1.0);
    mo.add(3,2,-2,up,1.0);
    mo.add(3,2,-2,dn,1.0);
    mo.add(3,2,-1,up,1.0);
    mo.add(3,2,-1,dn,1.0);
    mo.add(3,2, 0,up,1.0);
    mo.add(3,2, 0,dn,1.0);
    mo.add(3,2, 1,up,1.0);
    mo.add(3,2, 1,dn,1.0);
    mo.add(3,2, 2,up,1.0);
    mo.add(3,2, 2,dn,1.0);
    mo.add(4,0, 0,up,1.0);
    mo.add(4,0, 0,dn,1.0);
    mo.add(4,1,-1,up,1.0);
    mo.add(4,1,-1,dn,1.0);
    mo.add(4,1, 0,up,1.0);
    mo.add(4,1, 0,dn,1.0);
    mo.add(4,1, 1,up,1.0);
    mo.add(4,1, 1,dn,1.0);
    mo.add(4,2,-2,up,1.0);
    mo.add(4,2,-2,dn,1.0);
    mo.add(4,2,-1,up,1.0);
    mo.add(4,2,-1,dn,1.0);
    mo.add(4,2, 0,up,1.0);
    mo.add(4,2, 0,dn,1.0);
    mo.add(4,2, 1,up,1.0);
    mo.add(4,2, 1,dn,1.0);
    mo.add(4,2, 2,up,1.0);
    mo.add(4,2, 2,dn,1.0);
    mo.add(5,0, 0,up,1.0);
    mo.add(5,0, 0,dn,1.0);
    mo.add(5,1,-1,up,1.0);
    mo.add(5,1,-1,dn,1.0);
    mo.add(5,1, 0,up,1.0);
    mo.add(5,1, 0,dn,1.0);
    mo.add(5,1, 1,up,1.0);
    mo.add(5,1, 1,dn,1.0);
    mo.Nup += 27;
    mo.Ndown += 27;
    break;
  case(6):
    //add 86 orbitals (Rn)
    mo.add(1,0, 0,up,1.0);
    mo.add(1,0, 0,dn,1.0);
    mo.add(2,0, 0,up,1.0);
    mo.add(2,0, 0,dn,1.0);
    mo.add(2,1,-1,up,1.0);
    mo.add(2,1,-1,dn,1.0);
    mo.add(2,1, 0,up,1.0);
    mo.add(2,1, 0,dn,1.0);
    mo.add(2,1, 1,up,1.0);
    mo.add(2,1, 1,dn,1.0);
    mo.add(3,0, 0,up,1.0);
    mo.add(3,0, 0,dn,1.0);
    mo.add(3,1,-1,up,1.0);
    mo.add(3,1,-1,dn,1.0);
    mo.add(3,1, 0,up,1.0);
    mo.add(3,1, 0,dn,1.0);
    mo.add(3,1, 1,up,1.0);
    mo.add(3,1, 1,dn,1.0);
    mo.add(3,2,-2,up,1.0);
    mo.add(3,2,-2,dn,1.0);
    mo.add(3,2,-1,up,1.0);
    mo.add(3,2,-1,dn,1.0);
    mo.add(3,2, 0,up,1.0);
    mo.add(3,2, 0,dn,1.0);
    mo.add(3,2, 1,up,1.0);
    mo.add(3,2, 1,dn,1.0);
    mo.add(3,2, 2,up,1.0);
    mo.add(3,2, 2,dn,1.0);
    mo.add(4,0, 0,up,1.0);
    mo.add(4,0, 0,dn,1.0);
    mo.add(4,1,-1,up,1.0);
    mo.add(4,1,-1,dn,1.0);
    mo.add(4,1, 0,up,1.0);
    mo.add(4,1, 0,dn,1.0);
    mo.add(4,1, 1,up,1.0);
    mo.add(4,1, 1,dn,1.0);
    mo.add(4,2,-2,up,1.0);
    mo.add(4,2,-2,dn,1.0);
    mo.add(4,2,-1,up,1.0);
    mo.add(4,2,-1,dn,1.0);
    mo.add(4,2, 0,up,1.0);
    mo.add(4,2, 0,dn,1.0);
    mo.add(4,2, 1,up,1.0);
    mo.add(4,2, 1,dn,1.0);
    mo.add(4,2, 2,up,1.0);
    mo.add(4,2, 2,dn,1.0);
    mo.add(5,0, 0,up,1.0);
    mo.add(5,0, 0,dn,1.0);
    mo.add(5,1,-1,up,1.0);
    mo.add(5,1,-1,dn,1.0);
    mo.add(5,1, 0,up,1.0);
    mo.add(5,1, 0,dn,1.0);
    mo.add(5,1, 1,up,1.0);
    mo.add(5,1, 1,dn,1.0);
    mo.add(4,3,-3,up,1.0);
    mo.add(4,3,-3,dn,1.0);
    mo.add(4,3,-2,up,1.0);
    mo.add(4,3,-2,dn,1.0);
    mo.add(4,3,-1,up,1.0);
    mo.add(4,3,-1,dn,1.0);
    mo.add(4,3, 0,up,1.0);
    mo.add(4,3, 0,dn,1.0);
    mo.add(4,3, 1,up,1.0);
    mo.add(4,3, 1,dn,1.0);
    mo.add(4,3, 2,up,1.0);
    mo.add(4,3, 2,dn,1.0);
    mo.add(4,3, 3,up,1.0);
    mo.add(4,3, 3,dn,1.0);
    mo.add(5,2,-2,up,1.0);
    mo.add(5,2,-2,dn,1.0);
    mo.add(5,2,-1,up,1.0);
    mo.add(5,2,-1,dn,1.0);
    mo.add(5,2, 0,up,1.0);
    mo.add(5,2, 0,dn,1.0);
    mo.add(5,2, 1,up,1.0);
    mo.add(5,2, 1,dn,1.0);
    mo.add(5,2, 2,up,1.0);
    mo.add(5,2, 2,dn,1.0);
    mo.add(6,0, 0,up,1.0);
    mo.add(6,0, 0,dn,1.0);
    mo.add(6,1,-1,up,1.0);
    mo.add(6,1,-1,dn,1.0);
    mo.add(6,1, 0,up,1.0);
    mo.add(6,1, 0,dn,1.0);
    mo.add(6,1, 1,up,1.0);
    mo.add(6,1, 1,dn,1.0);
    mo.Nup += 43;
    mo.Ndown += 43;
    break;
  }
}

inline
void FillShellsHarmPot(SphericalOrbitalTraits::BasisSetType& mo,
                       int nmax, int np_offset=0)
{
  int up=1;
  int dn=-1;
  switch(nmax)
  {
  case(1):
    //add 2 orbitals
    mo.add(0,0,0,up,1.0);
    mo.add(0,0,0,dn,1.0);
    mo.Nup += 1;
    mo.Ndown += 1;
    break;
  case(2):
    //add 8 orbitals
    mo.add(np_offset,0,0,up,1.0);
    mo.add(np_offset,0,0,dn,1.0);
    mo.add(np_offset,1,-1,up,1.0);
    mo.add(np_offset,1,-1,dn,1.0);
    mo.add(np_offset,1, 0,up,1.0);
    mo.add(np_offset,1, 0,dn,1.0);
    mo.add(np_offset,1, 1,up,1.0);
    mo.add(np_offset,1, 1,dn,1.0);
    mo.Nup += 4;
    mo.Ndown += 4;
    break;
  case(3):
    //add 18 orbitals
    mo.add(np_offset,0,0,up,1.0);
    mo.add(np_offset,0,0,dn,1.0);
    mo.add(np_offset,1,-1,up,1.0);
    mo.add(np_offset,1,-1,dn,1.0);
    mo.add(np_offset,1, 0,up,1.0);
    mo.add(np_offset,1, 0,dn,1.0);
    mo.add(np_offset,1, 1,up,1.0);
    mo.add(np_offset,1, 1,dn,1.0);
    mo.add(np_offset,2,-2,up,1.0);
    mo.add(np_offset,2,-2,dn,1.0);
    mo.add(np_offset,2,-1,up,1.0);
    mo.add(np_offset,2,-1,dn,1.0);
    mo.add(np_offset,2, 0,up,1.0);
    mo.add(np_offset,2, 0,dn,1.0);
    mo.add(np_offset,2, 1,up,1.0);
    mo.add(np_offset,2, 1,dn,1.0);
    mo.add(np_offset,2, 2,up,1.0);
    mo.add(np_offset,2, 2,dn,1.0);
    mo.Nup += 9;
    mo.Ndown+= 9;
    break;
  case(4):
    //add 18 orbitals
    mo.add(np_offset,0,0,up,1.0);
    mo.add(np_offset,0,0,dn,1.0);
    mo.add(np_offset,1,-1,up,1.0);
    mo.add(np_offset,1,-1,dn,1.0);
    mo.add(np_offset,1, 0,up,1.0);
    mo.add(np_offset,1, 0,dn,1.0);
    mo.add(np_offset,1, 1,up,1.0);
    mo.add(np_offset,1, 1,dn,1.0);
    mo.add(np_offset,2,-2,up,1.0);
    mo.add(np_offset,2,-2,dn,1.0);
    mo.add(np_offset,2,-1,up,1.0);
    mo.add(np_offset,2,-1,dn,1.0);
    mo.add(np_offset,2, 0,up,1.0);
    mo.add(np_offset,2, 0,dn,1.0);
    mo.add(np_offset,2, 1,up,1.0);
    mo.add(np_offset,2, 1,dn,1.0);
    mo.add(np_offset,2, 2,up,1.0);
    mo.add(np_offset,2, 2,dn,1.0);
    mo.add(np_offset+1,0,0,up,1.0);
    mo.add(np_offset+1,0,0,dn,1.0);
    mo.Nup += 10;
    mo.Ndown+= 10;
    break;
  case(5):
    //add 34 orbitals
    mo.add(np_offset,0,0,up,1.0);
    mo.add(np_offset,0,0,dn,1.0);
    mo.add(np_offset,1,-1,up,1.0);
    mo.add(np_offset,1,-1,dn,1.0);
    mo.add(np_offset,1, 0,up,1.0);
    mo.add(np_offset,1, 0,dn,1.0);
    mo.add(np_offset,1, 1,up,1.0);
    mo.add(np_offset,1, 1,dn,1.0);
    mo.add(np_offset,2,-2,up,1.0);
    mo.add(np_offset,2,-2,dn,1.0);
    mo.add(np_offset,2,-1,up,1.0);
    mo.add(np_offset,2,-1,dn,1.0);
    mo.add(np_offset,2, 0,up,1.0);
    mo.add(np_offset,2, 0,dn,1.0);
    mo.add(np_offset,2, 1,up,1.0);
    mo.add(np_offset,2, 1,dn,1.0);
    mo.add(np_offset,2, 2,up,1.0);
    mo.add(np_offset,2, 2,dn,1.0);
    mo.add(np_offset+1,0,0,up,1.0);
    mo.add(np_offset+1,0,0,dn,1.0);
    mo.add(np_offset,3,-3,up,1.0);
    mo.add(np_offset,3,-3,dn,1.0);
    mo.add(np_offset,3,-2,up,1.0);
    mo.add(np_offset,3,-2,dn,1.0);
    mo.add(np_offset,3,-1,up,1.0);
    mo.add(np_offset,3,-1,dn,1.0);
    mo.add(np_offset,3, 0,up,1.0);
    mo.add(np_offset,3, 0,dn,1.0);
    mo.add(np_offset,3, 1,up,1.0);
    mo.add(np_offset,3, 1,dn,1.0);
    mo.add(np_offset,3, 2,up,1.0);
    mo.add(np_offset,3, 2,dn,1.0);
    mo.add(np_offset,3, 3,up,1.0);
    mo.add(np_offset,3, 3,dn,1.0);
    mo.add(np_offset+1,1,-1,up,1.0);
    mo.add(np_offset+1,1,-1,dn,1.0);
    mo.add(np_offset+1,1, 0,up,1.0);
    mo.add(np_offset+1,1, 0,dn,1.0);
    mo.add(np_offset+1,1, 1,up,1.0);
    mo.add(np_offset+1,1, 1,dn,1.0);
    mo.Nup += 17;
    mo.Ndown += 17;
    break;
  case(6):
    //not doing this for now
    //add 70 orbitals
    mo.add(0,0,0,up,1.0);
    mo.add(0,0,0,dn,1.0);
    mo.add(0,1,-1,up,1.0);
    mo.add(0,1,-1,dn,1.0);
    mo.add(0,1, 0,up,1.0);
    mo.add(0,1, 0,dn,1.0);
    mo.add(0,1, 1,up,1.0);
    mo.add(0,1, 1,dn,1.0);
    mo.add(0,2,-2,up,1.0);
    mo.add(0,2,-2,dn,1.0);
    mo.add(0,2,-1,up,1.0);
    mo.add(0,2,-1,dn,1.0);
    mo.add(0,2, 0,up,1.0);
    mo.add(0,2, 0,dn,1.0);
    mo.add(0,2, 1,up,1.0);
    mo.add(0,2, 1,dn,1.0);
    mo.add(0,2, 2,up,1.0);
    mo.add(0,2, 2,dn,1.0);
    mo.add(1,0, 0,up,1.0);
    mo.add(1,0, 0,dn,1.0);
    mo.add(0,3,-3,up,1.0);
    mo.add(0,3,-3,dn,1.0);
    mo.add(0,3,-2,up,1.0);
    mo.add(0,3,-2,dn,1.0);
    mo.add(0,3,-1,up,1.0);
    mo.add(0,3,-1,dn,1.0);
    mo.add(0,3, 0,up,1.0);
    mo.add(0,3, 0,dn,1.0);
    mo.add(0,3, 1,up,1.0);
    mo.add(0,3, 1,dn,1.0);
    mo.add(0,3, 2,up,1.0);
    mo.add(0,3, 2,dn,1.0);
    mo.add(0,3, 3,up,1.0);
    mo.add(0,3, 3,dn,1.0);
    mo.add(1,1,-1,up,1.0);
    mo.add(1,1,-1,dn,1.0);
    mo.add(1,1, 0,up,1.0);
    mo.add(1,1, 0,dn,1.0);
    mo.add(1,1, 1,up,1.0);
    mo.add(1,1, 1,dn,1.0);
    mo.add(0,4,-4,up,1.0);
    mo.add(0,4,-4,dn,1.0);
    mo.add(0,4,-3,up,1.0);
    mo.add(0,4,-3,dn,1.0);
    mo.add(0,4,-2,up,1.0);
    mo.add(0,4,-2,dn,1.0);
    mo.add(0,4,-1,up,1.0);
    mo.add(0,4,-1,dn,1.0);
    mo.add(0,4, 0,up,1.0);
    mo.add(0,4, 0,dn,1.0);
    mo.add(0,4, 1,up,1.0);
    mo.add(0,4, 1,dn,1.0);
    mo.add(0,4, 2,up,1.0);
    mo.add(0,4, 2,dn,1.0);
    mo.add(0,4, 3,up,1.0);
    mo.add(0,4, 3,dn,1.0);
    mo.add(0,4, 4,up,1.0);
    mo.add(0,4, 4,dn,1.0);
    mo.add(2,0, 0,up,1.0);
    mo.add(2,0, 0,dn,1.0);
    mo.add(1,2,-2,up,1.0);
    mo.add(1,2,-2,dn,1.0);
    mo.add(1,2,-1,up,1.0);
    mo.add(1,2,-1,dn,1.0);
    mo.add(1,2, 0,up,1.0);
    mo.add(1,2, 0,dn,1.0);
    mo.add(1,2, 1,up,1.0);
    mo.add(1,2, 1,dn,1.0);
    mo.add(1,2, 2,up,1.0);
    mo.add(1,2, 2,dn,1.0);
    mo.Nup += 35;
    mo.Ndown += 35;
    break;
  case(7):
    //add 112 orbitals
    mo.add(0,0,0,up,1.0);
    mo.add(0,0,0,dn,1.0);
    mo.add(0,1,-1,up,1.0);
    mo.add(0,1,-1,dn,1.0);
    mo.add(0,1, 0,up,1.0);
    mo.add(0,1, 0,dn,1.0);
    mo.add(0,1, 1,up,1.0);
    mo.add(0,1, 1,dn,1.0);
    mo.add(0,2,-2,up,1.0);
    mo.add(0,2,-2,dn,1.0);
    mo.add(0,2,-1,up,1.0);
    mo.add(0,2,-1,dn,1.0);
    mo.add(0,2, 0,up,1.0);
    mo.add(0,2, 0,dn,1.0);
    mo.add(0,2, 1,up,1.0);
    mo.add(0,2, 1,dn,1.0);
    mo.add(0,2, 2,up,1.0);
    mo.add(0,2, 2,dn,1.0);
    mo.add(1,0, 0,up,1.0);
    mo.add(1,0, 0,dn,1.0);
    mo.add(0,3,-3,up,1.0);
    mo.add(0,3,-3,dn,1.0);
    mo.add(0,3,-2,up,1.0);
    mo.add(0,3,-2,dn,1.0);
    mo.add(0,3,-1,up,1.0);
    mo.add(0,3,-1,dn,1.0);
    mo.add(0,3, 0,up,1.0);
    mo.add(0,3, 0,dn,1.0);
    mo.add(0,3, 1,up,1.0);
    mo.add(0,3, 1,dn,1.0);
    mo.add(0,3, 2,up,1.0);
    mo.add(0,3, 2,dn,1.0);
    mo.add(0,3, 3,up,1.0);
    mo.add(0,3, 3,dn,1.0);
    mo.add(1,1,-1,up,1.0);
    mo.add(1,1,-1,dn,1.0);
    mo.add(1,1, 0,up,1.0);
    mo.add(1,1, 0,dn,1.0);
    mo.add(1,1, 1,up,1.0);
    mo.add(1,1, 1,dn,1.0);
    mo.add(0,4,-4,up,1.0);
    mo.add(0,4,-4,dn,1.0);
    mo.add(0,4,-3,up,1.0);
    mo.add(0,4,-3,dn,1.0);
    mo.add(0,4,-2,up,1.0);
    mo.add(0,4,-2,dn,1.0);
    mo.add(0,4,-1,up,1.0);
    mo.add(0,4,-1,dn,1.0);
    mo.add(0,4, 0,up,1.0);
    mo.add(0,4, 0,dn,1.0);
    mo.add(0,4, 1,up,1.0);
    mo.add(0,4, 1,dn,1.0);
    mo.add(0,4, 2,up,1.0);
    mo.add(0,4, 2,dn,1.0);
    mo.add(0,4, 3,up,1.0);
    mo.add(0,4, 3,dn,1.0);
    mo.add(0,4, 4,up,1.0);
    mo.add(0,4, 4,dn,1.0);
    mo.add(2,0, 0,up,1.0);
    mo.add(2,0, 0,dn,1.0);
    mo.add(1,2,-2,up,1.0);
    mo.add(1,2,-2,dn,1.0);
    mo.add(1,2,-1,up,1.0);
    mo.add(1,2,-1,dn,1.0);
    mo.add(1,2, 0,up,1.0);
    mo.add(1,2, 0,dn,1.0);
    mo.add(1,2, 1,up,1.0);
    mo.add(1,2, 1,dn,1.0);
    mo.add(1,2, 2,up,1.0);
    mo.add(1,2, 2,dn,1.0);
    mo.add(0,5,-5,up,1.0);
    mo.add(0,5,-5,dn,1.0);
    mo.add(0,5,-4,up,1.0);
    mo.add(0,5,-4,dn,1.0);
    mo.add(0,5,-3,up,1.0);
    mo.add(0,5,-3,dn,1.0);
    mo.add(0,5,-2,up,1.0);
    mo.add(0,5,-2,dn,1.0);
    mo.add(0,5,-1,up,1.0);
    mo.add(0,5,-1,dn,1.0);
    mo.add(0,5, 0,up,1.0);
    mo.add(0,5, 0,dn,1.0);
    mo.add(0,5, 1,up,1.0);
    mo.add(0,5, 1,dn,1.0);
    mo.add(0,5, 2,up,1.0);
    mo.add(0,5, 2,dn,1.0);
    mo.add(0,5, 3,up,1.0);
    mo.add(0,5, 3,dn,1.0);
    mo.add(0,5, 4,up,1.0);
    mo.add(0,5, 4,dn,1.0);
    mo.add(0,5, 5,up,1.0);
    mo.add(0,5, 5,dn,1.0);
    mo.add(1,3,-3,up,1.0);
    mo.add(1,3,-3,dn,1.0);
    mo.add(1,3,-2,up,1.0);
    mo.add(1,3,-2,dn,1.0);
    mo.add(1,3,-1,up,1.0);
    mo.add(1,3,-1,dn,1.0);
    mo.add(1,3, 0,up,1.0);
    mo.add(1,3, 0,dn,1.0);
    mo.add(1,3, 1,up,1.0);
    mo.add(1,3, 1,dn,1.0);
    mo.add(1,3, 2,up,1.0);
    mo.add(1,3, 2,dn,1.0);
    mo.add(1,3, 3,up,1.0);
    mo.add(1,3, 3,dn,1.0);
    mo.add(2,1,-1,up,1.0);
    mo.add(2,1,-1,dn,1.0);
    mo.add(2,1, 0,up,1.0);
    mo.add(2,1, 0,dn,1.0);
    mo.add(2,1, 1,up,1.0);
    mo.add(2,1, 1,dn,1.0);
    mo.Nup += 56;
    mo.Ndown += 56;
    break;
  }
}

}
#endif
