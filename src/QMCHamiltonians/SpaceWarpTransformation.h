//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Raymond Clay
//
// File created by: Raymond Clay, Sandia National Laboratory, rclay@sandia.gov
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SPACEWARP_H
#define QMCPLUSPLUS_SPACEWARP_H

#include "Configuration.h"
#include "Particle/ParticleSet.h"

namespace qmcplusplus
{
/** @ingroup hamiltonian
 *\brief Calculates the AA Coulomb potential using PBCs
 *
 * Functionally identical to CoulombPBCAB but uses a templated version of
 * LRHandler.
 */
struct SpaceWarpTransformation : public QMCTraits
{
  typedef ParticleSet::ParticleScalar_t ParticleScalar_t;
  typedef ParticleSet::ParticlePos_t Force_t;
  typedef ParticleSet::ParticleGradient_t ParticleGradient_t;

  SpaceWarpTransformation(ParticleSet& elns, ParticleSet& ions);
  ~SpaceWarpTransformation(){};

  RealType f(RealType r, RealType a); 
  RealType df(RealType r, RealType a); 

  void setPow(RealType swpow_in){swpow=swpow_in;};
  void computeSWTIntermediates(ParticleSet& P, ParticleSet& ions); //This computes the intermediate matrices required to build all 
                                                      //space warp components and gradients.

  void getSWT(int iat, ParticleScalar_t& w, Force_t& grad_w); //For each ion component iat, this generates all the required space warp quantities to 
                                                              //generate the "space warp" contribution to the iat-th force component.
                                                              
  void computeSWT(ParticleSet& elec, ParticleSet& ions, Force_t& dEl, ParticleGradient_t& dlogpsi, Force_t& el_contribution, Force_t& psi_contribution); //Takes in precomputed grad(E_L) and grad(logPsi), and uses
                                                                                                                        //this to compute the ZV (el_contribution) and ZB (psi_contribution) 
                                                                                                                        //space warp contributions.                                                               

  int myTableIndex;
  int Nelec;
  int Nions;

  RealType swpow; //Power of space warp transformation.  Right now, r^{-swpow}.
  Matrix<RealType> warpval;  //Nelec x Nion matrix of F(|r_i-R_J|)
  Matrix<PosType>  gradval;  //Nelec x Nion matrix of \nabla_i F(|r_i-R_J|)
};
}
#endif


