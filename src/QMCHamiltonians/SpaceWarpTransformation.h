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
 *\brief This implements the differential space warp transformation for ZVZB estimators given by Sorella & Capriotti
 *       J. Chem. Phys., 133, 23411 (2010), https://doi.org/10.1063/1.3516208 
 *
 */
struct SpaceWarpTransformation : public QMCTraits
{
  typedef ParticleSet::ParticleScalar_t ParticleScalar_t;
  typedef ParticleSet::ParticlePos_t Force_t;
  typedef ParticleSet::ParticleGradient_t ParticleGradient_t;

  SpaceWarpTransformation(ParticleSet& elns, const ParticleSet& ions);
  ~SpaceWarpTransformation(){};


  RealType f(RealType r);  //space warp transformation function F.
  RealType df(RealType r); //The gradient of F(r)

  void setPow(RealType swpow_in)
  {
    swpow = swpow_in;
  }; //This sets the exponent for the power law space warp transformation.

  void computeSWTIntermediates(ParticleSet& P,
                               const ParticleSet& ions); //This computes the intermediate matrices required to build all
                                                         //space warp components and gradients.

  void getSWT(int iat,
              ParticleScalar_t& w,
              Force_t& grad_w); //For each ion component iat, this generates all the required space warp quantities to
                                //generate the "space warp" contribution to the iat-th force component.

  void computeSWT(ParticleSet& elec,
                  const ParticleSet& ions,
                  Force_t& dEl,
                  ParticleGradient_t& dlogpsi,
                  Force_t& el_contribution,
                  Force_t& psi_contribution); //Takes in precomputed grad(E_L) and grad(logPsi), and uses
                                              //this to compute the ZV (el_contribution) and ZB (psi_contribution)
                                              //space warp contributions.

  const int myTableIndex;
  const int Nelec;
  const int Nions;

  RealType swpow;           //Power of space warp transformation.  Right now, r^{-swpow}.
  Matrix<RealType> warpval; //Nelec x Nion matrix of F(|r_i-R_J|)
  Matrix<PosType> gradval;  //Nelec x Nion matrix of \nabla_i F(|r_i-R_J|)
};
} // namespace qmcplusplus
#endif
