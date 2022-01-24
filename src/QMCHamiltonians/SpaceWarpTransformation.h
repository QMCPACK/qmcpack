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
class SpaceWarpTransformation : public QMCTraits
{
  using ParticleScalar   = ParticleSet::ParticleScalar;
  using Force_t          = ParticleSet::ParticlePos;
  using ParticleGradient = ParticleSet::ParticleGradient;

public:
  SpaceWarpTransformation(ParticleSet& elns, const ParticleSet& ions);

  /** Space warp transformation function F(r).
   *
   * \param[in] r, the distance
   * \param[out] value of F(r)
   */
  RealType f(RealType r);

  /** Derivative of space warp transformation function F(r) w.r.t. r.
   *
   * \param[in] r, the distance
   * \param[out] value of F'(r)
   */
  RealType df(RealType r);

  /** Sets the exponent for power law space warp transformation
   *
   * \param[in] swpow_in the exponent
   */
  inline void setPow(RealType swpow_in) { swpow = swpow_in; };

  /** Generates required space warp quantities to generate the actual "Space Warp" contribution to the
   *  iat-th force component.
   *  \param[in] iat the ion index for the force.  
   *  \param[out] w, w_iat(r_i) for each i, where i is the electron index. 
   *  \param[out] grad_w,  grad_i w_iat(r_i) for each i, where i is the electron index.
   */
  void getSWT(int iat, ParticleScalar& w, Force_t& grad_w);

  /** Takes in precomputed grad(E_L) and grad(logPsi) and computes the ZV and ZB space warp contributions 
   *  to the force.
   *
   *  \param[in] elec, electron particle set. 
   *  \param[in] ions, ion particle set.
   *  \param[in] dEl, grad_i(E_L) for each electron i. E_L is the local energy.
   *  \param[in] dlogpsi, grad_i(logPsi) for each electron i.  
   *  \param[out] el_contribution, The zero-variance contribution from space warp.
   *  \param[out] psi_contribution, the zero-bias contribution from space warp.  Modifies the grad_I(logPsi) terms.
   */
  void computeSWT(ParticleSet& elec,
                  const ParticleSet& ions,
                  Force_t& dEl,
                  ParticleGradient& dlogpsi,
                  Force_t& el_contribution,
                  Force_t& psi_contribution);

private:
  /** Computes intermediate matrices required to build all space warp components and gradients.
   *  The intermediates calculated are "warpval" and "gradval".
   *
   * \param[in] P, the electron particle set.  
   * \param[in] ions, the ion particle set.  
   */
  void computeSWTIntermediates(ParticleSet& P, const ParticleSet& ions);


  ///The electron-ion table index in electron table.
  const int myTableIndex;
  const int Nelec;
  const int Nions;

  /// Power of space warp transformation.  Right now, r^{-swpow}.
  RealType swpow;
  /// Nelec x Nion matrix of F(|r_i-R_J|)
  Matrix<RealType> warpval;
  /// Nelec x Nion matrix of \nabla_i F(|r_i-R_J|)
  Matrix<PosType> gradval;
};
} // namespace qmcplusplus
#endif
