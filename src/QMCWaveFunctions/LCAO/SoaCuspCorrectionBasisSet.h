//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////


/** @file SoaCuspCorrectionBasisSet.h
 *
 * Convert CorrectingOrbitalBasisSet using MultiQuinticSpline1D<T>
 */
#ifndef QMCPLUSPLUS_SOA_CUSPCORRECTION_BASISSET_H
#define QMCPLUSPLUS_SOA_CUSPCORRECTION_BASISSET_H

#include "Configuration.h"
#include "QMCWaveFunctions/BasisSetBase.h"
#include "Particle/DistanceTable.h"
#include "MultiQuinticSpline1D.h"

namespace qmcplusplus
{
// /** Handles a set of correction orbitals per atom
//  *
//  * Reduction over the orbitals - beware of the reduction problem
//  */
template<typename T>
class CuspCorrectionAtomicBasis
{
  using QMCT = QMCTraits;
  typedef MultiQuinticSpline1D<T> RadialSetType;
  typedef ParticleSet::PosType PosType;

  QMCT::RealType r_max_ = 100;
  RadialSetType AOs;
  aligned_vector<size_t> ID;

public:
  CuspCorrectionAtomicBasis(){};

  /** copy constructor */
  CuspCorrectionAtomicBasis(const CuspCorrectionAtomicBasis& a) = default;

  inline void initializeRadialSet(LogGrid<T>& radial_grid, QMCT::IndexType orbital_set_size)
  {
    r_max_ = radial_grid.rmax();
    AOs.initialize(radial_grid, orbital_set_size);
  }

  template<class T1>
  inline void addSpline(int mo_idx, OneDimQuinticSpline<T1>& radial_spline)
  {
    AOs.add_spline(mo_idx, radial_spline);
  }

  inline void evaluate(const T r, T* restrict vals) const
  {
    //assume output vars are zero'd
    if (r >= r_max_)
      return;

    size_t nr = AOs.getNumSplines();
    //FIXME ad-hoc allocation for performance
    std::vector<T> phi(nr);

    AOs.evaluate(r, phi.data());
    for (size_t i = 0; i < nr; ++i)
      //vals[ID[i]]+=phi[i];
      vals[i] += phi[i];
  }

  inline void evaluate_vgl(const T r,
                           const PosType& dr,
                           T* restrict u,
                           T* restrict du_x,
                           T* restrict du_y,
                           T* restrict du_z,
                           T* restrict d2u) const
  {
    //assume output vars are zero'd
    if (r >= r_max_)
      return;

    size_t nr = AOs.getNumSplines();
    //FIXME ad-hoc allocation for performance
    std::vector<T> phi(nr);
    std::vector<T> dphi(nr);
    std::vector<T> d2phi(nr);

    AOs.evaluate(r, phi.data(), dphi.data(), d2phi.data());
    for (size_t i = 0; i < nr; ++i)
    {
      const size_t j = i; //ID[i];
      u[j] += phi[i];
      du_x[j] -= dphi[i] * dr[0] / r; // Displacements have opposite sign (relative to AOS)
      du_y[j] -= dphi[i] * dr[1] / r;
      du_z[j] -= dphi[i] * dr[2] / r;
      d2u[j] += d2phi[i] + 2 * dphi[i] / r;
    }
  }
};
} // namespace qmcplusplus
#endif
