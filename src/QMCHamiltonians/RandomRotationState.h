//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_RANDOM_ROTATION_STATE_H
#define QMCPLUSPLUS_RANDOM_ROTATION_STATE_H

#include "Configuration.h"
#include "Utilities/RandomGenerator.h"
#include <array>

namespace qmcplusplus
{

// The angular values for a random rotation matrix
// phi : 0 to 2*pi
// psi : 0 to 2*pi
// cth : -1 to 1  (cth = cos(theta))
class RandomRotationState
{
public:
  using RealType   = QMCTraits::RealType;
  using TensorType = QMCTraits::TensorType;

  RandomRotationState() : rng_vals_({0.0, 0.0, 0.0}) {}
  RandomRotationState(RandomGenerator& rng) : rng_vals_({rng(), rng(), rng()}) {}

  TensorType getRandomRotationMatrix() const
  {
    RealType phi(TWOPI * rng_vals_[0]);
    RealType psi(TWOPI * rng_vals_[1]);
    RealType cth(1.0 - 2 * rng_vals_[2]);
    RealType sph(std::sin(phi)), cph(std::cos(phi));
    RealType sth(std::sqrt(1.0 - cth * cth));
    RealType sps(std::sin(psi)), cps(std::cos(psi));
    // clang-format off
    return TensorType( cph * cth * cps - sph * sps,
                       sph * cth * cps + cph * sps,
                      -sth * cps,
                      -cph * cth * sps - sph * cps,
                      -sph * cth * sps + cph * cps,
                       sth * sps,
                       cph * sth,
                       sph * sth,
                       cth);
    // clang-format on
  }

private:
  std::array<RealType, 3> rng_vals_;
};

} // namespace qmcplusplus
#endif
