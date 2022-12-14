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

template<typename T>
inline Tensor<T, 3> getRotationMatrix(T rng1, T rng2, T rng3)
{
  using Tensor = QMCTraits::TensorType;
  // The angular values for a random rotation matrix
  // phi : 0 to 2*pi
  // psi : 0 to 2*pi
  // cth : -1 to 1  (cth = cos(theta))
  T phi(TWOPI * rng1);
  T psi(TWOPI * rng2);
  T cth(1.0 - 2 * rng3);
  T sph(std::sin(phi)), cph(std::cos(phi));
  T sth(std::sqrt(1.0 - cth * cth));
  T sps(std::sin(psi)), cps(std::cos(psi));
  // clang-format off
  return Tensor(cph * cth * cps - sph * sps,
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

inline QMCTraits::TensorType getRandomRotationMatrix(RandomGenerator &rng)
{
  return getRotationMatrix(rng(), rng(), rng());
}

} // namespace qmcplusplus
#endif
