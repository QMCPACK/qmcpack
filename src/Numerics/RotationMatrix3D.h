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


#ifndef QMCPLUSPLUS_ROTATION_MATRIX_3D_H
#define QMCPLUSPLUS_ROTATION_MATRIX_3D_H

#include "Containers/OhmmsPETE/Tensor.h"
#include "config/stdlib/Constants.h" // For TWOPI

namespace qmcplusplus
{

/// Create a random 3D rotation matrix from three random numbers.
/// Each input random number should be in the range [0,1).
/// See the script Numerics/tests/derive_random_rotation.py for more information about the algorithm.

template<typename T>
inline Tensor<T, 3> generateRotationMatrix(T rng1, T rng2, T rng3)
{
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
  return Tensor<T,3>( cph * cth * cps - sph * sps,
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

} // namespace qmcplusplus
#endif
