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


#ifndef QMCPLUSPLUS_RANDOM_ROTATION_MATRIX_H
#define QMCPLUSPLUS_RANDOM_ROTATION_MATRIX_H

#include "Configuration.h"
#include "Utilities/RandomGenerator.h"
#include "Numerics/RotationMatrix3D.h"
#include <array>

namespace qmcplusplus
{

/// Create a random 3D rotation matrix using a random generator
inline QMCTraits::TensorType generateRandomRotationMatrix(RandomBase<QMCTraits::FullPrecRealType>& rng)
{
  auto rng1 = rng();
  auto rng2 = rng();
  auto rng3 = rng();
  return generateRotationMatrix<QMCTraits::RealType>(rng1, rng2, rng3);
  // The order of evaluation of function arguments is unspecified by the standard.
  // The following code will cause failures in the deterministic tests.
  // return generateRotationMatrix<QMCTraits::RealType>(rng(), rng(), rng());
}

} // namespace qmcplusplus
#endif
