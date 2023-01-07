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
inline QMCTraits::TensorType generateRandomRotationMatrix(RandomGenerator& rng)
{
  return generateRotationMatrix(rng(), rng(), rng());
}

} // namespace qmcplusplus
#endif
