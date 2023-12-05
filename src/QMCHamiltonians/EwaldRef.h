//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

/**@file EwaldRef.h
 *
 * @brief Computes Ewald sums of the potential energy to a given
 *    tolerance for arbitrary collections of charges.
 *
 * The implementation follows formulas 6 and 7 from:
 *
 *   N. D. Drummond et al., Physical Review B 78 125106 (2008)
 *
 *   DOI:  https://doi.org/10.1103/PhysRevB.78.125106
 */

#ifndef QMCPLUSPLUS_EWALD_REF_H
#define QMCPLUSPLUS_EWALD_REF_H

#include <vector>
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/Tensor.h"

namespace qmcplusplus
{
namespace ewaldref
{
/// Reference Ewald implemented for 3D only
enum
{
  DIM = 3
};

/// Type for integers
using int_t = int;
/// Type for floating point numbers
using real_t = double;
/// Type for integer vectors of length DIM
using IntVec = TinyVector<int_t, DIM>;
/// Type for floating point vectors of length DIM
using RealVec = TinyVector<real_t, DIM>;
/// Type for floating point matrices of shape DIM,DIM
using RealMat = Tensor<real_t, DIM>;
/// Type for lists of particle positions
using PosArray = std::vector<RealVec>;
/// Type for lists of particle charges
using ChargeArray = std::vector<real_t>;

/** Compute the total Ewald potential energy for a collection of charges
 *
 *  Corresponds to the entirety of Drummond 2008 formula 5, but for
 *    arbitrary charges.
 *
 *  @param a: Real-space cell axes.
 *
 *  @param R: List of particle coordinates.
 *
 *  @param R: List of particle charges.
 *
 *  @param tol: Tolerance for the total potential energy in Ha.
 */
real_t ewaldEnergy(const RealMat& a, const PosArray& R, const ChargeArray& Q, real_t tol = 1e-10);

} // namespace ewaldref
} // namespace qmcplusplus

#endif
