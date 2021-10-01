//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "EstimatorTesting.h"

namespace qmcplusplus
{
namespace testing
{
Lattice makeTestLattice()
{
  Lattice lattice;
  lattice.BoxBConds = true; // periodic
  lattice.R = ParticleSet::Tensor_t(3.37316115, 3.37316115, 0.00000000, 0.00000000, 3.37316115, 3.37316115, 3.37316115,
                                    0.00000000, 3.37316115);
  lattice.reset();
  lattice.explicitly_defined = true;
  return lattice;
}
} // namespace testing
} // namespace qmcplusplus
