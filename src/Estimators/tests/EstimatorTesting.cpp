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
#include <SpeciesSet.h>
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

SpeciesSet makeSpeciesSet()
{
  SpeciesSet species_set;
  int ispecies                      = species_set.addSpecies("C");
  int iattribute                    = species_set.addAttribute("membersize");
  species_set(iattribute, ispecies) = 2;
  return species_set;
}
} // namespace testing
} // namespace qmcplusplus
