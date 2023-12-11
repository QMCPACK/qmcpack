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

SpeciesSet makeSpeciesSet(const SpeciesCases species_case)
{
  SpeciesSet species_set;
  switch (species_case)
  {
  case SpeciesCases::GOOD: {
    species_set.addSpecies("u");
    species_set.addSpecies("d");
    const int iattribute       = species_set.addAttribute("membersize");
    species_set(iattribute, 0) = 1;
    species_set(iattribute, 1) = 1;
    break;
  }
  case SpeciesCases::NO_MEMBERSIZE: {
    species_set.addSpecies("u");
    species_set.addSpecies("d");
    break;
  }
  }
  return species_set;
}

OEBAccessor::OEBAccessor(OperatorEstBase& oeb) : oeb_(oeb) {}

OEBAccessor::value_type& OEBAccessor::operator[](size_t pos) { return oeb_.data_[pos]; }

} // namespace testing
} // namespace qmcplusplus
