//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "SpinDensityInput.h"
#include "ValidSpinDensityInput.h"
#include "SpinDensityNew.h"
#include "Utilities/SpeciesSet.h"
#include "ParticleSet.h"

#include "OhmmsData/Libxml2Doc.h"

#include <stdio.h>
#include <sstream>

namespace qmcplusplus
{
TEST_CASE("SpinDensityNew", "[estimators]")
{
  Libxml2Document doc;
  bool okay = doc.parseFromString(testing::valid_spin_density_input_sections[0]);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();
  SpinDensityInput sdi;
  sdi.readXML(node);
  SpeciesSet species_set;
  int ispecies                      = species_set.addSpecies("C");
  int iattribute                    = species_set.addAttribute("membersize");
  species_set(iattribute, ispecies) = 2;
  SpinDensityNew(std::move(sdi), species_set);
}

TEST_CASE("SpinDensityNew::accumulate", "[estimators]")
{
  using MCPWalker = OperatorEstBase::MCPWalker;
  using QMCT      = QMCTraits;

  Libxml2Document doc;
  bool okay = doc.parseFromString(testing::valid_spin_density_input_sections[0]);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();
  SpinDensityInput sdi;
  sdi.readXML(node);
  SpeciesSet species_set;
  int ispecies = species_set.addSpecies("u");
  ispecies     = species_set.addSpecies("d");
  CHECK(ispecies == 1);
  int iattribute             = species_set.addAttribute("membersize");
  species_set(iattribute, 0) = 1;
  species_set(iattribute, 1) = 1;

  SpinDensityNew sdn(std::move(sdi), species_set);
  std::vector<MCPWalker> walkers;
  int nwalkers = 4;
  for (int iw = 0; iw < nwalkers; ++iw)
    walkers.emplace_back(2);

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> Lattice;
  Lattice.BoxBConds = true; // periodic
  Lattice.R = ParticleSet::Tensor_t(3.37316115, 3.37316115, 0.00000000, 0.00000000, 3.37316115, 3.37316115, 3.37316115,
                                    0.00000000, 3.37316115);
  Lattice.reset();

  std::vector<ParticleSet> psets;

  for (int iw = 0; iw < nwalkers; ++iw)
  {
    psets.emplace_back();
    ParticleSet& pset = psets.back();
    pset.create(2);
    pset.R[0] = ParticleSet::PosType(0.00000000, 0.00000000, 0.00000000);
    pset.R[1] = ParticleSet::PosType(1.68658058, 1.68658058, 1.68658058);
  }

  // walkers[0].R[0] = {0.5,0.5,0.5};
  // walkers[0].R[1] = {0.2,0.2,0.2};
  // walkers[1].R[0] = {0.5,0.5,0.5};
  // walkers[1].R[1] = {0.2,0.2,0.2};
  // walkers[2].R[0] = {0.5,0.5,0.5};
  // walkers[2].R[1] = {0.2,0.2,0.2};
  // walkers[3].R[0] = {0.5,0.5,0.5};
  // walkers[3].R[1] = {0.2,0.2,0.2};
  auto ref_walkers = makeRefVector<MCPWalker>(walkers);
  auto ref_psets   = makeRefVector<ParticleSet>(psets);

  sdn.accumulate(ref_walkers, ref_psets);

  std::vector<QMCT::RealType>& data_ref = sdn.get_data_ref();
  // There should be a check that the discretization of particle locations expressed in lattice coords
  // is correct.  This just checks it hasn't changed from how it was in SpinDensity which lacked testing.
  CHECK(data_ref[555] == 4);
  CHECK(data_ref[1777] == 4);
}

TEST_CASE("SpinDensityNew::collect", "[estimators]") {}
} // namespace qmcplusplus
