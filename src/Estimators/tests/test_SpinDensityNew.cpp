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
#include "RandomForTest.h"
#include "Utilities/SpeciesSet.h"
#include "ParticleSet.h"
#include "SpinDensityTesting.h"

#include "OhmmsData/Libxml2Doc.h"

#include <stdio.h>
#include <sstream>

namespace qmcplusplus
{
void accumulateFromPsets(int ncrowds, SpinDensityNew& sdn, UPtrVector<OperatorEstBase>& crowd_sdns)
{
  using QMCT = QMCTraits;

  for (int iops = 0; iops < ncrowds; ++iops)
  {
    std::vector<OperatorEstBase::MCPWalker> walkers;
    int nwalkers = 4;
    for (int iw = 0; iw < nwalkers; ++iw)
      walkers.emplace_back(2);

    std::vector<ParticleSet> psets;

    crowd_sdns.emplace_back(std::make_unique<SpinDensityNew>(sdn));
    SpinDensityNew& crowd_sdn = dynamic_cast<SpinDensityNew&>(*(crowd_sdns.back()));

    for (int iw = 0; iw < nwalkers; ++iw)
    {
      psets.emplace_back();
      ParticleSet& pset = psets.back();
      pset.create(2);
      pset.R[0] = ParticleSet::PosType(0.00000000, 0.00000000, 0.00000000);
      pset.R[1] = ParticleSet::PosType(0.68658058, 0.68658058, 0.68658058);
    }

    auto ref_walkers = makeRefVector<OperatorEstBase::MCPWalker>(walkers);
    auto ref_psets   = makeRefVector<ParticleSet>(psets);

    crowd_sdn.accumulate(ref_walkers, ref_psets);
  }
}

void randomUpdateAccumulate(testing::RandomForTest& rft, UPtrVector<OperatorEstBase>& crowd_sdns)
{
  using QMCT = QMCTraits;

  for (auto& uptr_crowd_sdn : crowd_sdns)
  {
    std::vector<OperatorEstBase::MCPWalker> walkers;
    int nwalkers = 4;
    for (int iw = 0; iw < nwalkers; ++iw)
      walkers.emplace_back(2);

    std::vector<ParticleSet> psets;

    SpinDensityNew& crowd_sdn = dynamic_cast<SpinDensityNew&>(*(uptr_crowd_sdn));

    std::vector<QMCT::RealType> rng_reals(nwalkers * QMCT::DIM * 2);
    rft.makeRngReals(rng_reals);
    auto it_rng_reals = rng_reals.begin();
    for (int iw = 0; iw < nwalkers; ++iw)
    {
      psets.emplace_back();
      ParticleSet& pset = psets.back();
      pset.create(2);
      pset.R[0] = ParticleSet::PosType(*it_rng_reals++, *it_rng_reals++, *it_rng_reals++);
      pset.R[1] = ParticleSet::PosType(*it_rng_reals++, *it_rng_reals++, *it_rng_reals++);
    }

    auto ref_walkers = makeRefVector<OperatorEstBase::MCPWalker>(walkers);
    auto ref_psets   = makeRefVector<ParticleSet>(psets);

    crowd_sdn.accumulate(ref_walkers, ref_psets);
  }
}

TEST_CASE("SpinDensityNew::SpinDensityNew(SPInput, SpeciesSet)", "[estimators]")
{
  Libxml2Document doc;
  bool okay = doc.parseFromString(testing::valid_spin_density_input_sections[testing::valid_spindensity_input_grid]);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();
  SpinDensityInput sdi;
  sdi.readXML(node);
  SpeciesSet species_set;
  int ispecies                      = species_set.addSpecies("C");
  int iattribute                    = species_set.addAttribute("membersize");
  species_set(iattribute, ispecies) = 2;
  SpinDensityInput sdi_copy         = sdi;
  SpinDensityNew(std::move(sdi), species_set);
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  CHECK_THROWS(SpinDensityNew(std::move(sdi_copy), lattice, species_set));
}

TEST_CASE("SpinDensityNew::SpinDensityNew(SPInput, Lattice, SpeciesSet)", "[estimators]")
{
  Libxml2Document doc;
  bool okay = doc.parseFromString(testing::valid_spin_density_input_sections[testing::valid_spindensity_input_no_cell]);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();
  SpinDensityInput sdi;
  sdi.readXML(node);
  SpeciesSet species_set;
  int ispecies                      = species_set.addSpecies("C");
  int iattribute                    = species_set.addAttribute("membersize");
  species_set(iattribute, ispecies) = 2;
  auto lattice                      = testing::makeTestLattice();
  SpinDensityNew(std::move(sdi), lattice, species_set);
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

  std::vector<ParticleSet> psets;

  for (int iw = 0; iw < nwalkers; ++iw)
  {
    psets.emplace_back();
    ParticleSet& pset = psets.back();
    pset.create(2);
    pset.R[0] = ParticleSet::PosType(0.00000000, 0.00000000, 0.00000000);
    pset.R[1] = ParticleSet::PosType(1.68658058, 1.68658058, 1.68658058);
  }

  auto ref_walkers = makeRefVector<MCPWalker>(walkers);
  auto ref_psets   = makeRefVector<ParticleSet>(psets);

  sdn.accumulate(ref_walkers, ref_psets);

  std::vector<QMCT::RealType>& data_ref = sdn.get_data_ref();
  // There should be a check that the discretization of particle locations expressed in lattice coords
  // is correct.  This just checks it hasn't changed from how it was in SpinDensity which lacked testing.
  CHECK(data_ref[555] == 4);
  CHECK(data_ref[1777] == 4);
}

TEST_CASE("SpinDensityNew::collect(DataLocality::crowd)", "[estimators]")
{
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

    UPtrVector<OperatorEstBase> crowd_sdns;
    int ncrowds = 2;

    accumulateFromPsets(ncrowds, sdn, crowd_sdns);

    RefVector<OperatorEstBase> crowd_oeb_refs = convertUPtrToRefVector(crowd_sdns);
    sdn.collect(crowd_oeb_refs);

    std::vector<QMCT::RealType>& data_ref = sdn.get_data_ref();
    // There should be a check that the discretization of particle locations expressed in lattice coords
    // is correct.  This just checks it hasn't changed from how it was in SpinDensity which lacked testing.
    CHECK(data_ref[555] == 4 * ncrowds);
    CHECK(data_ref[1666] == 4 * ncrowds);
  }
}

TEST_CASE("SpinDensityNew::collect(DataLocality::rank)", "[estimators]")
{
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

    SpinDensityNew sdn(std::move(sdi), species_set, DataLocality::rank);

    auto lattice = testing::makeTestLattice();

    UPtrVector<OperatorEstBase> crowd_sdns;
    int ncrowds = 2;

    accumulateFromPsets(ncrowds, sdn, crowd_sdns);

    RefVector<OperatorEstBase> crowd_oeb_refs = convertUPtrToRefVector(crowd_sdns);
    sdn.collect(crowd_oeb_refs);

    std::vector<QMCT::RealType>& data_ref = sdn.get_data_ref();
    // There should be a check that the discretization of particle locations expressed in lattice coords
    // is correct.  This just checks it hasn't changed from how it was in SpinDensity which lacked testing.
    CHECK(data_ref[555] == 4 * ncrowds);
    CHECK(data_ref[1666] == 4 * ncrowds);
  }
}

TEST_CASE("SpinDensityNew algorithm comparison", "[estimators]")
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

  SpinDensityInput sdi_copy = sdi;

  int ncrowds = 3;
  int nsteps  = 4;

  SpinDensityNew sdn_rank(std::move(sdi), species_set, DataLocality::rank);
  UPtrVector<OperatorEstBase> crowd_sdns_rank;
  accumulateFromPsets(ncrowds, sdn_rank, crowd_sdns_rank);
  testing::RandomForTest rng_for_test_rank;
  for (int i = 0; i < nsteps; ++i)
    randomUpdateAccumulate(rng_for_test_rank, crowd_sdns_rank);
  RefVector<OperatorEstBase> crowd_oeb_refs_rank = convertUPtrToRefVector(crowd_sdns_rank);
  sdn_rank.collect(crowd_oeb_refs_rank);
  std::vector<QMCT::RealType>& data_ref_rank = sdn_rank.get_data_ref();

  SpinDensityNew sdn_crowd(std::move(sdi), species_set, DataLocality::crowd);
  UPtrVector<OperatorEstBase> crowd_sdns_crowd;
  accumulateFromPsets(ncrowds, sdn_crowd, crowd_sdns_crowd);
  testing::RandomForTest rng_for_test_crowd;
  for (int i = 0; i < nsteps; ++i)
    randomUpdateAccumulate(rng_for_test_crowd, crowd_sdns_crowd);
  RefVector<OperatorEstBase> crowd_oeb_refs_crowd = convertUPtrToRefVector(crowd_sdns_crowd);
  sdn_crowd.collect(crowd_oeb_refs_crowd);
  std::vector<QMCT::RealType>& data_ref_crowd = sdn_crowd.get_data_ref();

  for (size_t i = 0; i < data_ref_rank.size(); ++i)
  {
    if (data_ref_crowd[i] != data_ref_rank[i])
      FAIL_CHECK("crowd local " << data_ref_crowd[i] << " != rank local " << data_ref_rank[i] << " at index " << i);
    break;
  }
}


} // namespace qmcplusplus
