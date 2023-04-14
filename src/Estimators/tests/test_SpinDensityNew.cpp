//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "SpinDensityInput.h"
#include "ValidSpinDensityInput.h"
#include "SpinDensityNew.h"
#include "hdf/hdf_archive.h"
#include "hdf/hdf_path.h"
#include "RandomForTest.h"
#include "ParticleSet.h"
#include "TrialWaveFunction.h"
#include "EstimatorTesting.h"
#include "NativeInitializerPrint.hpp"

#include "OhmmsData/Libxml2Doc.h"

#include <stdio.h>
#include <sstream>

namespace qmcplusplus
{

using QMCT         = QMCTraits;
using Real         = QMCTraits::RealType;
using FullPrecReal = QMCTraits::FullPrecRealType;

namespace testing
{
/** class to preserve access control in SpinDensityNew
 */
class SpinDensityNewTests
{
public:
  void testCopyConstructor(const SpinDensityNew& sdn)
  {
    SpinDensityNew sdn2(sdn);

    CHECK(sdn.species_size_ == sdn2.species_size_);
    CHECK(sdn.data_ != sdn2.data_);
  }

  const std::string& getName(const SpinDensityNew& sdn) { return sdn.my_name_; }
  const SpinDensityInput::DerivedParameters& getDerivedParameters(const SpinDensityNew& sdn)
  {
    return sdn.derived_parameters_;
  }
  static void report(const std::string& padding, const SpinDensityNew& sdn) { sdn.report(padding, std::cout); }
};
} // namespace testing

void accumulateFromPsets(int n_crowds,
                         int n_walkers,
                         SpinDensityNew& sdn,
                         UPtrVector<OperatorEstBase>& crowd_sdns,
                         int n_electrons,
                         PosUnit pos_unit = PosUnit::Lattice)
{
  // aka PtclOnLatticeTraits::ParticleLayout CrystalLattice<OHMMS_PRECISION, OHMMS_DIM>
  CrystalLattice<FullPrecReal, 3> lattice;
  lattice.set(Tensor<FullPrecReal, 3>(3.37316115, 3.37316115, 0.00000000, 0.00000000, 3.37316115, 3.37316115,
                                      3.37316115, 0.00000000, 3.37316115));
  const SimulationCell simulation_cell(lattice);
  for (int iops = 0; iops < n_crowds; ++iops)
  {
    std::vector<OperatorEstBase::MCPWalker> walkers;
    for (int iw = 0; iw < n_walkers; ++iw)
      walkers.emplace_back(n_electrons);

    std::vector<ParticleSet> psets;

    crowd_sdns.emplace_back(sdn.spawnCrowdClone());
    SpinDensityNew& crowd_sdn = dynamic_cast<SpinDensityNew&>(*(crowd_sdns.back()));

    for (int iw = 0; iw < n_walkers; ++iw)
    {
      psets.emplace_back(simulation_cell);
      ParticleSet& pset = psets.back();
      pset.create({n_electrons});
      pset.R.setUnit(pos_unit);
      Real inv_ne         = 1 / static_cast<Real>(n_electrons);
      const auto& lattice = simulation_cell.getLattice();
      std::cout << "inv_ne: " << inv_ne << '\n';
      // There is a strong assumption that the electron coordinates will only be from -0.5 to 0.5 in lattice coords
      for (int ie = 0; ie < n_electrons; ++ie)
      {
        ParticleSet::SingleParticlePos pos(ie * inv_ne, ie * inv_ne, ie * inv_ne);
        std::cout << "position" << pos[0] << " " << pos[1] << " " << pos[2] << '\n';
        ParticleSet::SingleParticlePos corner{-0.5, -0.5, -0.5};
        ParticleSet::SingleParticlePos adjusted_pos = pos - corner;
        switch (pos_unit)
        {
        case PosUnit::Lattice:
          pset.R[ie] = adjusted_pos;
          break;
        case PosUnit::Cartesian:
          std::cout << "Cartesian R: " << lattice.R << '\n';
          pset.R[ie] = lattice.toCart(adjusted_pos);
          std::cout << "cart position:" << pset.R[ie][0] << " " << pset.R[ie][1] << " " << pset.R[ie][2] << '\n';
          break;
        }
      }
    }

    std::vector<TrialWaveFunction> wfns;

    auto ref_walkers = makeRefVector<OperatorEstBase::MCPWalker>(walkers);
    auto ref_psets   = makeRefVector<ParticleSet>(psets);
    auto ref_wfns    = makeRefVector<TrialWaveFunction>(wfns);

    RandomGenerator rng;

    crowd_sdn.accumulate(ref_walkers, ref_psets, ref_wfns, rng);
  }
}

void randomUpdateAccumulate(int n_walkers,
                            testing::RandomForTest<QMCT::RealType>& rft,
                            UPtrVector<OperatorEstBase>& crowd_sdns,
                            int n_elecs)
{
  // aka PtclOnLatticeTraits::ParticleLayout CrystalLattice<OHMMS_PRECISION, OHMMS_DIM>
  CrystalLattice<FullPrecReal, 3> lattice;
  lattice.set(Tensor<FullPrecReal, 3>(3.37316115, 3.37316115, 0.00000000, 0.00000000, 3.37316115, 3.37316115,
                                      3.37316115, 0.00000000, 3.37316115));

  const SimulationCell simulation_cell(lattice);
  for (auto& uptr_crowd_sdn : crowd_sdns)
  {
    std::vector<OperatorEstBase::MCPWalker> walkers;
    for (int iw = 0; iw < n_walkers; ++iw)
      walkers.emplace_back(n_elecs);

    std::vector<ParticleSet> psets;

    SpinDensityNew& crowd_sdn = dynamic_cast<SpinDensityNew&>(*(uptr_crowd_sdn));

    std::vector<QMCT::RealType> rng_reals(n_walkers * QMCT::DIM * n_elecs);
    rft.fillVecRng(rng_reals);
    auto it_rng_reals = rng_reals.begin();
    for (int iw = 0; iw < n_walkers; ++iw)
    {
      psets.emplace_back(simulation_cell);
      ParticleSet& pset = psets.back();
      pset.create({n_elecs});
      for (int i = 0; i < n_elecs; ++i)
        pset.R[i] = ParticleSet::PosType(*it_rng_reals++, *it_rng_reals++, *it_rng_reals++);
    }

    std::vector<TrialWaveFunction> wfns;

    auto ref_walkers = makeRefVector<OperatorEstBase::MCPWalker>(walkers);
    auto ref_psets   = makeRefVector<ParticleSet>(psets);
    auto ref_wfns    = makeRefVector<TrialWaveFunction>(wfns);

    RandomGenerator rng;

    crowd_sdn.accumulate(ref_walkers, ref_psets, ref_wfns, rng);
  }
}

TEST_CASE("SpinDensityNew::SpinDensityNew(SPInput, SpeciesSet)", "[estimators]")
{
  Libxml2Document doc;
  bool okay = doc.parseFromString(testing::valid_spin_density_input_sections[testing::valid_spindensity_input_grid]);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();
  SpinDensityInput sdi(node);
  SpeciesSet species_set;
  int ispecies                      = species_set.addSpecies("C");
  int iattribute                    = species_set.addAttribute("membersize");
  species_set(iattribute, ispecies) = 2;
  SpinDensityInput sdi_copy         = sdi;
  SpinDensityNew(std::move(sdi), species_set);
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
}

TEST_CASE("SpinDensityNew::SpinDensityNew(SPInput, Lattice, SpeciesSet)", "[estimators]")
{
  Libxml2Document doc;
  bool okay = doc.parseFromString(testing::valid_spin_density_input_sections[testing::valid_spindensity_input_no_cell]);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();
  SpinDensityInput sdi(node);
  SpeciesSet species_set;
  int ispecies                      = species_set.addSpecies("C");
  int iattribute                    = species_set.addAttribute("membersize");
  species_set(iattribute, ispecies) = 2;
  auto lattice                      = testing::makeTestLattice();
  SpinDensityNew sdn(std::move(sdi), lattice, species_set);
  // make sure there is something in obdm's data
  using namespace testing;
  OEBAccessor oeba(sdn);
  oeba[0] = 1.0;
  SpinDensityNewTests sdnt;
  sdnt.testCopyConstructor(sdn);
}

TEST_CASE("SpinDensityNew::spawnCrowdClone()", "[estimators]")
{
  Libxml2Document doc;
  bool okay = doc.parseFromString(testing::valid_spin_density_input_sections[testing::valid_spindensity_input_no_cell]);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();
  SpinDensityInput sdi(node);
  SpeciesSet species_set;
  int ispecies                      = species_set.addSpecies("C");
  int iattribute                    = species_set.addAttribute("membersize");
  species_set(iattribute, ispecies) = 2;
  auto lattice                      = testing::makeTestLattice();
  SpinDensityNew original(std::move(sdi), lattice, species_set);
  auto clone = original.spawnCrowdClone();
  REQUIRE(clone != nullptr);
  REQUIRE(clone.get() != &original);
  REQUIRE(dynamic_cast<decltype(&original)>(clone.get()) != nullptr);
}

TEST_CASE("SpinDensityNew::accumulate", "[estimators]")
{
  using MCPWalker = OperatorEstBase::MCPWalker;
  using QMCT      = QMCTraits;

  Libxml2Document doc;
  bool okay = doc.parseFromString(testing::valid_spin_density_input_sections[0]);
  REQUIRE(okay);
  xmlNodePtr node = doc.getRoot();
  SpinDensityInput sdi(node);
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

  const SimulationCell simulation_cell;
  for (int iw = 0; iw < nwalkers; ++iw)
  {
    psets.emplace_back(simulation_cell);
    ParticleSet& pset = psets.back();
    pset.create({2});
    pset.R[0] = ParticleSet::PosType(0.00000000, 0.00000000, 0.00000000);
    pset.R[1] = ParticleSet::PosType(1.68658058, 1.68658058, 1.68658058);
  }

  std::vector<TrialWaveFunction> wfns;

  auto ref_walkers = makeRefVector<MCPWalker>(walkers);
  auto ref_psets   = makeRefVector<ParticleSet>(psets);
  auto ref_wfns    = makeRefVector<TrialWaveFunction>(wfns);

  RandomGenerator rng;

  sdn.accumulate(ref_walkers, ref_psets, ref_wfns, rng);

  std::vector<QMCT::RealType>& data_ref = sdn.get_data();
  // There should be a check that the discretization of particle locations expressed in lattice coords
  // is correct.  This just checks it hasn't changed from how it was in SpinDensity which lacked testing.
}

TEST_CASE("SpinDensityNew::collect(DataLocality::crowd)", "[estimators]")
{
  auto checkCollectForNParticles = [](int n_elec) {
    using MCPWalker = OperatorEstBase::MCPWalker;
    using QMCT      = QMCTraits;

    Libxml2Document doc;
    bool okay = doc.parseFromString(testing::valid_spin_density_input_sections[0]);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    SpinDensityInput sdi(node);
    SpeciesSet species_set;
    int ispecies = species_set.addSpecies("u");
    ispecies     = species_set.addSpecies("d");
    CHECK(ispecies == 1);
    int i_msize             = species_set.addAttribute("membersize");
    species_set(i_msize, 0) = n_elec / 2 + n_elec % 2;
    species_set(i_msize, 1) = n_elec / 2;

    SpinDensityNew sdn(std::move(sdi), species_set);

    UPtrVector<OperatorEstBase> crowd_sdns;
    int n_crowds  = 2;
    int n_walkers = 4;

    accumulateFromPsets(n_crowds, n_walkers, sdn, crowd_sdns, n_elec);

    RefVector<OperatorEstBase> crowd_oeb_refs = convertUPtrToRefVector(crowd_sdns);
    sdn.collect(crowd_oeb_refs);

    std::vector<QMCT::RealType>& data_ref = sdn.get_data();
    // There should be a check that the discretization of particle locations expressed in lattice coords
    // is correct.  This just checks it hasn't changed from how it was in SpinDensity which lacked testing.
    Real up_sum{0};
    Real dn_sum{0};
    for (int i = 0; i < 1000; ++i)
      up_sum += data_ref[i];
    for (int i = 0; i < 1000; ++i)
      dn_sum += data_ref[i + 1000];
    // \todo why is it 8?
    CHECK(up_sum == n_crowds * n_walkers * species_set(i_msize, 0));
    CHECK(dn_sum == n_crowds * n_walkers * species_set(i_msize, 1));
  };
  checkCollectForNParticles(1);
  checkCollectForNParticles(2);
  checkCollectForNParticles(3);
  checkCollectForNParticles(11);
}

TEST_CASE("SpinDensityNew::collect(DataLocality::rank)", "[estimators]")
{
  using MCPWalker = OperatorEstBase::MCPWalker;
  using QMCT      = QMCTraits;

  auto checkCollectAtPoint = [](int n_elec, PosUnit pos_unit) {
    Libxml2Document doc;
    bool okay = doc.parseFromString(testing::valid_spin_density_input_sections[0]);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    SpinDensityInput sdi(node);
    SpeciesSet species_set;
    int ispecies = species_set.addSpecies("u");
    ispecies     = species_set.addSpecies("d");
    CHECK(ispecies == 1);
    int i_msize = species_set.addAttribute("membersize");
    std::vector<int> species_size(species_set.size());
    species_size[0]         = n_elec / 2 + n_elec % 2;
    species_size[1]         = n_elec / 2;
    species_set(i_msize, 0) = species_size[0];
    species_set(i_msize, 1) = species_size[1];


    SpinDensityNew sdn(std::move(sdi), species_set, DataLocality::rank);

    testing::SpinDensityNewTests sdnt;
    const auto& derived_parameters = sdnt.getDerivedParameters(sdn);

    auto lattice = testing::makeTestLattice();

    UPtrVector<OperatorEstBase> crowd_sdns;
    int n_crowds  = 2;
    int n_walkers = 4;

    accumulateFromPsets(n_crowds, n_walkers, sdn, crowd_sdns, n_elec, pos_unit);

    RefVector<OperatorEstBase> crowd_oeb_refs = convertUPtrToRefVector(crowd_sdns);
    sdn.collect(crowd_oeb_refs);

    std::vector<QMCT::RealType>& data_ref = sdn.get_data();
    // There should be a check that the discretization of particle locations expressed in lattice coords
    // is correct.  This just checks it hasn't changed from how it was in SpinDensity which lacked testing.
    //
    INFO("n_elec: " << n_elec);
    CHECK(data_ref.size() == 2000);
    sdnt.report("", sdn);
    std::cout << NativePrint(data_ref) << '\n';


    size_t offset = 0;
    int p         = 0;
    Real inv_ne   = 1 / static_cast<Real>(n_elec);
    for (int s = 0; s < species_set.size(); ++s, offset += derived_parameters.npoints)
      for (int ps = 0; ps < species_size[s]; ++ps, ++p)
      {
        //QMCT::PosType u = lattice_.toUnit(pset.R[p] - dp_.corner);
        ParticleSet::SingleParticlePos pos(p * inv_ne, p * inv_ne, p * inv_ne);
        std::cout << "position" << pos[0] << " " << pos[1] << " " << pos[2] << '\n';
        //pos = lattice.toCart(pos);
        //std::cout << "position" << pos[0] << " " << pos[1] << " " << pos[2] <<'\n';

        //pset.R[ie] = pos;
        size_t point = offset;
        for (int d = 0; d < QMCT::DIM; ++d)
          point += derived_parameters.gdims[d] *
              ((int)(derived_parameters.grid[d] * (pos[d] - std::floor(pos[d])))); //periodic only
        INFO("species:" << s << " ps: " << ps << " point:" << point);
        CHECK(data_ref[point] == n_walkers * n_crowds);
      }
  };
  checkCollectAtPoint(2, PosUnit::Lattice);
  checkCollectAtPoint(3, PosUnit::Lattice);
  checkCollectAtPoint(2, PosUnit::Cartesian);
  checkCollectAtPoint(3, PosUnit::Cartesian);
}

TEST_CASE("SpinDensityNew algorithm comparison", "[estimators]")
{
  using MCPWalker = OperatorEstBase::MCPWalker;
  using QMCT      = QMCTraits;

  auto checkAlgorithmsForNParticles = [](int n_elec) {
    Libxml2Document doc;
    bool okay = doc.parseFromString(testing::valid_spin_density_input_sections[0]);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    SpinDensityInput sdi(node);
    SpeciesSet species_set;
    int ispecies = species_set.addSpecies("u");
    ispecies     = species_set.addSpecies("d");
    CHECK(ispecies == 1);
    int i_msize             = species_set.addAttribute("membersize");
    species_set(i_msize, 0) = n_elec / 2 + n_elec % 2;
    species_set(i_msize, 1) = n_elec / 2;

    SpinDensityInput sdi_copy = sdi;

    int n_crowds  = 3;
    int n_walkers = 4;
    int nsteps    = 4;

    SpinDensityNew sdn_rank(std::move(sdi), species_set, DataLocality::rank);
    UPtrVector<OperatorEstBase> crowd_sdns_rank;
    accumulateFromPsets(n_crowds, n_walkers, sdn_rank, crowd_sdns_rank, n_elec);
    testing::RandomForTest<QMCT::RealType> rng_for_test_rank;
    for (int i = 0; i < nsteps; ++i)
      randomUpdateAccumulate(n_walkers, rng_for_test_rank, crowd_sdns_rank, n_elec);
    RefVector<OperatorEstBase> crowd_oeb_refs_rank = convertUPtrToRefVector(crowd_sdns_rank);
    sdn_rank.collect(crowd_oeb_refs_rank);
    std::vector<QMCT::RealType>& data_ref_rank = sdn_rank.get_data();
    SpinDensityNew sdn_crowd(std::move(sdi), species_set, DataLocality::crowd);
    UPtrVector<OperatorEstBase> crowd_sdns_crowd;
    accumulateFromPsets(n_crowds, n_walkers, sdn_crowd, crowd_sdns_crowd, n_elec);
    testing::RandomForTest<QMCT::RealType> rng_for_test_crowd;
    for (int i = 0; i < nsteps; ++i)
      randomUpdateAccumulate(n_walkers, rng_for_test_crowd, crowd_sdns_crowd, n_elec);
    RefVector<OperatorEstBase> crowd_oeb_refs_crowd = convertUPtrToRefVector(crowd_sdns_crowd);
    sdn_crowd.collect(crowd_oeb_refs_crowd);
    std::vector<QMCT::RealType>& data_ref_crowd = sdn_crowd.get_data();

    for (size_t i = 0; i < data_ref_rank.size(); ++i)
    {
      if (data_ref_crowd[i] != data_ref_rank[i])
        FAIL_CHECK("crowd local " << data_ref_crowd[i] << " != rank local " << data_ref_rank[i] << " at index " << i);
      break;
    }

    auto sumAndCheck = [i_msize, nsteps, n_crowds, n_walkers, &species_set](const auto& data_ref) {
      Real up_sum{0};
      Real dn_sum{0};
      for (int i = 0; i < 1000; ++i)
        up_sum += data_ref[i];
      for (int i = 0; i < 1000; ++i)
        dn_sum += data_ref[i + 1000];
      // +1 is for the accumulateFromPsets
      CHECK(up_sum == (nsteps + 1) * n_crowds * n_walkers * species_set(i_msize, 0));
      CHECK(dn_sum == (nsteps + 1) * n_crowds * n_walkers * species_set(i_msize, 1));
    };
    sumAndCheck(data_ref_rank);
    sumAndCheck(data_ref_crowd);
  };

  checkAlgorithmsForNParticles(2);
  checkAlgorithmsForNParticles(3);
}

TEST_CASE("SpinDensityNew::registerAndWrite", "[estimators]")
{
  auto checkWriteForNParticles = [](int n_elec) {
    using MCPWalker = OperatorEstBase::MCPWalker;
    using QMCT      = QMCTraits;

    Libxml2Document doc;
    bool okay = doc.parseFromString(testing::valid_spin_density_input_sections[0]);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    SpinDensityInput sdi(node);
    SpeciesSet species_set;
    int ispecies = species_set.addSpecies("u");
    ispecies     = species_set.addSpecies("d");
    CHECK(ispecies == 1);
    int i_msize             = species_set.addAttribute("membersize");
    species_set(i_msize, 0) = n_elec / 2 + n_elec % 2;
    species_set(i_msize, 1) = n_elec / 2;

    SpinDensityNew sdn(std::move(sdi), species_set);

    UPtrVector<OperatorEstBase> crowd_sdns;
    int n_crowds  = 2;
    int n_walkers = 4;

    accumulateFromPsets(n_crowds, n_walkers, sdn, crowd_sdns, n_elec);

    RefVector<OperatorEstBase> crowd_oeb_refs = convertUPtrToRefVector(crowd_sdns);
    sdn.collect(crowd_oeb_refs);

    hdf_archive hd;
    std::string test_file{"spin_density_new_test.hdf"};
    okay = hd.create(test_file);
    REQUIRE(okay);

    sdn.registerOperatorEstimator(hd);
    sdn.write(hd);
    hd.close();

    hdf_archive read_hd;
    okay = read_hd.open(test_file);
    REQUIRE(okay);
    std::vector<FullPrecReal> up_data;
    std::vector<FullPrecReal> dn_data;

    testing::SpinDensityNewTests sdnt;
    hdf_path sdn_name_up{sdnt.getName(sdn)};

    sdn_name_up /= "u/value";
    std::cout << sdn_name_up.string() << '\n';
    std::array<int, 2> select_one_grid{0, -1};
    read_hd.readSlabSelection(up_data, select_one_grid, sdn_name_up.string()); //sdn_name_up.string());
    CHECK(up_data.size() == 1000);
    hdf_path sdn_name_dn{sdnt.getName(sdn)};
    sdn_name_dn /= "d/value";
    std::cout << sdn_name_dn.string() << '\n';
    read_hd.readSlabSelection(dn_data, select_one_grid, sdn_name_dn.string());
    CHECK(up_data.size() == 1000);

    auto sumAndCheckData = [n_crowds, n_walkers](const auto& data, int num_particles) {
      FullPrecReal sum = std::accumulate(data.begin(), data.end(), 0.0);
      // +1 is for the accumulateFromPsets
      INFO("num_particles: " << num_particles);
      CHECK(sum == n_crowds * n_walkers * num_particles);
    };

    std::cout << NativePrint(up_data);
    std::cout << NativePrint(dn_data);

    sumAndCheckData(up_data, species_set(i_msize, 0));
    sumAndCheckData(dn_data, species_set(i_msize, 1));
  };
  checkWriteForNParticles(2);
  checkWriteForNParticles(3);
}

} // namespace qmcplusplus
