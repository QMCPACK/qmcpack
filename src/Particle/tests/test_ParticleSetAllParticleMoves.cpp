//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//////////////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <stdexcept>

#include <catch2/catch_test_macros.hpp>
#include "Utilities/for_testing/Catch2Approx.h"

#include "Lattice/CrystalLattice.h"
#include "MinimalParticlePool.h"
#include "Particle/DistanceTable.h"
#include "Particle/LongRange/StructFact.h"
#include "Particle/PSdispatcher.h"
#include "Particle/ParticleSet.h"
#include "ResourceCollection.h"

namespace qmcplusplus
{
namespace
{
template<typename REAL>
REAL parseNonfinite(const char* value)
{
  if constexpr (sizeof(REAL) == sizeof(float))
    return std::strtof(value, nullptr);
  else
  {
    static_assert(sizeof(REAL) == sizeof(double));
    return std::strtod(value, nullptr);
  }
}

SimulationCell makeOpenCell()
{
  Lattice lattice;
  lattice.BoxBConds          = false;
  lattice.R                  = ParticleSet::Tensor_t(4.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 4.0);
  lattice.explicitly_defined = true;
  lattice.reset();
  return SimulationCell(lattice);
}

void checkPosition(const ParticleSet& pset, int iat, const ParticleSet::PosType& expected)
{
  for (int idim = 0; idim < OHMMS_DIM; ++idim)
  {
    CHECK(pset.R[iat][idim] == Approx(expected[idim]));
    CHECK(pset.getCoordinates().getAllParticlePos()[iat][idim] == Approx(expected[idim]));
  }
}

void checkDistances(const ParticleSet& elecs, int aa_id, int ab_id, const ParticleSet& ions)
{
  const auto distance = [](const ParticleSet::PosType& lhs, const ParticleSet::PosType& rhs) {
    const ParticleSet::PosType displacement = lhs - rhs;
    return std::sqrt(dot(displacement, displacement));
  };
  const auto& aa = elecs.getDistTableAA(aa_id).getDistances();
  const auto& ab = elecs.getDistTableAB(ab_id).getDistances();
  for (int iat = 0; iat < elecs.getTotalNum(); ++iat)
  {
    for (int jat = 0; jat < iat; ++jat)
      CHECK(aa[iat][jat] == Approx(distance(elecs.R[iat], elecs.R[jat])));
    for (int ion = 0; ion < ions.getTotalNum(); ++ion)
      CHECK(ab[iat][ion] == Approx(distance(elecs.R[iat], ions.R[ion])));
  }
}

class SingleWalkerAllParticleContext
{
public:
  SingleWalkerAllParticleContext(ParticleSet& pset, ParticleSet::Walker_t& walker)
      : p_list_(pset, {pset}), walkers_{walker}, resources_("single_walker_all_particle_resources")
  {
    pset.createResource(resources_);
    lock_ = std::make_unique<ResourceCollectionTeamLock<ParticleSet>>(resources_, p_list_);
  }

  template<CoordsType CT>
  bool makeMove(const MCCoords<CT>& displacements, bool skipSK = false)
  {
    std::vector<bool> valid(1);
    ParticleSet::mw_makeMoveAllParticles(p_list_, walkers_, displacements, valid, skipSK);
    return valid[0];
  }

  void resolve(bool accepted, bool skipSK = false)
  {
    ParticleSet::mw_accept_rejectMoveAllParticles(p_list_, walkers_, {accepted}, skipSK);
  }

private:
  RefVectorWithLeader<ParticleSet> p_list_;
  RefVector<ParticleSet::Walker_t> walkers_;
  ResourceCollection resources_;
  std::unique_ptr<ResourceCollectionTeamLock<ParticleSet>> lock_;
};

void exerciseDispatcherAllParticleMove(DynamicCoordinateKind kind)
{
  const SimulationCell cell = makeOpenCell();
  ParticleSet ions(cell);
  ions.setName("ion");
  ions.create({2});
  ions.R[0] = {1.0, 1.0, 1.0};
  ions.R[1] = {3.0, 3.0, 3.0};
  ions.update();

  ParticleSet p0(cell, kind);
  p0.setName("e");
  p0.create({2});
  p0.R[0]         = {0.5, 0.5, 0.5};
  p0.R[1]         = {1.5, 1.0, 0.5};
  const int ab_id = p0.addTable(ions);
  const int aa_id = p0.addTable(p0);
  p0.update();

  ParticleSet p1(p0);
  ParticleSet::Walker_t w0(2);
  ParticleSet::Walker_t w1(2);
  p0.saveWalker(w0);
  p1.saveWalker(w1);

  RefVectorWithLeader<ParticleSet> p_list(p0, {p0, p1});
  RefVector<ParticleSet::Walker_t> walkers{w0, w1};
  ResourceCollection resources("all_particle_move_resources");
  p0.createResource(resources);
  ResourceCollectionTeamLock<ParticleSet> lock(resources, p_list);

  // The flattened storage is particle-major: iat * nw + iw.
  MCCoords<CoordsType::POS> displacements(4);
  displacements.positions = {{0.2, 0.0, 0.0}, {0.1, 0.0, 0.0}, {0.0, 0.3, 0.0}, {0.0, -2.0, 0.0}};
  std::vector<bool> valid(2);
  PSdispatcher dispatcher(/*use_batch=*/true);
  dispatcher.flex_makeMoveAllParticles(p_list, walkers, displacements, valid);

  REQUIRE(valid == std::vector<bool>{true, false});
  checkPosition(p0, 0, {0.7, 0.5, 0.5});
  checkPosition(p0, 1, {1.5, 1.3, 0.5});
  checkPosition(p1, 0, w1.R[0]);
  checkPosition(p1, 1, w1.R[1]);
  CHECK(p0.getActivePtcl() == -1);
  CHECK(p1.getActivePtcl() == -1);
  checkDistances(p0, aa_id, ab_id, ions);
  checkDistances(p1, aa_id, ab_id, ions);

  dispatcher.flex_accept_rejectMoveAllParticles(p_list, walkers, {true, false});
  checkPosition(p0, 0, {0.7, 0.5, 0.5});
  checkPosition(p1, 1, w1.R[1]);
  checkDistances(p0, aa_id, ab_id, ions);
  checkDistances(p1, aa_id, ab_id, ions);

  // Saving the accepted transaction is required before a later rejection.
  dispatcher.flex_saveWalker(p_list, walkers);
  MCCoords<CoordsType::POS> next_displacements(4);
  next_displacements.positions = {{0.1, 0.0, 0.0}, {0.0, 0.0, 0.2}, {0.0, 0.1, 0.0}, {0.1, 0.0, 0.0}};
  dispatcher.flex_makeMoveAllParticles(p_list, walkers, next_displacements, valid);
  REQUIRE(valid == std::vector<bool>{true, true});
  dispatcher.flex_accept_rejectMoveAllParticles(p_list, walkers, {false, true});

  checkPosition(p0, 0, {0.7, 0.5, 0.5});
  checkPosition(p0, 1, {1.5, 1.3, 0.5});
  checkPosition(p1, 0, {0.5, 0.5, 0.7});
  checkPosition(p1, 1, {1.6, 1.0, 0.5});
  CHECK(p0.getActivePtcl() == -1);
  CHECK(p1.getActivePtcl() == -1);
  checkDistances(p0, aa_id, ab_id, ions);
  checkDistances(p1, aa_id, ab_id, ions);

  // Exercise all-accept and all-reject resolution, and walker synchronization.
  dispatcher.flex_saveWalker(p_list, walkers);
  for (size_t iw = 0; iw < p_list.size(); ++iw)
    for (int iat = 0; iat < p_list[iw].getTotalNum(); ++iat)
      CHECK(walkers[iw].get().R[iat] == p_list[iw].R[iat]);

  dispatcher.flex_makeMoveAllParticles(p_list, walkers, next_displacements, valid);
  REQUIRE(valid == std::vector<bool>{true, true});
  const ParticleSet::ParticlePos accepted_r0(p0.R);
  const ParticleSet::ParticlePos accepted_r1(p1.R);
  dispatcher.flex_accept_rejectMoveAllParticles(p_list, walkers, {true, true});
  for (int iat = 0; iat < p0.getTotalNum(); ++iat)
  {
    checkPosition(p0, iat, accepted_r0[iat]);
    checkPosition(p1, iat, accepted_r1[iat]);
  }

  dispatcher.flex_saveWalker(p_list, walkers);
  dispatcher.flex_makeMoveAllParticles(p_list, walkers, next_displacements, valid);
  REQUIRE(valid == std::vector<bool>{true, true});
  dispatcher.flex_accept_rejectMoveAllParticles(p_list, walkers, {false, false});
  for (int iat = 0; iat < p0.getTotalNum(); ++iat)
  {
    checkPosition(p0, iat, w0.R[iat]);
    checkPosition(p1, iat, w1.R[iat]);
  }
  checkDistances(p0, aa_id, ab_id, ions);
  checkDistances(p1, aa_id, ab_id, ions);
}
} // namespace

TEST_CASE("ParticleSet all-particle moves through the multiwalker dispatcher passthrough", "[particle]")
{
  for (const DynamicCoordinateKind kind : {DynamicCoordinateKind::DC_POS, DynamicCoordinateKind::DC_POS_OFFLOAD})
    exerciseDispatcherAllParticleMove(kind);
}

TEST_CASE("ParticleSet all-particle POS_SPIN rejection restores the resolved walker", "[particle]")
{
  ParticleSet pset(makeOpenCell());
  pset.create({2});
  pset.setSpinor(true);
  pset.R[0]  = {0.5, 0.5, 0.5};
  pset.R[1]  = {1.5, 1.0, 0.5};
  pset.spins = {0.25, -0.5};
  pset.update();

  ParticleSet::Walker_t walker(2);
  pset.saveWalker(walker);
  MCCoords<CoordsType::POS_SPIN> displacements(2);
  displacements.positions = {{0.1, 0.0, 0.0}, {0.0, 0.2, 0.0}};
  displacements.spins     = {0.5, -0.25};
  SingleWalkerAllParticleContext transaction(pset, walker);

  REQUIRE(transaction.makeMove(displacements));
  CHECK(pset.spins[0] == Approx(0.75));
  CHECK(pset.spins[1] == Approx(-0.75));
  transaction.resolve(false);

  checkPosition(pset, 0, walker.R[0]);
  checkPosition(pset, 1, walker.R[1]);
  CHECK(pset.spins[0] == Approx(walker.spins[0]));
  CHECK(pset.spins[1] == Approx(walker.spins[1]));
  CHECK(pset.getActivePtcl() == -1);

  displacements.positions = {{}, {}};
  displacements.spins     = {0.0, parseNonfinite<ParticleSet::RealType>("inf")};
  REQUIRE_FALSE(transaction.makeMove(displacements));
  CHECK(pset.spins[0] == Approx(walker.spins[0]));
  CHECK(pset.spins[1] == Approx(walker.spins[1]));
  checkPosition(pset, 0, walker.R[0]);
  checkPosition(pset, 1, walker.R[1]);
}

TEST_CASE("ParticleSet batched all-particle POS_SPIN resolves mixed outcomes", "[particle]")
{
  ParticleSet p0(makeOpenCell());
  p0.create({2});
  p0.R[0]  = {0.5, 0.5, 0.5};
  p0.R[1]  = {1.5, 1.0, 0.5};
  p0.spins = {0.25, -0.5};
  p0.update();
  ParticleSet p1(p0);

  ParticleSet::Walker_t w0(2);
  ParticleSet::Walker_t w1(2);
  p0.saveWalker(w0);
  p1.saveWalker(w1);
  RefVectorWithLeader<ParticleSet> p_list(p0, {p0, p1});
  RefVector<ParticleSet::Walker_t> walkers{w0, w1};
  ResourceCollection resources("all_particle_spin_move_resources");
  p0.createResource(resources);
  ResourceCollectionTeamLock<ParticleSet> lock(resources, p_list);

  MCCoords<CoordsType::POS_SPIN> displacements(4);
  displacements.positions = {{0.1, 0.0, 0.0}, {0.2, 0.0, 0.0}, {0.0, 0.1, 0.0}, {0.0, 0.2, 0.0}};
  displacements.spins     = {0.5, 1.0, -0.25, -0.75};
  std::vector<bool> valid(2);
  ParticleSet::mw_makeMoveAllParticles(p_list, walkers, displacements, valid);
  REQUIRE(valid == std::vector<bool>{true, true});
  CHECK(p0.spins[0] == Approx(0.75));
  CHECK(p0.spins[1] == Approx(-0.75));
  CHECK(p1.spins[0] == Approx(1.25));
  CHECK(p1.spins[1] == Approx(-1.25));

  ParticleSet::mw_accept_rejectMoveAllParticles(p_list, walkers, {true, false});
  CHECK(p0.spins[0] == Approx(0.75));
  CHECK(p0.spins[1] == Approx(-0.75));
  CHECK(p1.spins[0] == Approx(w1.spins[0]));
  CHECK(p1.spins[1] == Approx(w1.spins[1]));
  checkPosition(p1, 0, w1.R[0]);
  checkPosition(p1, 1, w1.R[1]);
}

TEST_CASE("ParticleSet all-particle periodic proposal remains unwrapped and supports absent SK", "[particle]")
{
  Lattice lattice;
  lattice.BoxBConds          = true;
  lattice.R                  = ParticleSet::Tensor_t(2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0);
  lattice.explicitly_defined = true;
  lattice.reset();
  ParticleSet pset{SimulationCell(lattice)};
  pset.create({1});
  pset.R[0] = {1.8, 0.5, 0.5};
  pset.update();
  REQUIRE_FALSE(pset.hasSK());

  ParticleSet::Walker_t walker(1);
  pset.saveWalker(walker);
  MCCoords<CoordsType::POS> displacement(1);
  displacement.positions = {{0.5, 0.0, 0.0}};
  SingleWalkerAllParticleContext transaction(pset, walker);
  REQUIRE(transaction.makeMove(displacement, true));
  checkPosition(pset, 0, {2.3, 0.5, 0.5});
  transaction.resolve(false, true);
  checkPosition(pset, 0, walker.R[0]);
}

TEST_CASE("ParticleSet all-particle move honors structure-factor skipping", "[particle]")
{
  auto particle_pool = MinimalParticlePool::make_diamondC_1x1x1(OHMMS::Controller);
  ParticleSet& pset  = *particle_pool.getParticleSet("e");
  REQUIRE(pset.hasSK());
  pset.update();

  ParticleSet::Walker_t walker(pset.getTotalNum());
  pset.saveWalker(walker);
  const auto resolved_rhok_r = pset.getSK().rhok_r;
  const auto resolved_rhok_i = pset.getSK().rhok_i;

  MCCoords<CoordsType::POS> displacement(pset.getTotalNum());
  std::fill(displacement.positions.begin(), displacement.positions.end(), ParticleSet::PosType{});
  displacement.positions[0] = {0.2, 0.1, 0.0};
  SingleWalkerAllParticleContext transaction(pset, walker);
  REQUIRE(transaction.makeMove(displacement, true));
  for (size_t i = 0; i < resolved_rhok_r.size(); ++i)
  {
    CHECK(pset.getSK().rhok_r.data()[i] == Approx(resolved_rhok_r.data()[i]));
    CHECK(pset.getSK().rhok_i.data()[i] == Approx(resolved_rhok_i.data()[i]));
  }

  transaction.resolve(false, true);
  REQUIRE(transaction.makeMove(displacement));
  bool any_rhok_changed = false;
  for (size_t i = 0; i < resolved_rhok_r.size(); ++i)
    any_rhok_changed = any_rhok_changed || pset.getSK().rhok_r.data()[i] != Approx(resolved_rhok_r.data()[i]) ||
        pset.getSK().rhok_i.data()[i] != Approx(resolved_rhok_i.data()[i]);
  CHECK(any_rhok_changed);

  transaction.resolve(false);
  for (int iat = 0; iat < pset.getTotalNum(); ++iat)
    checkPosition(pset, iat, walker.R[iat]);
  for (size_t i = 0; i < resolved_rhok_r.size(); ++i)
  {
    CHECK(pset.getSK().rhok_r.data()[i] == Approx(resolved_rhok_r.data()[i]));
    CHECK(pset.getSK().rhok_i.data()[i] == Approx(resolved_rhok_i.data()[i]));
  }
}

TEST_CASE("ParticleSet batched all-particle move updates and restores structure factors", "[particle]")
{
  auto particle_pool = MinimalParticlePool::make_diamondC_1x1x1(OHMMS::Controller);
  ParticleSet& p0    = *particle_pool.getParticleSet("e");
  REQUIRE(p0.hasSK());
  p0.update();
  ParticleSet p1(p0);

  ParticleSet::Walker_t w0(p0.getTotalNum());
  ParticleSet::Walker_t w1(p1.getTotalNum());
  p0.saveWalker(w0);
  p1.saveWalker(w1);
  const auto resolved_rhok_r = p1.getSK().rhok_r;
  const auto resolved_rhok_i = p1.getSK().rhok_i;

  RefVectorWithLeader<ParticleSet> p_list(p0, {p0, p1});
  RefVector<ParticleSet::Walker_t> walkers{w0, w1};
  ResourceCollection resources("all_particle_sk_move_resources");
  p0.createResource(resources);
  ResourceCollectionTeamLock<ParticleSet> lock(resources, p_list);

  const size_t num_walkers   = p_list.size();
  const size_t num_particles = p0.getTotalNum();
  MCCoords<CoordsType::POS> displacements(num_particles * num_walkers);
  std::fill(displacements.positions.begin(), displacements.positions.end(), ParticleSet::PosType{});
  displacements.positions[0] = {0.2, 0.1, 0.0};
  displacements.positions[1] = {-0.1, 0.2, 0.0};
  std::vector<bool> valid(num_walkers);
  ParticleSet::mw_makeMoveAllParticles(p_list, walkers, displacements, valid);
  REQUIRE(valid == std::vector<bool>{true, true});
  const auto proposed_rhok_r = p0.getSK().rhok_r;
  const auto proposed_rhok_i = p0.getSK().rhok_i;

  ParticleSet::mw_accept_rejectMoveAllParticles(p_list, walkers, {true, false});
  for (size_t i = 0; i < resolved_rhok_r.size(); ++i)
  {
    CHECK(p0.getSK().rhok_r.data()[i] == Approx(proposed_rhok_r.data()[i]));
    CHECK(p0.getSK().rhok_i.data()[i] == Approx(proposed_rhok_i.data()[i]));
    CHECK(p1.getSK().rhok_r.data()[i] == Approx(resolved_rhok_r.data()[i]));
    CHECK(p1.getSK().rhok_i.data()[i] == Approx(resolved_rhok_i.data()[i]));
  }
}

TEST_CASE("ParticleSet all-particle moves reject nonfinite proposals atomically", "[particle]")
{
  ParticleSet pset(makeOpenCell());
  pset.create({2});
  pset.R[0] = {0.5, 0.5, 0.5};
  pset.R[1] = {1.5, 0.5, 0.5};
  pset.update();
  ParticleSet::Walker_t walker(2);
  pset.saveWalker(walker);

  MCCoords<CoordsType::POS> displacements(2);
  displacements.positions = {{0.2, 0.0, 0.0}, {0.0, parseNonfinite<ParticleSet::RealType>("nan"), 0.0}};
  SingleWalkerAllParticleContext transaction(pset, walker);
  REQUIRE_FALSE(transaction.makeMove(displacements));
  checkPosition(pset, 0, walker.R[0]);
  checkPosition(pset, 1, walker.R[1]);

  displacements.positions[1] = {0.0, parseNonfinite<ParticleSet::RealType>("inf"), 0.0};
  REQUIRE_FALSE(transaction.makeMove(displacements));
  checkPosition(pset, 0, walker.R[0]);
  checkPosition(pset, 1, walker.R[1]);
}

#ifndef NDEBUG
TEST_CASE("ParticleSet all-particle moves reject malformed and active transactions", "[particle]")
{
  ParticleSet pset(makeOpenCell());
  pset.create({2});
  pset.R[0] = {0.5, 0.5, 0.5};
  pset.R[1] = {1.5, 0.5, 0.5};
  pset.update();
  ParticleSet::Walker_t walker(2);
  pset.saveWalker(walker);
  SingleWalkerAllParticleContext transaction(pset, walker);

  MCCoords<CoordsType::POS> malformed(1);
  malformed.positions = {{0.1, 0.0, 0.0}};
  REQUIRE_THROWS_AS(transaction.makeMove(malformed), std::runtime_error);

  MCCoords<CoordsType::POS> displacements(2);
  displacements.positions = {{0.1, 0.0, 0.0}, {0.0, 0.1, 0.0}};
  pset.makeMove(0, {0.1, 0.0, 0.0});
  REQUIRE_THROWS_AS(transaction.makeMove(displacements), std::runtime_error);

  ParticleSet malformed_pset(makeOpenCell());
  malformed_pset.create({1});
  ParticleSet::Walker_t malformed_walker(2);
  SingleWalkerAllParticleContext malformed_transaction(malformed_pset, malformed_walker);
  REQUIRE_THROWS_AS(malformed_transaction.resolve(false), std::runtime_error);
}

TEST_CASE("ParticleSet batched all-particle move rejects mismatched crowd topology before mutation", "[particle]")
{
  ParticleSet p0(makeOpenCell());
  p0.setName("e");
  p0.create({2});
  p0.R[0] = {0.5, 0.5, 0.5};
  p0.R[1] = {1.5, 0.5, 0.5};
  p0.update();
  ParticleSet p1(p0);
  p0.addTable(p0);

  ParticleSet::Walker_t w0(2);
  ParticleSet::Walker_t w1(2);
  p0.saveWalker(w0);
  p1.saveWalker(w1);
  RefVectorWithLeader<ParticleSet> p_list(p0, {p0, p1});
  RefVector<ParticleSet::Walker_t> walkers{w0, w1};
  MCCoords<CoordsType::POS> displacements(4);
  displacements.positions = {{0.1, 0.0, 0.0}, {0.2, 0.0, 0.0}, {0.0, 0.1, 0.0}, {0.0, 0.2, 0.0}};
  std::vector<bool> valid(2);

  REQUIRE_THROWS_AS(ParticleSet::mw_makeMoveAllParticles(p_list, walkers, displacements, valid), std::runtime_error);
  checkPosition(p0, 0, w0.R[0]);
  checkPosition(p0, 1, w0.R[1]);
  checkPosition(p1, 0, w1.R[0]);
  checkPosition(p1, 1, w1.R[1]);
}
#endif

} // namespace qmcplusplus
