//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include <functional>
#include <cmath>
#include <limits>
#include <catch2/catch_test_macros.hpp>
#include "Utilities/for_testing/Catch2Approx.h"

#include "test_WalkerControl.h"
#include "Message/Communicate.h"
#include "Message/UniformCommunicateError.h"
#include <catch2/matchers/catch_matchers_exception.hpp>
#include "ParticleSetPool.h"
#include "WaveFunctionPool.h"
#include "HamiltonianPool.h"
#include "QMCDrivers/MCPopulation.h"
#include "Utilities/MPIExceptionWrapper.hpp"
#include "Utilities/FairDivide.h"
#include "Utilities/for_testing/NativeInitializerPrint.hpp"


//#include "Concurrency/Info.hpp"
//#include "Concurrency/UtilityFunctions.hpp"

namespace qmcplusplus
{
namespace testing
{
UnifiedDriverWalkerControlMPITest::UnifiedDriverWalkerControlMPITest() : wc_(dpools_.comm, Random)
{
  int num_ranks = dpools_.comm.size();
  pop_ = std::make_unique<MCPopulation>(num_ranks, dpools_.comm.rank(), *dpools_.particle_pool->getParticleSet("e"),
                                        dpools_.wavefunction_pool->getWaveFunction().value(),
                                        dpools_.hamiltonian_pool->getHamiltonian().value());

  pop_->createWalkers(1, walker_confs);
}

void UnifiedDriverWalkerControlMPITest::reportWalkersPerRank(Communicate& c, MCPopulation& pop)
{
#if !defined(NDEBUG)
  std::vector<int> rank_walker_count(c.size(), 0);
  rank_walker_count[c.rank()] = pop.get_num_local_walkers();
  c.allreduce(rank_walker_count);

  const int current_population = std::accumulate(rank_walker_count.begin(), rank_walker_count.end(), 0);

  app_log() << "Walkers Per Rank (Total: " << current_population << ")" << std::endl;
  app_log() << "Rank   Count" << std::endl << "===========" << std::endl;
  for (int i = 0; i < rank_walker_count.size(); ++i)
    app_log() << std::setw(4) << i << "   " << rank_walker_count[i] << std::endl;
#endif
}

void UnifiedDriverWalkerControlMPITest::testNewDistribution(std::vector<int>& initial_num_per_rank,
                                                            std::vector<int>& minus,
                                                            std::vector<int>& plus)
{
  assert(initial_num_per_rank.size() == dpools_.comm.size());
  std::vector<int> num_per_rank = initial_num_per_rank;
  std::vector<int> fair_offset;
  wc_.determineNewWalkerPopulation(num_per_rank, fair_offset, minus, plus);
}

void UnifiedDriverWalkerControlMPITest::testInput()
{
  Libxml2Document omitted_doc;
  REQUIRE(omitted_doc.parseFromString(R"(<qmc method="dmc_batch"/>)"));
  REQUIRE(wc_.put(omitted_doc.getRoot()));
  CHECK(wc_.use_nonblocking_);
  CHECK_FALSE(wc_.debug_disable_branching_);

  const char* const configured_input = R"(
  <qmc method="dmc_batch">
    <parameter name="use_nonblocking">false</parameter>
    <parameter name="debug_disable_branching">true</parameter>
  </qmc>)";
  Libxml2Document configured_doc;
  REQUIRE(configured_doc.parseFromString(configured_input));
  REQUIRE(wc_.put(configured_doc.getRoot()));
  CHECK_FALSE(wc_.use_nonblocking_);
  CHECK(wc_.debug_disable_branching_);

  const char* const invalid_input = R"(
  <qmc method="dmc_batch">
    <parameter name="use_nonblocking">invalid</parameter>
  </qmc>)";
  Libxml2Document invalid_doc;
  REQUIRE(invalid_doc.parseFromString(invalid_input));
  CHECK_THROWS_MATCHES(wc_.put(invalid_doc.getRoot()), UniformCommunicateError,
                       Catch::Matchers::Message("WalkerControl::put parsing error. use_nonblocking only accepts 'yes'/'no'/'true'/'false' but the input value is 'invalid'."));
}

} // namespace testing

TEST_CASE("WalkerControl input", "[drivers][walker_control][input]")
{
  testing::UnifiedDriverWalkerControlMPITest test;
  test.testInput();
}

std::vector<int> testing::UnifiedDriverWalkerControlMPITest::comb(
    const std::vector<QMCTraits::FullPrecRealType>& weights,
    const std::vector<QMCTraits::FullPrecRealType>& rank_weights,
    int rank,
    int population,
    QMCTraits::FullPrecRealType offset,
    std::vector<int>& counts)
{
  UPtrVector<MCPopulation::MCPWalker> walkers;
  for (auto weight : weights)
  {
    walkers.push_back(std::make_unique<MCPopulation::MCPWalker>(1));
    walkers.back()->Weight = weight;
  }
  WalkerControl::assignStochasticComb(walkers, rank_weights, rank, population, offset, counts);
  std::vector<int> copies;
  for (const auto& walker : walkers)
    copies.push_back(static_cast<int>(walker->Multiplicity));
  return copies;
}

TEST_CASE("WalkerControl stochastic comb boundaries and offset sweep", "[drivers][stochastic_comb]")
{
  using Real    = QMCTraits::FullPrecRealType;
  using Fixture = testing::UnifiedDriverWalkerControlMPITest;
  std::vector<int> counts;
  CHECK(Fixture::comb({0.4, 1.6, 2.8, 3.2}, {8}, 0, 4, 0.25, counts) == std::vector<int>{0, 1, 2, 1});
  CHECK(counts == std::vector<int>{4});
  for (Real offset : {Real(0), Real(0.25), std::nextafter(Real(1), Real(0))})
  {
    CHECK(Fixture::comb({2, 2, 2, 2}, {8}, 0, 4, offset, counts) == std::vector<int>{1, 1, 1, 1});
    CHECK(Fixture::comb({0, 8, 0}, {8}, 0, 4, offset, counts) == std::vector<int>{0, 4, 0});
    CHECK(Fixture::comb({0, 0}, {0, 8}, 0, 4, offset, counts) == std::vector<int>{0, 0});
    CHECK(counts == std::vector<int>{0, 4});
    CHECK(Fixture::comb({3}, {3}, 0, 1, offset, counts) == std::vector<int>{1});
    std::vector<Real> rank_weights(16, 4);
    rank_weights.back() = 6;
    for (int rank = 0; rank < 16; ++rank)
    {
      const std::vector<Real> equal_weights(rank == 15 ? 3 : 2, 2);
      CHECK(Fixture::comb(equal_weights, rank_weights, rank, 33, offset, counts) ==
            std::vector<int>(equal_weights.size(), 1));
    }
  }
  CHECK(Fixture::comb({1.e-200, 1, 0}, {1}, 0, 4, 0, counts) == std::vector<int>{1, 3, 0});
  CHECK(Fixture::comb({1.e-200, 1, 0}, {1}, 0, 4, 0.25, counts) == std::vector<int>{0, 4, 0});
  CHECK(Fixture::comb({std::numeric_limits<Real>::min(), 1.e300}, {1.e300}, 0, 4, 0, counts) == std::vector<int>{1, 3});
  for (Real scale : {Real(1.e-300), Real(1), Real(1.e308)})
  {
    CHECK(Fixture::comb({0.1 * scale, 0.9 * scale}, {scale}, 0, 1, 0.05, counts) == std::vector<int>{1, 0});
    CHECK(Fixture::comb({0.1 * scale, 0.9 * scale}, {scale}, 0, 1, 0.5, counts) == std::vector<int>{0, 1});
  }

  const std::vector<Real> weights{0, 0.4, 1.6, 2.8, 3.2, 0};
  std::vector<Real> means(weights.size(), 0);
  constexpr int samples = 1024;
  for (int sample = 0; sample < samples; ++sample)
  {
    const auto copies = Fixture::comb(weights, {8}, 0, 4, Real(sample) / samples, counts);
    CHECK(std::accumulate(copies.begin(), copies.end(), 0) == 4);
    for (int iw = 0; iw < copies.size(); ++iw)
    {
      const Real expected = weights[iw] / 2;
      CHECK(copies[iw] >= std::floor(expected));
      CHECK(copies[iw] <= std::ceil(expected));
      means[iw] += Real(copies[iw]) / samples;
    }
  }
  for (int iw = 0; iw < weights.size(); ++iw)
    CHECK(std::abs(means[iw] - weights[iw] / 2) <= Real(1) / samples);

  for (const auto& rank_weights : std::vector<std::vector<Real>>{{0},
                                                                 {-1},
                                                                 {std::numeric_limits<Real>::infinity()},
                                                                 {std::numeric_limits<Real>::quiet_NaN()}})
    CHECK_THROWS_AS(Fixture::comb({1}, rank_weights, 0, 1, 0.25, counts), UniformCommunicateError);
  CHECK_THROWS_AS(Fixture::comb({1}, {1}, 0, 0, 0.25, counts), UniformCommunicateError);
  CHECK_THROWS_AS(Fixture::comb({1}, {1}, 0, 1, 1, counts), UniformCommunicateError);
}

void testing::UnifiedDriverWalkerControlMPITest::testFixedPopulation(int scenario, bool blocking)
{
  const int rank       = dpools_.comm.rank();
  const int ranks      = dpools_.comm.size();
  const int population = 2 * ranks + 1;
  std::vector<int> offsets(ranks + 1);
  FairDivideLow(population, ranks, offsets);
  pop_->createWalkers(offsets[rank + 1] - offsets[rank], walker_confs);
  pop_->syncWalkersPerRank(dpools_.comm);
  wc_.use_fixed_pop_   = true;
  wc_.use_nonblocking_ = !blocking;
  std::vector<double> global_weights(population);
  for (int iw = 0; iw < population; ++iw)
    global_weights[iw] = scenario == 0 ? 2 : scenario == 1 ? (iw == 0 ? 2 * population : 0) : (iw % 4) * 0.75;
  for (int iw = 0; iw < pop_->get_num_local_walkers(); ++iw)
  {
    auto& walker                                              = *pop_->get_walkers()[iw];
    walker.Weight                                             = global_weights[offsets[rank] + iw];
    walker.R[0][0]                                            = offsets[rank] + iw + 1;
    walker.Properties(WalkerProperties::Indexes::LOCALENERGY) = offsets[rank] + iw;
  }

  for (int step = 0; step < 3; ++step)
  {
    std::vector<double> histogram(population, 0);
    double local_weight = 0, local_energy = 0, local_energy2 = 0;
    int actual_movement  = 0;
    int positive_weights = 0;
    for (const auto& walker : pop_->get_walkers())
    {
      local_weight += walker->Weight;
      positive_weights += walker->Weight > 0;
      const double energy = walker->Properties(WalkerProperties::Indexes::LOCALENERGY);
      local_energy += walker->Weight * energy;
      local_energy2 += walker->Weight * energy * energy;
    }
    dpools_.comm.allreduce(local_weight);
    dpools_.comm.allreduce(local_energy);
    dpools_.comm.allreduce(local_energy2);
    dpools_.comm.allreduce(positive_weights);
    auto expected_rng = Random.makeClone();
    double offset     = (*expected_rng)();
    std::vector<RandomBase<double>::uint_type> expected_state, actual_state;
    if (rank == 0)
      expected_rng->save(expected_state);
    else
      Random.save(expected_state);
    dpools_.comm.bcast(offset);

    // Independent serial oracle uses explicit teeth; production never constructs this array.
    std::vector<int> expected(population, 0);
    long double cumulative = global_weights[0];
    int parent             = 0;
    for (int tooth = 0; tooth < population; ++tooth)
    {
      const long double position = (static_cast<long double>(tooth) + offset) * local_weight / population;
      while (parent + 1 < population && position >= cumulative)
        cumulative += global_weights[++parent];
      ++expected[parent];
    }
    wc_.branch(step, *pop_, false);
    Random.save(actual_state);
    CHECK(actual_state == expected_state);
    CHECK(pop_->get_num_local_walkers() == offsets[rank + 1] - offsets[rank]);
    CHECK(pop_->get_num_global_walkers() == population);
    CHECK(wc_.ensemble_property_.Weight == Approx(local_weight));
    CHECK(wc_.ensemble_property_.NumSamples == population);
    CHECK(wc_.ensemble_property_.LivingFraction == Approx(double(positive_weights) / population));
    CHECK(wc_.ensemble_property_.Energy == Approx(local_energy / local_weight));
    CHECK(wc_.ensemble_property_.Variance ==
          Approx(local_energy2 / local_weight - std::pow(local_energy / local_weight, 2)).margin(1.e-10));
    int logical_movement = 0;
    for (int ir = 0; ir < ranks; ++ir)
      logical_movement += std::max(0, wc_.num_per_rank_[ir] - offsets[ir + 1] + offsets[ir]);
    int payloads = wc_.saved_num_walkers_sent_;
    dpools_.comm.allreduce(payloads);
    CHECK(payloads <= logical_movement);
    if (scenario == 0)
      CHECK(payloads == 0);
    if (scenario == 1 && step == 0)
    {
      CHECK(logical_movement == population - (offsets[1] - offsets[0]));
      CHECK(payloads == ranks - 1);
    }
    for (const auto& walker : pop_->get_walkers())
    {
      CHECK(walker->Weight == Approx(local_weight / population).epsilon(1.e-12));
      CHECK(walker->Multiplicity == 1);
      CHECK(walker->getWalkerID() > 0);
      const int original = static_cast<int>(walker->R[0][0]) - 1;
      REQUIRE(original >= 0);
      REQUIRE(original < population);
      histogram[original] += 1;
      if (original < offsets[rank] || original >= offsets[rank + 1])
      {
        ++actual_movement;
        CHECK(walker->wasTouched);
      }
    }
    dpools_.comm.allreduce(actual_movement);
    CHECK(actual_movement == logical_movement);
    dpools_.comm.allreduce(histogram);
    for (int iw = 0; iw < population; ++iw)
      CHECK(histogram[iw] == expected[iw]);

    // Subsequent events carry the preceding non-unit W/N through multiplicative propagation.
    for (auto& walker : pop_->get_walkers())
      walker->Weight *= 1.25;
    // The rank-major parent ordering changes after balancing; gather the next oracle's marked weights.
    global_weights.assign(population, 0);
    for (int iw = 0; iw < pop_->get_num_local_walkers(); ++iw)
    {
      auto& walker                       = *pop_->get_walkers()[iw];
      walker.R[0][0]                     = offsets[rank] + iw + 1;
      global_weights[offsets[rank] + iw] = walker.Weight;
    }
    dpools_.comm.allreduce(global_weights);
  }
}

void testing::UnifiedDriverWalkerControlMPITest::testSkippedComb(bool debug)
{
  const int rank = dpools_.comm.rank();
  pop_->createWalkers(rank == 0 ? 3 : 1, walker_confs);
  pop_->syncWalkersPerRank(dpools_.comm);
  wc_.use_fixed_pop_           = true;
  wc_.debug_disable_branching_ = debug;
  std::vector<double> weights;
  for (int iw = 0; iw < pop_->get_num_local_walkers(); ++iw)
  {
    pop_->get_walkers()[iw]->Weight       = 0.5 + iw + rank;
    pop_->get_walkers()[iw]->Multiplicity = 7;
    weights.push_back(pop_->get_walkers()[iw]->Weight);
  }
  std::vector<RandomBase<double>::uint_type> before, after;
  Random.save(before);
  wc_.branch(0, *pop_, !debug);
  Random.save(after);
  CHECK(before == after);
  CHECK(pop_->get_num_local_walkers() == weights.size());
  for (int iw = 0; iw < weights.size(); ++iw)
  {
    CHECK(pop_->get_walkers()[iw]->Weight == weights[iw]);
    CHECK(pop_->get_walkers()[iw]->Multiplicity == 7);
  }
  CHECK(wc_.saved_num_walkers_sent_ == 0);
  CHECK(wc_.ensemble_property_.LivingFraction == 1);
}

void testing::UnifiedDriverWalkerControlMPITest::testInvalidWeights(QMCTraits::FullPrecRealType invalid_weight)
{
  wc_.use_fixed_pop_ = true;
  pop_->syncWalkersPerRank(dpools_.comm);
  pop_->get_walkers()[0]->Weight =
      invalid_weight == 0 || dpools_.comm.rank() == dpools_.comm.size() - 1 ? invalid_weight : 1;
  CHECK_THROWS_AS(wc_.branch(0, *pop_, false), UniformCommunicateError);
  CHECK(pop_->get_num_local_walkers() == 1);
  CHECK(pop_->get_walkers()[0]->Multiplicity == 1);
}

TEST_CASE("WalkerControl fixed population MPI and copying", "[drivers][stochastic_comb]")
{
  for (int scenario : {0, 1, 2})
    for (bool blocking : {false, true})
    {
      testing::UnifiedDriverWalkerControlMPITest test;
      test.testFixedPopulation(scenario, blocking);
    }
}

TEST_CASE("WalkerControl skipped comb and coordinated errors", "[drivers][stochastic_comb]")
{
  for (bool debug : {false, true})
  {
    testing::UnifiedDriverWalkerControlMPITest test;
    test.testSkippedComb(debug);
  }
  for (double invalid : {0., -1., std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
  {
    testing::UnifiedDriverWalkerControlMPITest test;
    test.testInvalidWeights(invalid);
  }
}

TEST_CASE("WalkerControl::determineNewWalkerPopulation", "[drivers][walker_control]")
{
  // These must be empty vectors!
  std::vector<int> minus;
  std::vector<int> plus;
  testing::UnifiedDriverWalkerControlMPITest test;
  int num_ranks = test.getNumRanks();
  std::vector num_per_rank(num_ranks, 1);
  num_per_rank[0] = num_ranks;
  test.testNewDistribution(num_per_rank, minus, plus);
  int rank = test.getRank();
  app_log() << "rank:" << rank << " minus: " << NativePrint(minus) << std::endl;

  app_log() << "rank:" << rank << " plus: " << NativePrint(plus) << std::endl;
  CHECK(minus.size() == num_ranks - 1);
  CHECK(plus.size() == num_ranks - 1);
  app_log() << "rank:" << rank << " plus: " << NativePrint(num_per_rank) << std::endl;
}

void testing::UnifiedDriverWalkerControlMPITest::testPopulationDiff(std::vector<int>& rank_counts_before,
                                                                    std::vector<int>& rank_counts_after)
{
  using MCPWalker = MCPopulation::MCPWalker;

  int rank               = dpools_.comm.rank();
  auto num_local_walkers = pop_->get_num_local_walkers();

  // if rank counts_before is different manipulation multiplicities on walkers to get the rank count we want.
  if (num_local_walkers != rank_counts_before[rank])
  {
    auto walker_diff       = rank_counts_before[rank] - num_local_walkers;
    MCPWalker& last_walker = *(pop_->get_walkers()[num_local_walkers - 1]);
    last_walker.Multiplicity += walker_diff;
  }

  std::vector<int> fair_offset;
  std::vector<int> minus;
  std::vector<int> plus;

  wc_.setNumPerRank(rank_counts_before);
  // note this ordering of calls is from the main code
  wc_.swapWalkersSimple(*pop_);
  wc_.killDeadWalkersOnRank(*pop_);
  pop_->fissionHighMultiplicityWalkers();
  reportWalkersPerRank(dpools_.comm, *pop_);
  CHECK(pop_->get_num_local_walkers() == rank_counts_after[rank]);
}

void testing::UnifiedDriverWalkerControlMPITest::testWalkerIDs(std::vector<std::vector<int>> walker_ids_after,
                                                               std::vector<std::vector<int>> parent_ids_after)
{
  int rank = dpools_.comm.rank();
  std::vector<int> walker_ids;
  std::vector<int> parent_ids;
#ifndef NDEBUG
  for (int iw = 0; iw < pop_->get_walkers().size(); ++iw)
  {
    walker_ids.push_back(pop_->get_walkers()[iw]->getWalkerID());
    parent_ids.push_back(pop_->get_walkers()[iw]->getParentID());
  }
  app_log() << "rank: " << rank << "  walker ids: " << NativePrint(walker_ids)
            << " parent ids: " << NativePrint(parent_ids) << std::endl;
#endif
  for (int iw = 0; iw < walker_ids_after[rank].size(); ++iw)
  {
    CHECK(pop_->get_walkers()[iw]->getWalkerID() == walker_ids_after[rank][iw]);
  }
  for (int iw = 0; iw < parent_ids_after[rank].size(); ++iw)
  {
    CHECK(pop_->get_walkers()[iw]->getParentID() == parent_ids_after[rank][iw]);
  }
}

#define CALC_ID_LAMBDA \
  [num_ranks](int rank, int index_walker_created) { return index_walker_created * num_ranks + rank + 1; }

TEST_CASE("MPI WalkerControl population swap walkers", "[drivers][walker_control]")
{
  auto test_func = []() {
    outputManager.pause();
    testing::UnifiedDriverWalkerControlMPITest test;
    outputManager.resume();

    SECTION("Balanced")
    {
      std::vector<int> count_before(test.getNumRanks(), 1);
      std::vector<int> count_after(test.getNumRanks(), 1);
      // One walker on every node, should be no swapping
      test.testPopulationDiff(count_before, count_after);
    }

    SECTION("LoadBalance")
    {
      auto num_ranks = test.getNumRanks();
      std::vector<int> count_before(num_ranks, 1);
      count_before[0] = num_ranks;
      std::vector<int> count_after(num_ranks, 2);
      count_after[0] = 1;
      app_log() << "count_before: " << NativePrint(count_before) << "  count_after: " << NativePrint(count_after)
                << std::endl;
      test.testPopulationDiff(count_before, count_after);
      std::vector<std::vector<int>> ar_wids;
      std::vector<std::vector<int>> ar_pids;
      auto calcID = CALC_ID_LAMBDA;
      for (int ir = 0; ir < num_ranks; ++ir)
      {
        ar_wids.push_back({ir + 1});
        if (ir > 0)
          ar_wids.back().push_back(calcID(ir, 1));

        ar_pids.push_back({0});
        if (ir > 0)
        {
          ar_pids.back().push_back(1);
        }
      }
      test.testWalkerIDs(ar_wids, ar_pids);
    }

    SECTION("LoadBalance Multiple Copy Optimization")
    {
      int multiple   = 3;
      auto num_ranks = test.getNumRanks();
      std::vector<int> walker_multiplicity_total(num_ranks, 1);
      walker_multiplicity_total[0] = num_ranks * multiple;
      int total_walkers = std::accumulate(walker_multiplicity_total.begin(), walker_multiplicity_total.end(), 0);
      std::vector<int> count_after = fairDivide(total_walkers, num_ranks);
      std::reverse(count_after.begin(), count_after.end());
      app_log() << "walker_multiplicity_before: " << NativePrint(walker_multiplicity_total)
                << "count_after: " << NativePrint(count_after) << std::endl;
      test.testPopulationDiff(walker_multiplicity_total, count_after);
      std::vector<std::vector<int>> ar_wids;
      std::vector<std::vector<int>> ar_pids;
      auto calcID = CALC_ID_LAMBDA;

      for (int ir = 0; ir < num_ranks; ++ir)
      {
        ar_wids.push_back({calcID(ir, 0)});
        for (int iw = 1; iw < count_after[ir]; ++iw)
          ar_wids.back().push_back(calcID(ir, iw));

        ar_pids.push_back({0});
        if (ir == 0)
        {
          for (int iw = 1; iw < count_after[ir]; ++iw)
            ar_pids.back().push_back(calcID(0, 0));
        }
        else
        {
          ar_pids.back().push_back(1);
          for (int iw = 2; iw < count_after[ir]; ++iw)
            ar_pids.back().push_back(calcID(ir, 1));
        }
      }

      app_log() << std::endl;
      test.testWalkerIDs(ar_wids, ar_pids);
    }

    SECTION("MultiStage LoadBalance")
    {
      auto num_ranks = test.getNumRanks();
      std::vector<int> walker_multiplicity_total(num_ranks, 1);
      walker_multiplicity_total[0] = num_ranks;
      std::vector<int> count_after(num_ranks, 2);
      count_after[0] = 1;
      test.testPopulationDiff(walker_multiplicity_total, count_after);
      std::vector<std::vector<int>> ar_wids;
      std::vector<std::vector<int>> ar_pids;
      auto calcID = CALC_ID_LAMBDA;
      for (int ir = 0; ir < num_ranks; ++ir)
      {
        ar_wids.push_back({ir + 1});
        if (ir > 0)
          ar_wids.back().push_back(calcID(ir, 1));

        ar_pids.push_back({0});
        if (ir > 0)
        {
          ar_pids.back().push_back(1);
        }
      }
      test.testWalkerIDs(ar_wids, ar_pids);
      if (num_ranks == 1)
      {
        walker_multiplicity_total = count_after;
      }
      else
      {
        // Next Stage
        walker_multiplicity_total = count_after;
        walker_multiplicity_total[num_ranks - 1] += 2;
        count_after[0] += 1;
        count_after.back() += 1;
        app_log() << "walker_multiplicity_total: " << NativePrint(walker_multiplicity_total)
                  << "  count_after: " << NativePrint(count_after) << std::endl;
        test.testPopulationDiff(walker_multiplicity_total, count_after);
        ar_wids.back().push_back(calcID(num_ranks - 1, 2));
        ar_wids[0].push_back(calcID(0, 1));
        ar_pids.back().push_back(calcID(num_ranks - 1, 1));
        ar_pids[0].push_back(calcID(num_ranks - 1, 1));
        test.testWalkerIDs(ar_wids, ar_pids);

        // Walker dropped
        walker_multiplicity_total = count_after;
        walker_multiplicity_total[num_ranks - 1] -= 1;
        count_after[num_ranks - 1] -= 1;
        test.testPopulationDiff(walker_multiplicity_total, count_after);
        ar_wids.back().pop_back();
        ar_pids.back().pop_back();
        app_log() << "walker_multiplicity_total: " << NativePrint(walker_multiplicity_total)
                  << "  count_after: " << NativePrint(count_after) << std::endl;
        test.testWalkerIDs(ar_wids, ar_pids);

        // Walkers added again
        walker_multiplicity_total = count_after;
        walker_multiplicity_total[num_ranks - 1] += 2;
        count_after[num_ranks - 2] += 1;
        count_after.back() += 1;
        app_log() << "walker_multiplicity_total: " << NativePrint(walker_multiplicity_total)
                  << "  count_after: " << NativePrint(count_after) << std::endl;
        test.testPopulationDiff(walker_multiplicity_total, count_after);
        ar_wids.back().push_back(calcID(num_ranks - 1, 3));
        ar_wids[num_ranks - 2].push_back(calcID(num_ranks - 2, 2));
        ar_pids[num_ranks - 1].push_back(calcID(num_ranks - 1, 1));
        ar_pids[num_ranks - 2].push_back(calcID(num_ranks - 1, 1));
        test.testWalkerIDs(ar_wids, ar_pids);
      }
    }
  };
  MPIExceptionWrapper mew;
  mew(test_func);
}

} // namespace qmcplusplus
