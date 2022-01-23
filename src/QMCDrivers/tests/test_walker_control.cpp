//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "Utilities/RandomGenerator.h"
#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCDrivers/WalkerControlBase.h"
#ifdef HAVE_MPI
#include "QMCDrivers/DMC/WalkerControlMPI.h"
#include "QMCDrivers/DMC/WalkerReconfigurationMPI.h"
#endif

#include <stdio.h>
#include <string>
#include <random>

using std::string;

namespace qmcplusplus
{
void output_vector(const std::string& name, std::vector<int>& vec)
{
  std::cout << name;
  for (int i = 0; i < vec.size(); i++)
  {
    std::cout << vec[i] << " ";
  }
  std::cout << std::endl;
}

// uncomment the std::cout and output_vector lines to see the walker assignments
TEST_CASE("Walker control assign walkers", "[drivers][walker_control]")
{
  int Cur_pop     = 8;
  int NumContexts = 4;
  std::vector<int> FairOffset(NumContexts + 1);

  // The in loop copy is necessary to support the assert at the end of the swaps.
  // This was important for debugging but will go in a future PR as part of cleaning
  // update determineNewWalkerPopulation
  std::vector<int> NumPerRank = {4, 4, 0, 0};
  std::vector<int> NewNum     = NumPerRank;
  for (int me = 0; me < NumContexts; me++)
  {
    std::vector<int> minus;
    std::vector<int> plus;
    std::vector<int> num_per_rank = NumPerRank;

    //std::cout << "For processor number " << me << std::endl;
    WalkerControlMPI::determineNewWalkerPopulation(Cur_pop, NumContexts, me, num_per_rank, FairOffset, minus, plus);

    REQUIRE(minus.size() == plus.size());
    //output_vector("  Minus: ", minus);

    //output_vector("  Plus: ", plus);

    for (int i = 0; i < plus.size(); i++)
    {
      if (me == plus[i])
        NewNum[plus[i]]--;
    }
    for (int i = 0; i < minus.size(); i++)
    {
      if (me == minus[i])
        NewNum[minus[i]]++;
    }
  }
  //output_vector("New num per node: ", NewNum);

  for (int i = 0; i < NewNum.size(); i++)
  {
    int num = FairOffset[i + 1] - FairOffset[i];
    REQUIRE(NewNum[i] == num);
  }
}

TEST_CASE("WalkerControl round trip index conversions", "[drivers][walker_control]")
{
  std::vector<int> walker_counts; //= {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 23, 37};
  for (int i = 0; i < 1000; ++i)
    walker_counts.push_back(i);
  for (int i = 0; i < 1000000; i += 1000)
    walker_counts.push_back(i);
  for (int i = 0; i < 100000000; i += 100000)
    walker_counts.push_back(i);

  std::vector<QMCTraits::FullPrecRealType> fp_counts(walker_counts.size());
  std::vector<int> walker_count_results(walker_counts.size());
  for (int iw = 0; iw < walker_counts.size(); ++iw)
  {
    fp_counts[iw] = walker_counts[iw];
  }
  for (int iw = 0; iw < fp_counts.size(); ++iw)
  {
    walker_count_results[iw] = static_cast<int>(fp_counts[iw]);
  }
  bool all_pass = true;
  for (int iw = 0; iw < walker_counts.size(); ++iw)
  {
    all_pass &= (walker_counts[iw] == walker_count_results[iw]);
  }
  REQUIRE(all_pass);
}

// uncomment the std::cout and output_vector lines to see the walker assignments
TEST_CASE("Walker control assign walkers odd ranks", "[drivers][walker_control]")
{
  int Cur_pop                 = 9;
  int NumContexts             = 3;
  std::vector<int> NumPerRank = {5, 2, 2};
  std::vector<int> FairOffset(NumContexts + 1);

  std::vector<int> NewNum = NumPerRank;
  for (int me = 0; me < NumContexts; me++)
  {
    std::vector<int> minus;
    std::vector<int> plus;
    std::vector<int> num_per_rank = NumPerRank;
    //std::cout << "For processor number " << me << std::endl;
    WalkerControlMPI::determineNewWalkerPopulation(Cur_pop, NumContexts, me, num_per_rank, FairOffset, minus, plus);

    REQUIRE(minus.size() == plus.size());
    //output_vector("  Minus: ", minus);

    //output_vector("  Plus: ", plus);

    for (int i = 0; i < plus.size(); i++)
    {
      if (me == plus[i])
        NewNum[plus[i]]--;
    }
    for (int i = 0; i < minus.size(); i++)
    {
      if (me == minus[i])
        NewNum[minus[i]]++;
    }
  }
  //output_vector("New num per node: ", NewNum);

  for (int i = 0; i < NewNum.size(); i++)
  {
    int num = FairOffset[i + 1] - FairOffset[i];
    REQUIRE(NewNum[i] == num);
  }
}

#ifdef PROPERTY_TESTING
// Eventually will create some way to build and run property-based tests, which
// are tests that use random inputs and verify that certain properties hold true.
//
// In this case, the test constructs random number of processors with a random
// number of walkers on each node. The test checks two properties - that the
// plus and minus lists have the same size, and that each node ultimately has
// the number of walkers chosen by FairDivideLow.

TEST_CASE("Walker control assign walkers many", "[drivers][walker_control][property]")
{
  // Use random device for seed for coverage
  //std::random_device rd;
  //std::mt19937 mt(rd());
  // Use fixed seed for reproducibility
  std::mt19937 mt(100);
  std::uniform_int_distribution<int> NumNodes(1, 1000);
  std::uniform_int_distribution<int> TotalPop(0, 1000);

  for (int nt = 0; nt < 1000; nt++)
  {
    int NumContexts = NumNodes(mt);
    int Cur_pop     = TotalPop(mt);
    std::uniform_int_distribution<int> WalkerPop(0, 2 * Cur_pop / NumContexts);
    std::vector<int> NumPerRank(NumContexts);
    int current_pop = Cur_pop;
    for (int i = 0; i < NumContexts; i++)
    {
      int p = WalkerPop(mt);
      p     = std::min(current_pop, p);
      current_pop -= p;
      // Make sure all walkers are accounted for on the last node
      if (i == NumContexts - 1 && current_pop > 0)
      {
        p += current_pop;
      }
      NumPerRank[i] = p;
    }
    //std::cout << "NumNodes = " << NumContexts << std::endl;
    //std::cout << "TotalPop = " << Cur_pop << std::endl;
    //output_vector("Start: ",NumPerRank);

    std::vector<int> NewNum = NumPerRank;
    for (int me = 0; me < NumContexts; me++)
    {
      std::vector<int> minus;
      std::vector<int> plus;

      determineNewWalkerPopulation(Cur_pop, NumContexts, me, NumPerRank, minus, plus);
      REQUIRE(minus.size() == plus.size());

      for (int i = 0; i < plus.size(); i++)
      {
        if (me == plus[i])
          NewNum[plus[i]]--;
      }
      for (int i = 0; i < minus.size(); i++)
      {
        if (me == minus[i])
          NewNum[minus[i]]++;
      }
    }

    std::vector<int> FairOffset;
    FairDivideLow(Cur_pop, NumContexts, FairOffset);
    for (int i = 0; i < NewNum.size(); i++)
    {
      int num = FairOffset[i + 1] - FairOffset[i];
      REQUIRE(NewNum[i] == num);
    }
  }
}
#endif

#ifdef HAVE_MPI

struct WalkerControlMPITest
{
  /** Currently only passes for 1,2, or 3 ranks
   */
  void operator()(bool use_nonblocking)
  {
    Communicate* c = OHMMS::Controller;
    WalkerControlMPI wc(c);

    wc.use_nonblocking = use_nonblocking;

    const SimulationCell simulation_cell;
    MCWalkerConfiguration W(simulation_cell); // Unused in the function
    using Walker_t = MCWalkerConfiguration::Walker_t;

    //UPtrVector<Walker_t> walkers;
    // Set up Cur_pop
    // Set up good_w and bad_w

    wc.Cur_pop = c->size();
    for (int i = 0; i < c->size(); i++)
    {
      wc.NumPerRank[i] = 1;
    }
    // One walker on every node, should be no swapping
    wc.good_w.push_back(std::make_unique<Walker_t>());
    wc.good_w[0]->ID = c->rank();
    wc.ncopy_w.push_back(0);

    wc.swapWalkersSimple(W);

    REQUIRE(wc.good_w.size() == 1);
    REQUIRE(wc.bad_w.size() == 0);


    // 3 walkers on rank 0, 1 walker on others - should redistribute if
    //   there is more than one rank
    if (c->size() > 1)
    {
      if (c->rank() == 0)
      {
        wc.good_w.push_back(std::make_unique<Walker_t>());
        wc.good_w.push_back(std::make_unique<Walker_t>());

        // Use the ID variable to check that the walker content was transmitted
        wc.good_w[1]->ID = c->size();
        wc.good_w[2]->ID = c->size() + 1;

        wc.ncopy_w.push_back(0);
        wc.ncopy_w.push_back(0);
      }
      wc.NumPerRank[0] = 3;
      wc.Cur_pop += 2;

      wc.swapWalkersSimple(W);

      //std::cout << " Rank = " << c->rank() << " good size = " << wc.good_w.size() <<
      //          " ID = " << wc.good_w[0]->ID << std::endl;

      if (c->rank() == c->size() - 2)
      {
        REQUIRE(wc.good_w.size() == 2);
        // This check is a bit too restrictive - no guarantee the last walker was the
        //  one transmitted
        bool okay1 = wc.good_w[1]->ID == c->size() || wc.good_w[1]->ID == c->size() + 1;
        REQUIRE(okay1);
      }
      else if (c->rank() == c->size() - 1)
      {
        REQUIRE(wc.good_w.size() == 2);
        bool okay2 = wc.good_w[1]->ID == c->size() || wc.good_w[1]->ID == c->size() + 1;
        REQUIRE(okay2);
      }
      else
      {
        REQUIRE(wc.good_w.size() == 1);
        REQUIRE(wc.good_w[0]->ID == c->rank());
      }
      wc.NumPerRank[0]             = 1;
      wc.NumPerRank[c->size() - 1] = 2;
      wc.NumPerRank[c->size() - 2] = 2;
    }


    // And now the strange case
    // 6 walkers on rank0, 2 on rank1, 2 on rank2
    if (c->size() > 2)
    {
      if (c->rank() == 0)
      {
        wc.good_w.push_back(std::make_unique<Walker_t>());
        wc.good_w.push_back(std::make_unique<Walker_t>());
        // wc.good_w.push_back(new Walker_t());
        // wc.good_w.push_back(new Walker_t());
        int nwalkers_rank                = wc.good_w.size();
        wc.good_w[nwalkers_rank - 1]->ID = c->size() + 5;
        wc.good_w[nwalkers_rank - 2]->ID = c->size() + 4;
        // wc.good_w[nwalkers_rank - 3]->ID = c->size() + 3;
        // wc.good_w[nwalkers_rank - 4]->ID = c->size() + 2;

        wc.ncopy_w.push_back(2);
        wc.ncopy_w.push_back(1);
        // wc.ncopy_w.push_back(0);
        // wc.ncopy_w.push_back(0);
      }
      else if (c->rank() == 1)
      {
        //wc.bad_w.push_back(wc.good_w[0]);
        //wc.bad_w.push_back(wc.good_w[1]);
        int nwalkers_rank = wc.good_w.size();
        //wc.good_w.pop_back();
        //wc.good_w.pop_back();
        //wc.ncopy_w.pop_back();
        //wc.ncopy_w.pop_back();
      }
      wc.NumPerRank[0] = 6;
      wc.Cur_pop += 5;

      reportWalkersPerRank(c, wc);

      wc.swapWalkersSimple(W);

      reportWalkersPerRank(c, wc);


      // These are unique walkers
      if (c->rank() == c->size() - 2)
      {
        CHECK(wc.good_w.size() == 3);
      }
      else if (c->rank() == c->size() - 1)
      {
        CHECK(wc.good_w.size() == 3);
      }
      else
      {
        CHECK(wc.good_w.size() == 2);
      }

      int walker_count = wc.copyWalkers(W);

      reportWalkersPerRank(c, wc);

      if (c->rank() == c->size() - 2)
      {
        CHECK(walker_count == 3);
      }
      else if (c->rank() == c->size() - 1)
      {
        CHECK(walker_count == 4);
      }
      else
      {
        CHECK(walker_count == 3);
      }
    }
  }

private:
  void reportWalkersPerRank(Communicate* c, WalkerControlMPI& wc)
  {
    std::vector<int> rank_walker_count(c->size(), 0);
    rank_walker_count[c->rank()] = wc.good_w.size();
    c->allreduce(rank_walker_count);
    if (c->rank() == 0)
    {
      std::cout << "Walkers Per Rank (Total: " << wc.Cur_pop << ")\n";
      for (int i = 0; i < rank_walker_count.size(); ++i)
      {
        std::cout << " " << i << "  " << rank_walker_count[i] << "\n";
      }
    }
  }
};

void test_swap_walkers(bool use_nonblocking)
{
  WalkerControlMPITest test;
  test(use_nonblocking);
}

TEST_CASE("Walker control swap walkers blocking", "[drivers][walker_control]")
{
  SECTION("blocking")
  {
    bool non_blocking = false;
    test_swap_walkers(non_blocking);
  }
}

TEST_CASE("Walker control swap walkers nonblocking", "[drivers][walker_control]")
{
  SECTION("nonblocking")
  {
    bool non_blocking = true;
    test_swap_walkers(non_blocking);
  }
}

TEST_CASE("Walker control reconfiguration", "[drivers][walker_control]")
{
  Communicate* c = OHMMS::Controller;
  WalkerReconfigurationMPI wr(c);

  wr.dN.resize(c->size());

  wr.dN[c->rank()] = 0;
  if (c->size() > 1)
  {
    wr.dN[0] = 1;
    wr.dN[1] = -1;
  }

  const SimulationCell simulation_cell;
  MCWalkerConfiguration W(simulation_cell);
  W.createWalkers(1);
  using Walker_t = MCWalkerConfiguration::Walker_t;

  if (c->rank() == 0)
  {
    W[0]->ID = 100;
  }
  else
  {
    W[0]->ID = 1;
  }

  std::vector<int> plus, minus;
  if (c->size() > 1)
  {
    if (c->rank() == 0)
      plus.push_back(0);
    if (c->rank() == 1)
      minus.push_back(0);
  }

  if (plus.size() > 0)
  {
    wr.sendWalkers(W, plus);
    REQUIRE(W.WalkerList.size() == 1);
    REQUIRE(W[0]->ID == 100);
  }

  if (minus.size() > 0)
  {
    wr.recvWalkers(W, minus);
    REQUIRE(W.WalkerList.size() == 1);
    // Use ID and ParentID to check the walker was transmitted
    REQUIRE(W[0]->ID == 1 + c->size());
    REQUIRE(W[0]->ParentID == 100);
  }
}

#endif
} // namespace qmcplusplus
