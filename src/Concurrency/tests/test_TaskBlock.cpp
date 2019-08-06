//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include <atomic>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "Concurrency/TaskBlock.hpp"

namespace qmcplusplus
{

class IncrementOnDestruction
{
  std::atomic<int>& counter_;
public:
  IncrementOnDestruction(std::atomic<int>& counter) : counter_(counter) {}
  ~IncrementOnDestruction() { counter_++; }
};

template<Threading TT>
void TestTaskWithBarrier(const int ip,
			 TaskBlockBarrier<TT>& barrier,
			 int seconds,
			 std::atomic<int>& counter)
{
  IncrementOnDestruction inc(counter);
  sleep(seconds + 2 * ip);
  REQUIRE(counter == 0);
  barrier.wait();
}

template<Threading TT>
void TestTask(const int ip,
	      std::atomic<int>& counter)
{
  counter++;
}



TEST_CASE("TaskBlock simplecase openmp", "[concurrency]")
{
  int num_threads = 8;
  TaskBlock<Threading::OPENMP> test_block(num_threads);
  std::atomic<int> count(0);
  test_block(TestTask<Threading::OPENMP>,
	     std::ref(count));
  REQUIRE(count == 8);
}

TEST_CASE("TaskBlock barrier case openmp", "[concurrency]")
{
  int num_threads = 2;
  TaskBlock<Threading::OPENMP> test_block(num_threads);
  TaskBlockBarrier<Threading::OPENMP> barrier(num_threads);
  std::atomic<int> count(0);
  test_block(TestTaskWithBarrier<Threading::OPENMP>,
	     barrier,
	     1,
	     std::ref(count));
  REQUIRE(count == 2);
}

TEST_CASE("TaskBlock simple case std::thread", "[concurrency]")
{
  int num_threads = 8;
  TaskBlock<Threading::STD> test_block(num_threads);
  std::atomic<int> count(0);
  test_block(TestTask<Threading::STD>,
	     std::ref(count));
  REQUIRE(count == 8);
}

TEST_CASE("TaskBlock barrier case std::thread", "[concurrency]")
{
  int num_threads = 2;
  TaskBlock<Threading::STD> test_block(num_threads);
  TaskBlockBarrier<Threading::STD> barrier(num_threads);
  std::atomic<int> count(0);
  test_block(TestTaskWithBarrier<Threading::STD>,
	     std::ref(barrier),
	     1,
	     std::ref(count));
  REQUIRE(count == 2);
}

}
