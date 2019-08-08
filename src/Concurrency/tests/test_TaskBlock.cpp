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
#include <chrono>
#include <thread>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "Concurrency/TaskBlock.hpp"

namespace qmcplusplus
{

template<Threading TT>
class IncrementOnDestruction
{
  std::atomic<int>& counter_;
public:
  IncrementOnDestruction(std::atomic<int>& counter) : counter_(counter) {}
  ~IncrementOnDestruction() { counter_++; }
};

/** Openmp is not necessarily ok with std::atomic
 *  So cannot use the general test
 *
 *  Since in general I think std::atomic should be usable
 *  I leave this in as the general case for other models
 *  And as a reminder that openmp is really a C oriented runtime
 *  That doesn't support stdlibc++
 */
class IncrementOnDestructionOMP
{
  int& counter_;
public:
  IncrementOnDestructionOMP(int& counter) : counter_(counter) {}
  ~IncrementOnDestructionOMP()
  {
      #pragma omp atomic update
      counter_++;
  }
};

template<Threading TT>
void TestTaskWithBarrier(const int ip,
			 TaskBlockBarrier<TT>& barrier,
			 int milliseconds,
			 std::atomic<int>& counter)
{
  IncrementOnDestruction<TT> inc(counter);
  std::this_thread::sleep_for(std::chrono::milliseconds(milliseconds + 100 * ip));
  REQUIRE(counter == 0);
  barrier.wait();
}

void TestTaskWithBarrierOMP(const int ip,
			 TaskBlockBarrier<Threading::OPENMP>& barrier,
			 int milliseconds,
			 int& counter)
{
  IncrementOnDestructionOMP inc(counter);
  std::this_thread::sleep_for(std::chrono::milliseconds(milliseconds + 100 * ip));
  REQUIRE(counter == 0);
  barrier.wait();
}

template<Threading TT>
void TestTask(const int ip,
	      std::atomic<int>& counter)
{
  counter++;
}

void TestTaskOMP(const int ip,
				 int& counter)
{
  #pragma omp atomic update
  counter++;
}



TEST_CASE("TaskBlock simplecase openmp", "[concurrency]")
{
  int num_threads = 8;
  TaskBlock<Threading::OPENMP> test_block(num_threads);
  int count(0);
  test_block(TestTaskOMP,
	     count);
  REQUIRE(count == 8);
}

TEST_CASE("TaskBlock barrier case openmp", "[concurrency]")
{
  int num_threads = 2;
  TaskBlock<Threading::OPENMP> test_block(num_threads);
  TaskBlockBarrier<Threading::OPENMP> barrier(num_threads);
  int count(0);
  test_block(TestTaskWithBarrierOMP,
	     std::ref(barrier),
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
