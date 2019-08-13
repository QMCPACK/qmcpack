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

#include "Concurrency/TasksOneToOne.hpp"

namespace qmcplusplus
{

void TestTask(const int ip,
              std::atomic<int>& counter)
{
  counter++;
}

/** Openmp generally works but is not guaranteed with std::atomic
 */
void TestTaskOMP(const int ip,
           	 int& counter)
{
  #pragma omp atomic update
  counter++;
}


TEST_CASE("TaskBlock function case openmp", "[concurrency]")
{
  int num_threads = 8;
  TasksOneToOne<Threading::OPENMP> test_block(num_threads);
  int count(0);
  test_block(TestTaskOMP,
	     std::ref(count));
  REQUIRE(count == 8);
}

TEST_CASE("TaskBlock lambda case openmp", "[concurrency]")
{
  int num_threads = 8;
  TasksOneToOne<Threading::OPENMP> test_block(num_threads);
  int count(0);
  test_block([](int id, int& c){
#pragma omp atomic update
		 c++;
	     },count);
  REQUIRE(count == 8);
}

TEST_CASE("TaskBlock simple case std::thread", "[concurrency]")
{
  int num_threads = 8;
  TasksOneToOne<Threading::STD> test_block(num_threads);
  std::atomic<int> count(0);
  test_block(TestTask,
	     std::ref(count));
  REQUIRE(count == 8);
}

}
