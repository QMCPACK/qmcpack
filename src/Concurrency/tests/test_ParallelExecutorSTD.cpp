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
#include <thread>
#include <functional>

#include "catch.hpp"

#include "Concurrency/ParallelExecutor.hpp"

namespace qmcplusplus
{
void TestTask(const int ip, std::atomic<int>& counter) { ++counter; }

TEST_CASE("ParallelExecutor<STD> function case", "[concurrency]")
{
  int num_threads = 8;
  ParallelExecutor<Executor::STD_THREADS> test_block;
  std::atomic<int> count(0);
  test_block(num_threads, TestTask, std::ref(count));
  REQUIRE(count == 8);
}

TEST_CASE("ParallelExecutor<STD> lambda case", "[concurrency]")
{
  int num_threads = 8;
  ParallelExecutor<Executor::STD_THREADS> test_block;
  std::atomic<int> count(0);
  test_block(
      num_threads, [](int id, std::atomic<int>& my_count) { ++my_count; }, std::ref(count));
  REQUIRE(count == 8);
}

TEST_CASE("ParallelExecutor<STD> nested case", "[concurrency]")
{
  int num_threads = 8;
  ParallelExecutor<Executor::STD_THREADS> test_block;
  std::atomic<int> count(0);
  test_block(
      num_threads,
      [num_threads](int task_id, std::atomic<int>& my_count) {
        ParallelExecutor<Executor::STD_THREADS> test_block2;
        test_block2(num_threads, TestTask, std::ref(my_count));
      },
      std::ref(count));
  REQUIRE(count == 64);
}

} // namespace qmcplusplus
