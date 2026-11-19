//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "Platforms/Host/OutputManager.h"

#include "Concurrency/ParallelExecutor.hpp"

namespace qmcplusplus
{
/** Openmp generally works but is not guaranteed with std::atomic
 */
void TestTaskOMP(const int /*ip*/, int& counter)
{
#pragma omp atomic update
  counter++;
}

TEST_CASE("ParallelExecutor<OPENMP> function case", "[concurrency]")
{
  const int num_threads = omp_get_max_threads();
  ParallelExecutor<Executor::OPENMP> test_block;
  int count(0);
  test_block(num_threads, TestTaskOMP, count);
  REQUIRE(count == num_threads);
}

TEST_CASE("ParallelExecutor<OPENMP> lambda case", "[concurrency]")
{
  const int num_threads = omp_get_max_threads();
  app_log() << "omp_get_max_threads() == " << num_threads << '\n';
  ParallelExecutor<Executor::OPENMP> test_block;
  int count(0);
  test_block(
      num_threads,
      [](int /*id*/, int& c) {
#pragma omp atomic update
        c++;
      },
      count);
  REQUIRE(count == num_threads);
}

TEST_CASE("ParallelExecutor<OPENMP> nested case", "[concurrency]")
{
  int num_threads = 1;
  ParallelExecutor<Executor::OPENMP> test_block;
  int count(0);
  auto nested_tasks = [num_threads](int /*task_id*/, int& my_count) {
    ParallelExecutor<Executor::OPENMP> test_block2;
    test_block2(num_threads, TestTaskOMP, my_count);
  };
#ifdef _OPENMP
  REQUIRE_THROWS_WITH(test_block(num_threads, nested_tasks, count),
                      Catch::Matchers::ContainsSubstring("Unexpected exception thrown in threaded section: task_id 0") &&
                          Catch::Matchers::ContainsSubstring("ParallelExecutor should not be used for nested openmp threading"));
#endif
}

TEST_CASE("ParallelExecutor<OPENMP> exception handling and iteration skipping", "[concurrency]")
{
  ParallelExecutor<Executor::OPENMP> test_block;
  const int total_tasks = 100;
  int tasks_executed    = 0;

  auto failing_task = [&tasks_executed](int task_id) {
#pragma omp atomic update
    tasks_executed++;
    if (task_id == 0)
      throw std::runtime_error("Task 0 failed deliberately");
  };

  // A single thread runs task 0 first, making the skipping of the remaining iterations deterministic.
  const int saved_max_threads = omp_get_max_threads();
  omp_set_num_threads(1);
  REQUIRE_THROWS_WITH(test_block(total_tasks, failing_task),
                      Catch::Matchers::ContainsSubstring("Unexpected exception thrown in threaded section: task_id 0") &&
                          Catch::Matchers::ContainsSubstring("Task 0 failed deliberately"));
  omp_set_num_threads(saved_max_threads);
#ifdef _OPENMP
  CHECK(tasks_executed == 1);
#endif
}

TEST_CASE("ParallelExecutor<OPENMP> multiple concurrent exceptions", "[concurrency]")
{
  ParallelExecutor<Executor::OPENMP> test_block;
  const int total_tasks = 100;

  // Every task throws. Only the first one is rethrown, the others go to std::cerr.
  auto failing_task = [](int /*task_id*/) { throw std::runtime_error("Task failed deliberately"); };

  REQUIRE_THROWS_WITH(test_block(total_tasks, failing_task),
                      Catch::Matchers::ContainsSubstring("Unexpected exception thrown in threaded section: task_id") &&
                          Catch::Matchers::ContainsSubstring("Task failed deliberately"));
}

TEST_CASE("ParallelExecutor<OPENMP> non-std exception", "[concurrency]")
{
  ParallelExecutor<Executor::OPENMP> test_block;

  auto failing_task = [](int task_id) {
    if (task_id == 0)
      throw 42;
  };

  REQUIRE_THROWS_WITH(test_block(10, failing_task),
                      Catch::Matchers::ContainsSubstring("Unexpected exception thrown in threaded section: task_id 0") &&
                          Catch::Matchers::ContainsSubstring("not derived from std::exception"));
}

} // namespace qmcplusplus
