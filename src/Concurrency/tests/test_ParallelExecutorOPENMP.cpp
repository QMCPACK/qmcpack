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

#include "catch.hpp"

#include "Concurrency/ParallelExecutor.hpp"

namespace qmcplusplus
{
/** Openmp generally works but is not guaranteed with std::atomic
 */
void TestTaskOMP(const int ip, int& counter)
{
#pragma omp atomic update
  counter++;
}

TEST_CASE("ParallelExecutor<OPENMP> function case", "[concurrency]")
{
  const int num_threads = omp_get_max_threads();
  ParallelExecutor<Executor::OPENMP> test_block;
  int count(0);
  test_block(num_threads, TestTaskOMP, std::ref(count));
  REQUIRE(count == num_threads);
}

TEST_CASE("ParallelExecutor<OPENMP> lambda case", "[concurrency]")
{
  const int num_threads = omp_get_max_threads();
  std::cout << "omp_get_max_threads() == " << num_threads << '\n';
  ParallelExecutor<Executor::OPENMP> test_block;
  int count(0);
  test_block(
      num_threads,
      [](int id, int& c) {
#pragma omp atomic update
        c++;
      },
      std::ref(count));
  REQUIRE(count == num_threads);
}

TEST_CASE("ParallelExecutor<OPENMP> nested case", "[concurrency]")
{
  int num_threads = 1;
  ParallelExecutor<Executor::OPENMP> test_block;
  int count(0);
  auto nested_tasks = [num_threads](int task_id, int& my_count) {
    ParallelExecutor<Executor::OPENMP> test_block2;
    test_block2(num_threads, TestTaskOMP, std::ref(my_count));
  };
#ifdef _OPENMP
  REQUIRE_THROWS_WITH(test_block(num_threads, nested_tasks, std::ref(count)),
                      Catch::Contains("ParallelExecutor should not be used for nested openmp threading"));
#endif
}

} // namespace qmcplusplus
