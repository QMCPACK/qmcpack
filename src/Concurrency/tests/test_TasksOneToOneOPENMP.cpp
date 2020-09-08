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

#include "catch.hpp"

#include "Concurrency/TasksOneToOne.hpp"

namespace qmcplusplus
{
/** Openmp generally works but is not guaranteed with std::atomic
 */
void TestTaskOMP(const int ip, int& counter)
{
#pragma omp atomic update
  counter++;
}

TEST_CASE("TasksOneToOne<OPENMP> function case", "[concurrency]")
{
  const int num_threads = omp_get_max_threads();
  TasksOneToOne<Executor::OPENMP> test_block(num_threads);
  int count(0);
  test_block(TestTaskOMP, std::ref(count));
  REQUIRE(count == num_threads);
}

TEST_CASE("TasksOneToOne<OPENMP> lambda case", "[concurrency]")
{
  const int num_threads = omp_get_max_threads();
  TasksOneToOne<Executor::OPENMP> test_block(num_threads);
  int count(0);
  test_block(
      [](int id, int& c) {
#pragma omp atomic update
        c++;
      },
      std::ref(count));
  REQUIRE(count == num_threads);
}

// Sadly this test case is not straight forward with openmp
// Wouldn't be with std::thread either both call terminate when
// the exception is not
//
TEST_CASE("TasksOneToOne<OPENMP> nested case", "[concurrency]")
{
  int num_threads = 1;
  TasksOneToOne<Executor::OPENMP> test_block(num_threads);
  int count(0);
  auto nested_tasks = [num_threads](int task_id, int& my_count) {
    TasksOneToOne<Executor::OPENMP> test_block2(num_threads);
    test_block2(TestTaskOMP, std::ref(my_count));
  };
  REQUIRE_THROWS_WITH(test_block(nested_tasks, std::ref(count)),
                      Catch::Contains("TasksOneToOne should not be used for nested openmp threading"));
}

} // namespace qmcplusplus
