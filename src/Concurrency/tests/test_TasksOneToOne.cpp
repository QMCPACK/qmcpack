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

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "Concurrency/TasksOneToOne.hpp"

namespace qmcplusplus
{
void TestTask(const int ip, std::atomic<int>& counter) { counter++; }

/** Openmp generally works but is not guaranteed with std::atomic
 */
void TestTaskOMP(const int ip, int& counter)
{
#pragma omp atomic update
  counter++;
}

TEST_CASE("TasksOneToOne<OPENMP> function case", "[concurrency]")
{
  int num_threads = 8;
  TasksOneToOne<Threading::OPENMP> test_block(num_threads);
  int count(0);
  test_block(TestTaskOMP, std::ref(count));
  REQUIRE(count == 8);
}

TEST_CASE("TasksOneToOne<OPENMP> lambda case", "[concurrency]")
{
  int num_threads = 8;
  TasksOneToOne<Threading::OPENMP> test_block(num_threads);
  int count(0);
  test_block(
      [](int id, int& c) {
#pragma omp atomic update
        c++;
      },
      std::ref(count));
  REQUIRE(count == 8);
}

// Sadly this test case is not straight forward with openmp
// Wouldn't be with std::thread either both call terminate when
// the exception is not
//
TEST_CASE("TasksOneToOne<OPENMP> nested case", "[concurrency]")
{
  int num_threads = 1;
  TasksOneToOne<Threading::OPENMP> test_block(num_threads);
  int count(0);
  auto nested_tasks = [num_threads](int task_id, int& my_count) {
    TasksOneToOne<Threading::OPENMP> test_block2(num_threads);
    test_block2(TestTaskOMP, std::ref(my_count));
  };
  bool threw = false;
  REQUIRE_THROWS_WITH(test_block(nested_tasks, std::ref(count)),
                      Catch::Contains("TasksOneToOne should not be used for nested openmp threading"));
}

TEST_CASE("TasksOneToOne<STD> function case", "[concurrency]")
{
  int num_threads = 8;
  TasksOneToOne<Threading::STD> test_block(num_threads);
  std::atomic<int> count(0);
  test_block(TestTask, std::ref(count));
  REQUIRE(count == 8);
}

TEST_CASE("TasksOneToOne<STD> lambda case", "[concurrency]")
{
  int num_threads = 8;
  TasksOneToOne<Threading::STD> test_block(num_threads);
  std::atomic<int> count(0);
  test_block([](int id, auto& my_count) { my_count++; }, std::ref(count));
  REQUIRE(count == 8);
}

TEST_CASE("TasksOneToOne<STD> nested case", "[concurrency]")
{
  int num_threads = 8;
  TasksOneToOne<Threading::STD> test_block(num_threads);
  int count(0);
  test_block(
      [num_threads](int task_id, int& my_count) {
        TasksOneToOne<Threading::STD> test_block2(num_threads);
        test_block2(TestTaskOMP, std::ref(my_count));
      },
      std::ref(count));
  REQUIRE(count == 64);
}

} // namespace qmcplusplus
