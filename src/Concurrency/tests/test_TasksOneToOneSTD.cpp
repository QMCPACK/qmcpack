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

#include "Concurrency/TasksOneToOne.hpp"

namespace qmcplusplus
{
void TestTask(const int ip, std::atomic<int>& counter) { ++counter; }

TEST_CASE("TasksOneToOne<STD> function case", "[concurrency]")
{
  int num_threads = 8;
  TasksOneToOne<Executor::STD> test_block(num_threads);
  std::atomic<int> count(0);
  test_block(TestTask, std::ref(count));
  REQUIRE(count == 8);
}

TEST_CASE("TasksOneToOne<STD> lambda case", "[concurrency]")
{
  int num_threads = 8;
  TasksOneToOne<Executor::STD> test_block(num_threads);
  std::atomic<int> count(0);
  test_block([](int id, std::atomic<int>& my_count) { ++my_count; }, std::ref(count));
  REQUIRE(count == 8);
}

TEST_CASE("TasksOneToOne<STD> nested case", "[concurrency]")
{
  int num_threads = 8;
  TasksOneToOne<Executor::STD> test_block(num_threads);
  std::atomic<int> count(0);
  test_block(
      [num_threads](int task_id, std::atomic<int>& my_count) {
        TasksOneToOne<Executor::STD> test_block2(num_threads);
        test_block2(TestTask, std::ref(my_count));
      },
      std::ref(count));
  REQUIRE(count == 64);
}

} // namespace qmcplusplus
