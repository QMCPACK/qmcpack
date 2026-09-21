//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
#include <catch2/catch_test_macros.hpp>
#include "Message/Communicate.h"
#include "Message/CommOperators.h"

namespace qmcplusplus
{
TEST_CASE("test_communicate_split_one", "[message]")
{
  Communicate* c = OHMMS::Controller;

  auto c2 = std::make_unique<Communicate>(*c, 1);

  REQUIRE(c2->size() == c->size());
  REQUIRE(c2->rank() == c->rank());

  REQUIRE(c2->isGroupLeader() == (c->rank() == 0));

  auto& inter_group_comm = c2->getInterGroupComm();
  REQUIRE(inter_group_comm.size() == 1);
  REQUIRE(inter_group_comm.rank() == 0);

  std::string real_name = c->getName();

  std::string name = "myname";
  c->setName(name);
  REQUIRE(c->getName() == name);

  std::string other_name = "myothername";
  c->setName(other_name.data(), other_name.size());
  REQUIRE(c->getName() == other_name);

  c->setName(real_name);
}

TEST_CASE("test_communicate_split_two", "[message]")
{
  Communicate* c = OHMMS::Controller;
  if (c->size() >= 2)
  {
    auto c2 = std::make_unique<Communicate>(*c, 2);

    std::vector<int> new_size(2);
    new_size[0] = c->size() / 2;
    new_size[1] = c->size() / 2;

    int midpoint     = c->size() / 2;
    int new_group_id = c->rank() < midpoint ? 0 : 1;
    int new_rank     = c->rank();
    if (c->rank() >= midpoint)
      new_rank -= midpoint;

    // Adjust for odd size - the last group has the extra process
    if (c->size() % 2 == 1)
      new_size[1] = new_size[1] + 1;

    REQUIRE(c2->size() == new_size[new_group_id]);
    REQUIRE(c2->rank() == new_rank);

    REQUIRE(c2->isGroupLeader() == (c->rank() == 0 || c->rank() == midpoint));

    auto& inter_group_comm = c2->getInterGroupComm();
    REQUIRE(inter_group_comm.size() == (c->rank() < midpoint * 2 ? 2 : 1));
    if (c->rank() < midpoint * 2)
      REQUIRE(inter_group_comm.rank() == new_group_id);
    else
      REQUIRE(inter_group_comm.rank() == 0);
  }
}

TEST_CASE("test_communicate_split_four", "[message]")
{
  Communicate* c = OHMMS::Controller;
  // For simplicity, only test the case where the number of processes is divisible by 4.
  if (c->size() % 4 == 0)
  {
    auto c2 = std::make_unique<Communicate>(*c, 4);

    const int group_size = c->size() / 4;
    const int new_rank   = c->rank() % group_size;
    REQUIRE(c2->size() == group_size);
    REQUIRE(c2->rank() == new_rank);
    REQUIRE(c2->isGroupLeader() == (new_rank == 0));
    auto& inter_group_comm = c2->getInterGroupComm();
    REQUIRE(inter_group_comm.size() == 4);
    REQUIRE(inter_group_comm.rank() == c->rank() / group_size);
  }
}

TEST_CASE("test_communicate_split_two_stripe_three", "[message]")
{
  Communicate* c = OHMMS::Controller;
  // For simplicity, only test the case where the number of processes is divisible by 6.
  if (c->size() % 6 != 0)
    return;

  auto c2              = std::make_unique<Communicate>(*c, 2, 3);
  const int group_size = c->size() / 2;
  const int new_rank   = c->rank() % 3 + c->rank() / 6 * 3;
  REQUIRE(c2->size() == group_size);
  REQUIRE(c2->rank() == new_rank);
  REQUIRE(c2->getGroupID() == (c->rank() / 3 % 2));
}

TEST_CASE("test_communicate_complex_pointer_bcast", "[message]")
{
  Communicate* c = OHMMS::Controller;
  std::vector<std::complex<double>> values{{-11.0, -12.0}, {-13.0, -14.0}};
  if (c->rank() == 0)
    values = {{1.5, -2.5}, {3.5, -4.5}};

  c->bcast(values.data(), values.size());

  REQUIRE(values[0] == std::complex<double>{1.5, -2.5});
  REQUIRE(values[1] == std::complex<double>{3.5, -4.5});
}


TEST_CASE("test_communicate_complex_pointer_reduce_in_place", "[message]")
{
  Communicate* c = OHMMS::Controller;
  std::vector<std::complex<double>> values{{(double)(c->rank() + 1), (double)(c->rank() + 2)},
                                           {(double)(c->rank() + 3), (double)(c->rank() + 4)}};

  c->reduce_in_place(values.data(), values.size());

  if (c->rank() == 0)
  {
    double expected_real_0 = 0.0;
    double expected_imag_0 = 0.0;
    double expected_real_1 = 0.0;
    double expected_imag_1 = 0.0;
    for (int i = 0; i < c->size(); i++)
    {
      expected_real_0 += (double)(i + 1);
      expected_imag_0 += (double)(i + 2);
      expected_real_1 += (double)(i + 3);
      expected_imag_1 += (double)(i + 4);
    }
    REQUIRE(values[0] == std::complex<double>{expected_real_0, expected_imag_0});
    REQUIRE(values[1] == std::complex<double>{expected_real_1, expected_imag_1});
  }
}

TEST_CASE("test_communicate_complex_pointer_allgather", "[message]")
{
  Communicate* c = OHMMS::Controller;
  std::vector<std::complex<double>> sb{{(double)(c->rank() + 1), (double)(c->rank() + 2)},
                                       {(double)(c->rank() + 3), (double)(c->rank() + 4)}};
  std::vector<std::complex<double>> rb(c->size() * sb.size());

  c->allgather(sb.data(), rb.data(), sb.size());

  for (int i = 0; i < c->size(); i++)
  {
    double expected_real_0 = (double)(i + 1);
    double expected_imag_0 = (double)(i + 2);
    double expected_real_1 = (double)(i + 3);
    double expected_imag_1 = (double)(i + 4);
    REQUIRE(rb[i * 2 + 0] == std::complex<double>{expected_real_0, expected_imag_0});
    REQUIRE(rb[i * 2 + 1] == std::complex<double>{expected_real_1, expected_imag_1});
  }
}

} // namespace qmcplusplus
