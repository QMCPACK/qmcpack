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


#include "catch.hpp"
#include "Message/Communicate.h"

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

} // namespace qmcplusplus
