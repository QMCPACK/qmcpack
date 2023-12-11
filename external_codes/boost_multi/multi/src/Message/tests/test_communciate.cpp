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

  if (c->rank() == 0)
  {
    REQUIRE(c2->isGroupLeader() == true);
    auto GroupLeaderComm = c2->getGroupLeaderComm();
    REQUIRE(GroupLeaderComm != nullptr);
    REQUIRE(GroupLeaderComm->size() == 1);
    REQUIRE(GroupLeaderComm->rank() == 0);
  }
  else
  {
    REQUIRE(c2->isGroupLeader() == false);
    REQUIRE(c2->getGroupLeaderComm() == nullptr);
  }

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
    {
      new_rank -= midpoint;
    }
    // Adjust for odd size - the last group has the extra process
    if (c->size() % 2 == 1)
    {
      new_size[1] = new_size[1] + 1;
    }
    REQUIRE(c2->size() == new_size[new_group_id]);
    REQUIRE(c2->rank() == new_rank);

    if (c->rank() == 0 || c->rank() == midpoint)
    {
      REQUIRE(c2->isGroupLeader() == true);
      auto GroupLeaderComm = c2->getGroupLeaderComm();
      REQUIRE(GroupLeaderComm != nullptr);
      REQUIRE(GroupLeaderComm->size() == 2);
      if (c->rank() == 0)
      {
        REQUIRE(GroupLeaderComm->rank() == 0);
      }
      else
      {
        REQUIRE(GroupLeaderComm->rank() == 1);
      }
    }
    else
    {
      REQUIRE(c2->isGroupLeader() == false);
      REQUIRE(c2->getGroupLeaderComm() == nullptr);
    }
  }
}

TEST_CASE("test_communicate_split_four", "[message]")
{
  Communicate* c = OHMMS::Controller;
  // For simplicity, only test the case where the number of processes is divisible by 4.
  if (c->size() % 4 == 0)
  {
    auto c2 = std::make_unique<Communicate>(*c, 4);

    REQUIRE(c2->size() == c->size() / 4);
    int group_size = c->size() / 4;
    int new_rank   = c->rank() % group_size;
    REQUIRE(c2->rank() == new_rank);

    if (new_rank == 0)
    {
      REQUIRE(c2->isGroupLeader() == true);
      auto GroupLeaderComm = c2->getGroupLeaderComm();
      REQUIRE(GroupLeaderComm != nullptr);
      REQUIRE(GroupLeaderComm->size() == 4);
      REQUIRE(GroupLeaderComm->rank() == c->rank() / group_size);
    }
    else
    {
      REQUIRE(c2->isGroupLeader() == false);
      REQUIRE(c2->getGroupLeaderComm() == nullptr);
    }
  }
}

} // namespace qmcplusplus
