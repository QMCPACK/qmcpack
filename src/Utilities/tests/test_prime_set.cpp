//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, mdewing@anl.gov Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include <iostream>
#include "Message/Communicate.h"
#include "Utilities/PrimeNumberSet.h"
#include <stdio.h>
#include <string>
#include <vector>

using std::string;
using std::vector;
namespace qmcplusplus {

TEST_CASE("prime number set 32 bit", "[utilities]")
{
  PrimeNumberSet<uint32_t> pns;
  //std::cout << "32 bit size = "<< pns.size() << std::endl;
  REQUIRE(pns.size() == 4097);
  REQUIRE(pns[0] == 3);

  std::vector<uint32_t> more_primes;
  // get prime numbers already in the list
  pns.get(1, 2, more_primes);
  REQUIRE(more_primes.size() == 2);
  REQUIRE(more_primes[0] == 5);

  // generate additional prime numbers
  pns.get(4098, 2, more_primes);
  REQUIRE(more_primes.size() == 4);
  REQUIRE(more_primes[2] > pns[4096]);

}

TEST_CASE("prime number set 64 bit", "[utilities]")
{
  PrimeNumberSet<uint64_t> pns;
  //std::cout << "64 bit size = "<< pns.size() << std::endl;
  REQUIRE(pns.size() == 55109);
  REQUIRE(pns[0] == 3);

  std::vector<uint64_t> more_primes;
  // get prime numbers already in the list
  pns.get(1, 2, more_primes);
  REQUIRE(more_primes.size() == 2);
  REQUIRE(more_primes[0] == 5);

  // generate additional prime numbers
  pns.get(55110, 2, more_primes);
  REQUIRE(more_primes.size() == 4);
  REQUIRE(more_primes[2] > pns[55108]);
}

}
