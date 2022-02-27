//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Utilities/FakeRandom.h"
#include "Utilities/StdRandom.h"
#include "Message/Communicate.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/TinyVector.h"
#include <vector>

using std::cerr;
using std::endl;
using std::vector;
#include "ParticleBase/RandomSeqGenerator.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
TEST_CASE("gaussian random array length 1", "[particle_base]")
{
  FakeRandom rg;
  std::vector<double> a(2);
  assignGaussRand(a.data(), 1, rg);

  // assuming RNG input is 0.5
  REQUIRE(a[0] == Approx(-1.1774100224305424));
  REQUIRE(a[1] == Approx(0.0)); // ensure no overflow
}

TEST_CASE("gaussian random array length 2", "[particle_base]")
{
  FakeRandom rg;
  std::vector<double> a(3);
  assignGaussRand(a.data(), 2, rg);

  // assuming RNG input is 0.5
  REQUIRE(a[0] == Approx(-1.1774100224305424));
  REQUIRE(a[1] == Approx(1.4419114152535772e-16));
  REQUIRE(a[2] == Approx(0.0)); // ensure no overflow
}

TEST_CASE("gaussian random array length 3", "[particle_base]")
{
  FakeRandom rg;
  std::vector<double> a(4);
  assignGaussRand(a.data(), 3, rg);

  // assuming RNG input is 0.5
  REQUIRE(a[0] == Approx(-1.1774100224305424));
  REQUIRE(a[1] == Approx(1.4419114152535772e-16));
  REQUIRE(a[2] == Approx(-1.1774100224305424));
  REQUIRE(a[3] == Approx(0.0)); // ensure no overflow
}

TEST_CASE("gaussian random particle attrib array length 1", "[particle_base]")
{
  FakeRandom rg;
  ParticleAttrib<TinyVector<double, 1>> PA;
  PA.resize(1);
  makeGaussRandomWithEngine(PA, rg);

  // assuming RNG input is 0.5
  REQUIRE(PA[0][0] == Approx(-1.1774100224305424));
}

TEST_CASE("gaussian random input one", "[particle_base]")
{
  FakeRandom rg;
  rg.set_value(1.0);
  std::vector<double> a(2);
  assignGaussRand(a.data(), 2, rg);

  // uniform RNG input is 1.0
  // most uniform RNGs do not produce 1.0 exactly (open interval),
  // but the code is there to prevent it, so a test.
  REQUIRE(a[0] == Approx(8.49042441685));
  REQUIRE(a[1] == Approx(0.0)); // ensure no overflow
}

TEST_CASE("gaussian random input zero", "[particle_base]")
{
  FakeRandom rg;
  rg.set_value(0.0);
  std::vector<double> a(2);
  assignGaussRand(a.data(), 2, rg);

  // uniform RNG input is 0.0
  REQUIRE(a[0] == Approx(0.0));
  REQUIRE(a[1] == Approx(0.0));
}

TEST_CASE("makeGaussRandomWithEngine(MCCoords...)", "[particle_base]")
{
  int size_test = 7;
  std::vector<double> gauss_random_vals(size_test * 3 + (size_test * 3) % 2 + size_test);
  {
    StdRandom<double> rng;
    makeGaussRandomWithEngine(gauss_random_vals, rng);
  }

  auto checkRs = [&](auto& rs) {
    for (int i = 0; i < size_test; ++i)
    {
      CHECK(Approx(gauss_random_vals[3 * i]) == rs[i][0]);
      CHECK(Approx(gauss_random_vals[3 * i + 1]) == rs[i][1]);
      CHECK(Approx(gauss_random_vals[3 * i + 2]) == rs[i][2]);
    }
  };

  MCCoords<CoordsType::POS> mc_coords_rs(size_test);
  {
    StdRandom<double> rng;
    makeGaussRandomWithEngine(mc_coords_rs, rng);
    checkRs(mc_coords_rs.positions);
  }
  MCCoords<CoordsType::POS_SPIN> mc_coords_rsspins(size_test);
  {
    StdRandom<double> rng;
    makeGaussRandomWithEngine(mc_coords_rsspins, rng);
    checkRs(mc_coords_rsspins.positions);
    // Mod 2 is result of how gaussianDistribution is generated.
    int offset_for_rs = (3 * size_test) + (3 * size_test) % 2;
    for (int i = 0; i < size_test; ++i)
      CHECK(Approx(gauss_random_vals[offset_for_rs + i]) == mc_coords_rsspins.spins[i]);
  }
}


} // namespace qmcplusplus
