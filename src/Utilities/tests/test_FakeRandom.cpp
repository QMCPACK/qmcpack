//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Steven Hahn, hahnse@ornl.gov, Oak Ridge National lab
//
// File created by: Steven Hahn, hahnse@ornl.gov, Oak Ridge National lab
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Utilities/FakeRandom.h"

namespace qmcplusplus
{
TEST_CASE("FakeRandom determinism", "[utilities]")
{
  double expected = 0.5;
  FakeRandom<double> our_rng;
  for (auto i = 0; i < 3; ++i)
    CHECK(our_rng() == Approx(expected));
}

TEST_CASE("FakeRandom set_value", "[utilities]")
{
  FakeRandom<double> our_rng;
  double expected = 0.25;
  our_rng.set_value(expected);
  for (auto i = 0; i < 3; ++i)
    CHECK(our_rng() == Approx(expected));
}

TEST_CASE("FakeRandom read write", "[utilities]")
{
  using DoubleRNG = FakeRandom<double>;
  DoubleRNG rng;

  rng.set_value(1.0);

  std::stringstream stream;
  rng.write(stream);

  DoubleRNG rng2;

  rng2.read(stream);

  for (auto i = 0; i < 3; ++i)
    CHECK(rng2() == rng());
}

TEST_CASE("FakeRandom noops", "[utilities]")
{
  using DoubleRNG = FakeRandom<double>;
  DoubleRNG rng;
  double expected{2.5};

  rng.set_value(expected);
  CHECK(rng() == Approx(expected));

  rng.init(0);
  CHECK(rng() == Approx(expected));

  rng.seed(1);
  CHECK(rng() == Approx(expected));

  std::vector<DoubleRNG::uint_type> state;
  rng.load(state);
  CHECK(state.empty());
  CHECK(rng() == Approx(expected));

  rng.save(state);
  CHECK(rng() == Approx(expected));
}

TEST_CASE("FakeRandom clone", "[utilities]")
{
  using DoubleRNG = FakeRandom<double>;
  DoubleRNG rng;

  std::stringstream stream1;
  rng.write(stream1);

  auto rng2 = rng.makeClone();

  std::stringstream stream2;
  rng2->write(stream2);

  CHECK(stream1.str() == stream2.str());

  for (auto i = 0; i < 3; ++i)
    CHECK((*rng2)() == rng());
}

} // namespace qmcplusplus
