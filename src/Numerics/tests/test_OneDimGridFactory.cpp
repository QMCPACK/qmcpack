//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: QMCPACK developers
//////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include "Utilities/for_testing/Catch2Approx.h"

#include <cmath>

#include "Numerics/OneDimGridFactory.h"
#include "OhmmsData/Libxml2Doc.h"

namespace qmcplusplus
{
namespace
{
std::unique_ptr<OneDimGridFactory::GridType> makeGrid(const char* grid_xml)
{
  Libxml2Document doc;
  REQUIRE(doc.parseFromString(grid_xml));
  return OneDimGridFactory::createGrid(doc.getRoot());
}
} // namespace

TEST_CASE("OneDimGridFactory input controls", "[numerics]")
{
  SECTION("omitted-input defaults")
  {
    auto grid = OneDimGridFactory::createGrid(nullptr);

    REQUIRE(grid != nullptr);
    CHECK(grid->getGridTag() == LOG_1DGRID);
    CHECK(grid->size() == 1001);
    CHECK(grid->rmin() == Approx(1.0e-5));
    CHECK(grid->rmax() == Approx(100.0));
  }

  SECTION("bounded logarithmic grid and retired identity attributes")
  {
    auto grid = makeGrid(R"(<grid type="log" ri="0.25" rf="4.0" npts="5"
                                  id="unused" name="unused" ref="unused"/>)");

    REQUIRE(grid != nullptr);
    CHECK(grid->getGridTag() == LOG_1DGRID);
    CHECK(grid->size() == 5);
    CHECK(grid->rmin() == Approx(0.25));
    CHECK(grid->rmax() == Approx(4.0));
    CHECK(grid->r(2) == Approx(1.0));
  }

  SECTION("origin logarithmic grid with ascale and astep")
  {
    auto grid = makeGrid(R"(<grid type="log" ascale="2.0" astep="0.25" npts="4"/>)");

    REQUIRE(grid != nullptr);
    CHECK(grid->getGridTag() == LOGZERO_1DGRID);
    CHECK(grid->size() == 4);
    CHECK(grid->r(0) == Approx(0.0));
    CHECK(grid->r(1) == Approx(2.0 * (std::exp(0.25) - 1.0)));
  }

  SECTION("origin logarithmic grid with scale and step aliases")
  {
    auto grid = makeGrid(R"(<grid type="log" scale="3.0" step="0.5" npts="3"/>)");

    REQUIRE(grid != nullptr);
    CHECK(grid->getGridTag() == LOGZERO_1DGRID);
    CHECK(grid->size() == 3);
    CHECK(grid->r(1) == Approx(3.0 * (std::exp(0.5) - 1.0)));
  }

  SECTION("linear grid")
  {
    auto grid = makeGrid(R"(<grid type="linear" ri="0.5" rf="2.5" npts="5"/>)");

    REQUIRE(grid != nullptr);
    CHECK(grid->getGridTag() == LINEAR_1DGRID);
    CHECK(grid->size() == 5);
    CHECK(grid->rmin() == Approx(0.5));
    CHECK(grid->rmax() == Approx(2.5));
    CHECK(grid->r(2) == Approx(1.5));
  }
}
} // namespace qmcplusplus
