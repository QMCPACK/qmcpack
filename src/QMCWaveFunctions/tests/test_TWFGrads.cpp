//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "TWFGrads.hpp"

namespace qmcplusplus
{
TEST_CASE("TWFGrads", "[QMCWaveFunctions]")
{
  {
    constexpr auto CT = CoordsType::POS;
    auto grads        = TWFGrads<CT>(3);
    REQUIRE(grads.grads_positions.size() == 3);

    grads.grads_positions[0] = QMCTraits::GradType({0.1, 0.2, 0.3});
    grads.grads_positions[1] = QMCTraits::GradType({0.4, 0.5, 0.6});
    grads.grads_positions[2] = QMCTraits::GradType({0.7, 0.8, 0.9});

    auto shift               = TWFGrads<CT>(3);
    shift.grads_positions[0] = QMCTraits::GradType({1.0, 1.0, 1.0});
    shift.grads_positions[1] = QMCTraits::GradType({1.0, 1.0, 1.0});
    shift.grads_positions[2] = QMCTraits::GradType({1.0, 1.0, 1.0});

    grads += shift;
#ifdef QMC_COMPLEX
    CHECK(grads.grads_positions[0][0] == ComplexApprox(QMCTraits::ValueType(1.1, 0.0)));
    CHECK(grads.grads_positions[0][1] == ComplexApprox(QMCTraits::ValueType(1.2, 0.0)));
    CHECK(grads.grads_positions[0][2] == ComplexApprox(QMCTraits::ValueType(1.3, 0.0)));
    CHECK(grads.grads_positions[1][0] == ComplexApprox(QMCTraits::ValueType(1.4, 0.0)));
    CHECK(grads.grads_positions[1][1] == ComplexApprox(QMCTraits::ValueType(1.5, 0.0)));
    CHECK(grads.grads_positions[1][2] == ComplexApprox(QMCTraits::ValueType(1.6, 0.0)));
    CHECK(grads.grads_positions[2][0] == ComplexApprox(QMCTraits::ValueType(1.7, 0.0)));
    CHECK(grads.grads_positions[2][1] == ComplexApprox(QMCTraits::ValueType(1.8, 0.0)));
    CHECK(grads.grads_positions[2][2] == ComplexApprox(QMCTraits::ValueType(1.9, 0.0)));
#else
    CHECK(grads.grads_positions[0][0] == Approx(1.1));
    CHECK(grads.grads_positions[0][1] == Approx(1.2));
    CHECK(grads.grads_positions[0][2] == Approx(1.3));
    CHECK(grads.grads_positions[1][0] == Approx(1.4));
    CHECK(grads.grads_positions[1][1] == Approx(1.5));
    CHECK(grads.grads_positions[1][2] == Approx(1.6));
    CHECK(grads.grads_positions[2][0] == Approx(1.7));
    CHECK(grads.grads_positions[2][1] == Approx(1.8));
    CHECK(grads.grads_positions[2][2] == Approx(1.9));
#endif
  }
  {
    constexpr auto CT = CoordsType::POS_SPIN;
    auto grads        = TWFGrads<CT>(2);
    REQUIRE(grads.grads_positions.size() == 2);
    REQUIRE(grads.grads_spins.size() == 2);
    grads.grads_positions[0] = QMCTraits::GradType({0.1, 0.2, 0.3});
    grads.grads_positions[1] = QMCTraits::GradType({0.4, 0.5, 0.6});
    grads.grads_spins[0]     = QMCTraits::ComplexType(0.7, 0.8);
    grads.grads_spins[1]     = QMCTraits::ComplexType(0.9, 1.0);

    auto shift               = TWFGrads<CT>(2);
    shift.grads_positions[0] = QMCTraits::GradType({1.0, 1.0, 1.0});
    shift.grads_positions[1] = QMCTraits::GradType({1.0, 1.0, 1.0});
    shift.grads_spins[0]     = QMCTraits::ComplexType(1.0, 0.0);
    shift.grads_spins[1]     = QMCTraits::ComplexType(1.0, 0.0);

    grads += shift;
#ifdef QMC_COMPLEX
    CHECK(grads.grads_positions[0][0] == ComplexApprox(QMCTraits::ValueType(1.1, 0.0)));
    CHECK(grads.grads_positions[0][1] == ComplexApprox(QMCTraits::ValueType(1.2, 0.0)));
    CHECK(grads.grads_positions[0][2] == ComplexApprox(QMCTraits::ValueType(1.3, 0.0)));
    CHECK(grads.grads_positions[1][0] == ComplexApprox(QMCTraits::ValueType(1.4, 0.0)));
    CHECK(grads.grads_positions[1][1] == ComplexApprox(QMCTraits::ValueType(1.5, 0.0)));
    CHECK(grads.grads_positions[1][2] == ComplexApprox(QMCTraits::ValueType(1.6, 0.0)));
#else
    CHECK(grads.grads_positions[0][0] == Approx(1.1));
    CHECK(grads.grads_positions[0][1] == Approx(1.2));
    CHECK(grads.grads_positions[0][2] == Approx(1.3));
    CHECK(grads.grads_positions[1][0] == Approx(1.4));
    CHECK(grads.grads_positions[1][1] == Approx(1.5));
    CHECK(grads.grads_positions[1][2] == Approx(1.6));
#endif
    CHECK(grads.grads_spins[0] == ComplexApprox(QMCTraits::ComplexType(1.7, 0.8)));
    CHECK(grads.grads_spins[1] == ComplexApprox(QMCTraits::ComplexType(1.9, 1.0)));
  }
}

} // namespace qmcplusplus
