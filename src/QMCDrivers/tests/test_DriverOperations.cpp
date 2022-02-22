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
#include "DriverOperations.hpp"

namespace qmcplusplus
{
TEST_CASE("DriverOperations", "[QMCDrivers]")
{
  auto tau       = 1.0;
  auto invmass   = 0.2;
  auto spin_mass = 0.5;
  {
    TauParams<QMCTraits::RealType, CoordsType::POS> taus(tau, invmass, spin_mass);
    constexpr auto mct = CoordsType::POS;
    auto mc_coords     = MCCoords<mct>();
    mc_coords.resize(3);
    QMCTraits::PosType p1({0.1, 0.2, 0.3});
    QMCTraits::PosType p2({0.4, 0.5, 0.6});
    QMCTraits::PosType p3({0.7, 0.8, 0.9});
    mc_coords.positions = {p1, p2, p3};

    mc_coords = mc_coords + scaleBySqrtTau(taus, mc_coords);
    CHECK(Approx(mc_coords.positions[0][0]) == 0.14472135955);
    CHECK(Approx(mc_coords.positions[0][1]) == 0.28944271910);
    CHECK(Approx(mc_coords.positions[0][2]) == 0.43416407865);
    CHECK(Approx(mc_coords.positions[1][0]) == 0.57888543820);
    CHECK(Approx(mc_coords.positions[1][1]) == 0.72360679775);
    CHECK(Approx(mc_coords.positions[1][2]) == 0.86832815730);
    CHECK(Approx(mc_coords.positions[2][0]) == 1.01304951685);
    CHECK(Approx(mc_coords.positions[2][1]) == 1.15777087640);
    CHECK(Approx(mc_coords.positions[2][2]) == 1.30249223595);

    std::vector<QMCTraits::RealType> loggf, loggb;
    loggf.resize(3);
    loggb.resize(3);

    updateForwardLogGreensFunction(mc_coords, loggf);
    CHECK(Approx(loggf[0]) == -0.146609903370);
    CHECK(Approx(loggf[1]) == -0.806354468535);
    CHECK(Approx(loggf[2]) == -2.031594375270);

    updateReverseLogGreensFunction(mc_coords, taus, loggb);
    CHECK(Approx(loggb[0]) == -0.733049516850);
    CHECK(Approx(loggb[1]) == -4.031772342675);
    CHECK(Approx(loggb[2]) == -10.15797187635);
  }

  {
    TauParams<QMCTraits::RealType, CoordsType::POS_SPIN> taus(tau, invmass, spin_mass);
    constexpr auto mct = CoordsType::POS_SPIN;
    auto mc_coords     = MCCoords<mct>();
    mc_coords.resize(3);
    QMCTraits::PosType p1({-0.1, -0.2, -0.3});
    QMCTraits::PosType p2({-0.4, -0.5, -0.6});
    QMCTraits::PosType p3({-0.7, -0.8, -0.9});
    mc_coords.positions = {p1, p2, p3};
    mc_coords.spins     = {0.1, 0.2, 0.3};

    mc_coords = mc_coords + scaleBySqrtTau(taus, mc_coords);
    CHECK(Approx(mc_coords.positions[0][0]) == -0.14472135955);
    CHECK(Approx(mc_coords.positions[0][1]) == -0.28944271910);
    CHECK(Approx(mc_coords.positions[0][2]) == -0.43416407865);
    CHECK(Approx(mc_coords.positions[1][0]) == -0.57888543820);
    CHECK(Approx(mc_coords.positions[1][1]) == -0.72360679775);
    CHECK(Approx(mc_coords.positions[1][2]) == -0.86832815730);
    CHECK(Approx(mc_coords.positions[2][0]) == -1.01304951685);
    CHECK(Approx(mc_coords.positions[2][1]) == -1.15777087640);
    CHECK(Approx(mc_coords.positions[2][2]) == -1.30249223595);

    CHECK(Approx(mc_coords.spins[0]) == 0.163245553203);
    CHECK(Approx(mc_coords.spins[1]) == 0.326491106407);
    CHECK(Approx(mc_coords.spins[2]) == 0.489736659610);

    std::vector<QMCTraits::RealType> loggf, loggb;
    loggf.resize(3);
    loggb.resize(3);

    updateForwardLogGreensFunction(mc_coords, loggf);
    CHECK(Approx(loggf[0]) == -0.159934458690);
    CHECK(Approx(loggf[1]) == -0.859652589816);
    CHECK(Approx(loggf[2]) == -2.151515373153);

    updateReverseLogGreensFunction(mc_coords, taus, loggb);
    CHECK(Approx(loggb[0]) == -0.766360905151);
    CHECK(Approx(loggb[1]) == -4.165017895878);
    CHECK(Approx(loggb[2]) == -10.457774371057);
  }
}
} // namespace qmcplusplus
