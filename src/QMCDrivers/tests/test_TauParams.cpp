//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "TauParams.hpp"

namespace qmcplusplus
{

TEST_CASE("Taus", "[Particle]")
{
  auto tau       = 1.0;
  auto invmass   = 0.2;
  auto spin_mass = 0.5;

  TauParams<QMCTraits::RealType, CoordsType::POS> taus_rs(tau, invmass, spin_mass);
  CHECK(Approx(taus_rs.tauovermass) == 0.2);
  CHECK(Approx(taus_rs.oneover2tau) == 2.5);
  CHECK(Approx(taus_rs.sqrttau) == 0.447213595499957927703605);

  TauParams<QMCTraits::RealType, CoordsType::POS_SPIN> taus_rsspins(tau, invmass, spin_mass);
  CHECK(Approx(taus_rsspins.spin_tauovermass) == 0.4);
  CHECK(Approx(taus_rsspins.spin_oneover2tau) == 1.25);
  CHECK(Approx(taus_rsspins.spin_sqrttau) == 0.632455532033675882352952);
}

} // namespace qmcplusplus
