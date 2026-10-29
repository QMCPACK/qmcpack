//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:  Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
#include <catch2/catch_test_macros.hpp>
#include "Utilities/for_testing/Catch2Approx.h"

#include <memory>
#include "Particle/ParticleSet.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCDrivers/DriftOperators.h"
#include "QMCDrivers/GreenFunctionModifiers/DriftModifierBuilder.h"
#include "QMCDrivers/GreenFunctionModifiers/DriftModifierUNR.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{
TEST_CASE("drift pbyp and node correction real", "[drivers][drift]")
{
  const SimulationCell simulation_cell;
  MCWalkerConfiguration elec(simulation_cell);

  elec.setName("elec");
  std::vector<int> agroup(1);
  agroup[0] = 1;
  elec.create(agroup);

  ParticleSet::RealType tau  = 0.5;
  ParticleSet::RealType mass = 0.85;
  std::vector<ParticleSet::RealType> massinv(1, 1. / mass);
  ParticleSet::ParticlePos drift(1);

  // check from -xtot/2 to xtot/2 in step size of dx i.e. np.arange(-xtot/2,xtot/2,dx)
  double xtot  = 10.;
  int nx       = 100;
  double gradx = -xtot / 2.;
  double dx    = xtot / nx;

  //app_log() << " begin printing" << std::endl;
  for (int ix = 0; ix < nx; ix++)
  {
    elec.G[0][0] = gradx;
    setScaledDriftPbyPandNodeCorr(tau, massinv, elec.G, drift);
    double dval = drift[0][0];

    double scale_factor = (-1. + std::sqrt(1. + 2. * gradx * gradx * tau / mass)) / (gradx * gradx * tau / mass);
    CHECK(dval == Approx(scale_factor * gradx * tau / mass));

    //app_log() << gradx << " " << dval << std::endl;
    gradx += dx;
  }
  //app_log() << " end printing." << std::endl;
}

#ifdef QMC_COMPLEX
TEST_CASE("drift pbyp and node correction complex", "[drivers][drift]")
{ // basically copy and pasted from real test, except "myi"
  const SimulationCell simulation_cell;
  MCWalkerConfiguration elec(simulation_cell);

  elec.setName("elec");
  std::vector<int> agroup(1);
  agroup[0] = 1;
  elec.create(agroup);

  ParticleSet::RealType tau  = 0.5;
  ParticleSet::RealType mass = 0.85;
  std::vector<ParticleSet::RealType> massinv(1, 1. / mass);
  ParticleSet::ParticlePos drift(1);

  // check from -xtot/2 to xtot/2 in step size of dx i.e. np.arange(-xtot/2,xtot/2,dx)
  double xtot  = 10.;
  int nx       = 100;
  double gradx = -xtot / 2.;
  double dx    = xtot / nx;

  // imaginary component of wf gradient should NOT affect drift
  std::complex<double> myi(0, 1.9);
  for (int ix = 0; ix < nx; ix++)
  {
    elec.G[0][0] = gradx + myi;
    setScaledDriftPbyPandNodeCorr(tau, massinv, elec.G, drift);
    double dval = drift[0][0];

    double scale_factor = (-1. + std::sqrt(1. + 2. * gradx * gradx * tau / mass)) / (gradx * gradx * tau / mass);
    CHECK(dval == Approx(scale_factor * gradx * tau / mass));

    gradx += dx;
  }
}
#endif

TEST_CASE("get scaled drift real", "[drivers][drift]")
{
  const SimulationCell simulation_cell;
  MCWalkerConfiguration elec(simulation_cell);

  elec.setName("elec");
  std::vector<int> agroup(1);
  agroup[0] = 1;
  elec.create(agroup);

  ParticleSet::RealType tau  = 0.5;
  ParticleSet::RealType mass = 0.85;
  std::vector<ParticleSet::RealType> massinv(1, 1. / mass);
  ParticleSet::PosType drift;

  // check from -xtot/2 to xtot/2 in step size of dx i.e. np.arange(-xtot/2,xtot/2,dx)
  double xtot  = 10.;
  int nx       = 100;
  double gradx = -xtot / 2.;
  double dx    = xtot / nx;

  DriftModifierUNR DM;
  for (int ix = 0; ix < nx; ix++)
  {
    elec.G[0][0] = gradx;
    DM.getDrift(tau / mass, elec.G[0], drift);
    double dval = drift[0];

    double scale_factor = (-1. + std::sqrt(1. + 2. * gradx * gradx * tau / mass)) / (gradx * gradx * tau / mass);
    CHECK(dval == Approx(scale_factor * gradx * tau / mass));

    gradx += dx;
  }
}

#ifdef QMC_COMPLEX
TEST_CASE("get scaled drift complex", "[drivers][drift]")
{
  const SimulationCell simulation_cell;
  MCWalkerConfiguration elec(simulation_cell);

  elec.setName("elec");
  std::vector<int> agroup(1);
  agroup[0] = 1;
  elec.create(agroup);

  ParticleSet::RealType tau  = 0.5;
  ParticleSet::RealType mass = 0.85;
  std::vector<ParticleSet::RealType> massinv(1, 1. / mass);
  ParticleSet::PosType drift;

  // check from -xtot/2 to xtot/2 in step size of dx i.e. np.arange(-xtot/2,xtot/2,dx)
  double xtot  = 10.;
  int nx       = 100;
  double gradx = -xtot / 2.;
  double dx    = xtot / nx;

  DriftModifierUNR DM;

  // imaginary component of wf gradient should NOT affect drift
  std::complex<double> myi(0, 1.9);
  for (int ix = 0; ix < nx; ix++)
  {
    elec.G[0][0] = gradx + myi;
    DM.getDrift(tau / mass, elec.G[0], drift);
    double dval = drift[0];

    double scale_factor = (-1. + std::sqrt(1. + 2. * gradx * gradx * tau / mass)) / (gradx * gradx * tau / mass);
    CHECK(dval == Approx(scale_factor * gradx * tau / mass));

    gradx += dx;
  }
}
#endif

TEST_CASE("DriftModifierUNR custom parameter and builder", "[drivers][drift]")
{
  using RealType    = QMCTraits::RealType;
  using PosType     = QMCTraits::PosType;
  using GradType    = QMCTraits::GradType;
  using ComplexType = QMCTraits::ComplexType;

  const RealType tau       = 0.5;
  const RealType mass      = 0.85;
  const RealType spin_mass = 2.0;
  const RealType unr_a     = 0.5;

  // Test string/value builder and unknown modifier exception
  CHECK_THROWS_AS(createDriftModifier("invalid_modifier", unr_a), std::runtime_error);
  std::unique_ptr<DriftModifierBase> dm_from_args(createDriftModifier("UNR", unr_a));

  // Test XML builder and parseXML
  const char* xml_str = R"(
  <qmc method="vmc" move="pbyp">
    <parameter name="drift_modifier"> UNR </parameter>
    <parameter name="drift_UNR_a"> 0.5 </parameter>
  </qmc>
)";
  Libxml2Document doc;
  REQUIRE(doc.parseFromString(xml_str));
  xmlNodePtr node = doc.getRoot();
  std::unique_ptr<DriftModifierBase> dm_from_xml(createDriftModifier(node, *OHMMS::Controller));
  REQUIRE(dm_from_xml->parseXML(node));

  const RealType tauovermass = tau / mass;
  const GradType qf_pos(1.5, -2.0, 0.5);
  const RealType vsq = 1.5 * 1.5 + 2.0 * 2.0 + 0.5 * 0.5;
  const RealType sc  = (-1.0 + std::sqrt(1.0 + 2.0 * unr_a * tauovermass * vsq)) / (unr_a * vsq);

  PosType drift_args;
  PosType drift_xml;
  dm_from_args->getDrift(tauovermass, qf_pos, drift_args);
  dm_from_xml->getDrift(tauovermass, qf_pos, drift_xml);

  for (int idim = 0; idim < 3; ++idim)
  {
    CHECK(drift_args[idim] == Approx(sc * std::real(qf_pos[idim])));
    CHECK(drift_xml[idim] == Approx(sc * std::real(qf_pos[idim])));
  }

  // Test zero-gradient branch
  const GradType qf_zero(0.0, 0.0, 0.0);
  PosType drift_zero;
  dm_from_args->getDrift(tauovermass, qf_zero, drift_zero);
  for (int idim = 0; idim < 3; ++idim)
    CHECK(drift_zero[idim] == Approx(0.0));

  // Test scalar (spin) getDrift and getDrifts with POS_SPIN
  const RealType spin_tauovermass = tauovermass / spin_mass;
  const ComplexType qf_spin(1.2, -3.4);
  const RealType spin_vsq = 1.2 * 1.2;
  const RealType spin_sc =
      (-1.0 + std::sqrt(1.0 + 2.0 * unr_a * spin_tauovermass * spin_vsq)) / (unr_a * spin_vsq);

  ParticleSet::Scalar_t spin_drift = 0.0;
  dm_from_args->getDrift(spin_tauovermass, qf_spin, spin_drift);
  CHECK(spin_drift == Approx(spin_sc * std::real(qf_spin)));

  ParticleSet::Scalar_t spin_drift_zero = 1.0;
  dm_from_args->getDrift(spin_tauovermass, ComplexType(0.0, 2.5), spin_drift_zero);
  CHECK(spin_drift_zero == Approx(0.0));

  TauParams<RealType, CoordsType::POS_SPIN> taus(tau, 1.0 / mass, spin_mass);
  TWFGrads<CoordsType::POS_SPIN> grads(2);
  grads.grads_positions[0] = qf_pos;
  grads.grads_positions[1] = qf_zero;
  grads.grads_spins[0]     = qf_spin;
  grads.grads_spins[1]     = ComplexType(0.0, 1.0);

  MCCoords<CoordsType::POS_SPIN> drifts(2);
  dm_from_args->getDrifts(taus, grads, drifts);

  for (int idim = 0; idim < 3; ++idim)
  {
    CHECK(drifts.positions[0][idim] == Approx(sc * std::real(qf_pos[idim])));
    CHECK(drifts.positions[1][idim] == Approx(0.0));
  }
  CHECK(drifts.spins[0] == Approx(spin_sc * std::real(qf_spin)));
  CHECK(drifts.spins[1] == Approx(0.0));
}

} // namespace qmcplusplus
