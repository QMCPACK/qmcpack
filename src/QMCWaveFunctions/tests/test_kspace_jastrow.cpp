//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Yubo "Paul Yang", yubo.paul.yang@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Yubo "Paul Yang", yubo.paul.yang@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "OptimizableObject.h"
#include <catch2/catch_test_macros.hpp>
#include "Utilities/for_testing/Catch2Approx.h"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Particle/ParticleSet.h"
#include "Particle/ParticleSetPool.h"
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/Jastrow/kSpaceJastrow.h"
#include "QMCWaveFunctions/Jastrow/kSpaceJastrowBuilder.h"
#include "ParticleIO/LatticeIO.h"

#include <stdio.h>
#include <array>
#include <cmath>
#include <string>

using std::string;

namespace qmcplusplus
{
TEST_CASE("kspace jastrow", "[wavefunction]")
{
  Communicate& c(*OHMMS::Controller);

  // initialize simulationcell for kvectors
  const char* xmltext = R"(<tmp>
  <simulationcell>
     <parameter name="lattice" units="bohr">
              6.00000000        0.00000000        0.00000000
              0.00000000        6.00000000        0.00000000
              0.00000000        0.00000000        6.00000000
     </parameter>
     <parameter name="bconds">
        p p p
     </parameter>
     <parameter name="LR_dim_cutoff"> 15 </parameter>
  </simulationcell>
</tmp>)";
  Libxml2Document doc;
  REQUIRE(doc.parseFromString(xmltext));

  xmlNodePtr root  = doc.getRoot();
  xmlNodePtr part1 = xmlFirstElementChild(root);

  // read lattice
  Lattice lattice;
  LatticeParser lp(lattice);
  lp.put(part1);
  lattice.print(app_log(), 0);

  const SimulationCell simulation_cell(lattice);
  ParticleSet ions_(simulation_cell);
  ParticleSet elec_(simulation_cell);

  ions_.setName("ion");
  ions_.create({1});
  ions_.R[0] = {0.0, 0.0, 0.0};
  elec_.setName("elec");
  elec_.create({2, 0});
  elec_.R[0]                   = {-0.28, 0.0225, -2.709};
  elec_.R[1]                   = {-1.08389, 1.9679, -0.0128914};
  SpeciesSet& tspecies         = elec_.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int downIdx                  = tspecies.addSpecies("d");
  int chargeIdx                = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;
  // initialize SK
  elec_.createSK();

  const char* particles = R"(<tmp>
<jastrow name="Jk" type="kSpace" source="ion">
  <correlation kc="1.5" type="Two-Body" symmetry="isotropic">
    <coefficients id="cG2" type="Array">
      -100. -50.
    </coefficients>
  </correlation>
</jastrow>
</tmp>
)";
  REQUIRE(doc.parseFromString(particles));

  root            = doc.getRoot();
  xmlNodePtr jas1 = xmlFirstElementChild(root);

  kSpaceJastrowBuilder jastrow(c, elec_, ions_);
  std::unique_ptr<WaveFunctionComponent> jas(jastrow.buildComponent(jas1));

  // update all distance tables
  elec_.update();

  double logpsi_real = std::real(jas->evaluateLog(elec_, elec_.G, elec_.L));
  CHECK(logpsi_real == Approx(-4.4088303951));

  // Cached accepted/rejected moves must match independent full evaluations.
#if defined(MIXED_PRECISION)
  const double tolerance = 1e-5;
#else
  const double tolerance = 1e-11;
#endif
  ParticleSet reference(simulation_cell);
  reference.setName("reference");
  reference.create({2, 0});
  reference.R = elec_.R;
  reference.update();
  kSpaceJastrowBuilder reference_builder(c, reference, ions_);
  auto reference_jas = reference_builder.buildComponent(jas1);
  // Independently expand |rho_G|^2 = 2 + 2 cos(G . (r_0 - r_1)).
  // The cutoff includes the first two shells, with one vector from each +/- pair.
  const std::array<std::array<int, 3>, 9> reciprocal_indices = {
      {{0, 0, 1}, {0, 1, -1}, {0, 1, 0}, {0, 1, 1}, {1, -1, 0}, {1, 0, -1}, {1, 0, 0}, {1, 0, 1}, {1, 1, 0}}};
  const auto check_full_evaluation = [&](double log_value) {
    const double wave_number = 2.0 * M_PI / 6.0;
    double expected_log = 0.0, expected_laplacian = 0.0;
    std::array<double, 3> expected_gradient{};
    for (const auto& index : reciprocal_indices)
    {
      const int shell          = index[0] * index[0] + index[1] * index[1] + index[2] * index[2];
      const double coefficient = (shell == 1 ? -100.0 : -50.0) / 216.0;
      double phase             = 0.0;
      for (int dim = 0; dim < 3; ++dim)
        phase += wave_number * index[dim] * (reference.R[0][dim] - reference.R[1][dim]);
      expected_log += coefficient * (2.0 + 2.0 * std::cos(phase));
      expected_laplacian -= 2.0 * coefficient * wave_number * wave_number * shell * std::cos(phase);
      for (int dim = 0; dim < 3; ++dim)
        expected_gradient[dim] -= 2.0 * coefficient * wave_number * index[dim] * std::sin(phase);
    }
    CHECK(log_value == Approx(expected_log).margin(tolerance));
    for (int particle = 0; particle < 2; ++particle)
    {
      const double sign = particle == 0 ? 1.0 : -1.0;
      for (int dim = 0; dim < 3; ++dim)
        CHECK(std::real(reference.G[particle][dim]) == Approx(sign * expected_gradient[dim]).margin(tolerance));
      CHECK(std::real(reference.L[particle]) == Approx(expected_laplacian).margin(tolerance));
    }
  };
  CHECK(std::real(jas->get_log_value()) == Approx(logpsi_real));
  for (int step = 0; step < 8; ++step)
  {
    CAPTURE(step);
    const int iat = step % 2;
    reference.R   = elec_.R;
    reference.update();
    reference.G          = 0;
    reference.L          = 0;
    const double old_log = std::real(reference_jas->evaluateLog(reference, reference.G, reference.L));
    check_full_evaluation(old_log);
    const ParticleSet::PosType displacement(0.11 + 0.02 * step, -0.09, 0.07);
    elec_.makeMove(iat, displacement);
    reference.R[iat] += displacement;
    reference.update();
    reference.G          = 0;
    reference.L          = 0;
    const double new_log = std::real(reference_jas->evaluateLog(reference, reference.G, reference.L));
    check_full_evaluation(new_log);
    const double expected_ratio = std::exp(new_log - old_log);
    if (step % 2 == 0)
      CHECK(std::real(jas->ratio(elec_, iat)) == Approx(expected_ratio).epsilon(tolerance));
    else
    {
      WaveFunctionComponent::GradType proposed_gradient;
      proposed_gradient = 0;
      CHECK(std::real(jas->ratioGrad(elec_, iat, proposed_gradient)) == Approx(expected_ratio).epsilon(tolerance));
      for (int dim = 0; dim < 3; ++dim)
        CHECK(std::real(proposed_gradient[dim]) == Approx(std::real(reference.G[iat][dim])).margin(tolerance));
    }
    if (step % 3 == 2)
    {
      jas->restore(iat);
      elec_.rejectMove(iat);
      CHECK(std::real(jas->get_log_value()) == Approx(old_log).margin(tolerance));
    }
    else
    {
      jas->acceptMove(elec_, iat);
      elec_.acceptMove(iat);
      CHECK(std::real(jas->get_log_value()) == Approx(new_log).margin(tolerance));
    }
  }
}

TEST_CASE("kspace jastrow derivatives", "[wavefunction]")
{
  Communicate& c(*OHMMS::Controller);

  // initialize simulationcell for kvectors
  const char* xmltext = R"(<tmp>
  <simulationcell>
     <parameter name="lattice" units="bohr">
              1.00000000        0.00000000        0.00000000
             -0.27834620        1.98053614        0.00000000
              0.31358539        1.32438627        2.67351177
     </parameter>
     <parameter name="bconds">
        p p p
     </parameter>
     <parameter name="LR_dim_cutoff"> 15 </parameter>
  </simulationcell>
</tmp>)";
  Libxml2Document doc;
  REQUIRE(doc.parseFromString(xmltext));

  xmlNodePtr root  = doc.getRoot();
  xmlNodePtr part1 = xmlFirstElementChild(root);

  // read lattice
  Lattice lattice;
  LatticeParser lp(lattice);
  lp.put(part1);
  lattice.print(app_log(), 0);

  const SimulationCell simulation_cell(lattice);
  ParticleSet ions_(simulation_cell);
  ParticleSet elec_(simulation_cell);

  ions_.setName("ion");
  ions_.create({1});
  ions_.R[0] = {0.7, 0.8, 0.9};
  elec_.setName("elec");
  elec_.create({2, 0});
  elec_.R[0]                   = {0.1, 0.2, 0.3};
  elec_.R[1]                   = {0.4, 0.5, 0.6};
  SpeciesSet& tspecies         = elec_.getSpeciesSet();
  int upIdx                    = tspecies.addSpecies("u");
  int downIdx                  = tspecies.addSpecies("d");
  int chargeIdx                = tspecies.addAttribute("charge");
  tspecies(chargeIdx, upIdx)   = -1;
  tspecies(chargeIdx, downIdx) = -1;
  // initialize SK
  elec_.createSK();


  const char* jk1input = R"(<tmp>
<jastrow name="Jk" type="kSpace" source="ion0">
  <correlation kc="5.0" type="One-Body" symmetry="crystal">
    <coefficients id="cG1" type="Array">
0.8137632591137914 0.44470998182408406  0.5401857355121155 0.25752329923480577  0.8232477489081528 0.44767385249164526 0.6464929258779923 0.4929171993220475  0.7846576485476944 0.662660459182017
    </coefficients>
  </correlation>
</jastrow>
</tmp>
)";
  REQUIRE(doc.parseFromString(jk1input));

  root            = doc.getRoot();
  xmlNodePtr jas1 = xmlFirstElementChild(root);

  kSpaceJastrowBuilder jastrow(c, elec_, ions_);
  std::unique_ptr<WaveFunctionComponent> wfc(jastrow.buildComponent(jas1));
  kSpaceJastrow* jas = dynamic_cast<kSpaceJastrow*>(wfc.get());

  // update all distance tables
  elec_.update();

  double logpsi_real = std::real(jas->evaluateLog(elec_, elec_.G, elec_.L));
  CHECK(logpsi_real == Approx(0.7137163755813973));
  CHECK(std::real(elec_.G[0][0]) == Approx(0.0));
  CHECK(std::real(elec_.G[0][1]) == Approx(-0.35067397));
  CHECK(std::real(elec_.G[0][2]) == Approx(-1.18746358));
  CHECK(std::real(elec_.G[1][0]) == Approx(0.0));
  CHECK(std::real(elec_.G[1][1]) == Approx(-0.5971649));
  CHECK(std::real(elec_.G[1][2]) == Approx(-1.30622405));
  CHECK(std::real(elec_.L[0]) == Approx(-9.23735526));
  CHECK(std::real(elec_.L[1]) == Approx(2.37396869));

  OptVariables opt_vars;
  jas->checkInVariablesExclusive(opt_vars);
  opt_vars.resetIndex();
  jas->checkOutVariables(opt_vars);
  jas->resetParametersExclusive(opt_vars);

  const int nopt = opt_vars.size();
  Vector<ParticleSet::ValueType> dlogpsi(nopt);
  Vector<ParticleSet::ValueType> dhpsioverpsi(nopt);

  std::vector<double> refvals1 = {0.17404927, 0.30881786,  0.03441212, 0.31140348,  0.3374628,
                                  0.14393708, -0.17959687, 0.12910817, -0.14896848, 0.2460836};
  std::vector<double> refvals2 = {-0.433200292526692,  1.3469998192405797, -0.5970340626326731, 1.6549156762746264,
                                  2.155972109951776,   0.7377064401570266, -2.33721875453341,   0.20791372534745278,
                                  -3.0523924081649056, 1.7868809500942717};


  dlogpsi = 0;
  jas->evaluateDerivativesWF(elec_, opt_vars, dlogpsi);
  for (int i = 0; i < nopt; i++)
    CHECK(std::real(dlogpsi[i]) == Approx(refvals1[i]));

  dlogpsi = 0;
  jas->evaluateDerivatives(elec_, opt_vars, dlogpsi, dhpsioverpsi);
  for (int i = 0; i < nopt; i++)
    CHECK(std::real(dlogpsi[i]) == Approx(refvals1[i]));
  for (int i = 0; i < nopt; i++)
    CHECK(std::real(dhpsioverpsi[i]) == Approx(refvals2[i]));

  // Twobody check
  const char* jk2input = R"(<tmp>
<jastrow name="Jk" type="kSpace" source="ion0">
  <correlation kc="5.0" type="Two-Body" symmetry="crystal">
    <coefficients id="cG2" type="Array">
    0.9531019536367156  0.18148850587408794  0.20772539137877666  0.9340655612889098  0.50544913283957
    </coefficients>
  </correlation>
</jastrow>
</tmp>
)";
  REQUIRE(doc.parseFromString(jk2input));

  root            = doc.getRoot();
  xmlNodePtr jas2 = xmlFirstElementChild(root);

  kSpaceJastrowBuilder jastrow2(c, elec_, ions_);
  std::unique_ptr<WaveFunctionComponent> wfc2(jastrow2.buildComponent(jas2));
  kSpaceJastrow* j2 = dynamic_cast<kSpaceJastrow*>(wfc2.get());

  // update all distance tables
  elec_.update();
  elec_.G = 0;
  elec_.L = 0;
  //reference values from python code in test directory
  logpsi_real = std::real(j2->evaluateLog(elec_, elec_.G, elec_.L));
  CHECK(logpsi_real == Approx(1.3399793683));
  CHECK(std::real(elec_.G[0][0]) == Approx(0.0));
  CHECK(std::real(elec_.G[0][1]) == Approx(1.37913407));
  CHECK(std::real(elec_.G[0][2]) == Approx(2.47457664));
  CHECK(std::real(elec_.G[1][0]) == Approx(0.0));
  CHECK(std::real(elec_.G[1][1]) == Approx(-1.37913407));
  CHECK(std::real(elec_.G[1][2]) == Approx(-2.47457664));
  CHECK(std::real(elec_.L[0]) == Approx(-1.13586493));
  CHECK(std::real(elec_.L[1]) == Approx(-1.13586493));

  OptVariables opt_vars2;
  j2->checkInVariablesExclusive(opt_vars2);
  opt_vars2.resetIndex();
  j2->checkOutVariables(opt_vars2);
  j2->resetParametersExclusive(opt_vars2);

  const int nopt2 = opt_vars2.size();
  Vector<ParticleSet::ValueType> dlogpsi2(nopt2);
  Vector<ParticleSet::ValueType> dhpsioverpsi2(nopt2);

  std::vector<double> refvals3 = {0.66537659, 0.51973641, 0.71270004, 0.25905169, 0.4381535};
  std::vector<double> refvals4 = {-1.2583630528695267, -2.8959032092323866, 4.02906470009512, -11.046492939481567,
                                  -7.338195726624974};

  dlogpsi2 = 0.0;
  j2->evaluateDerivativesWF(elec_, opt_vars2, dlogpsi2);
  for (int i = 0; i < nopt2; i++)
    CHECK(std::real(dlogpsi2[i]) == Approx(refvals3[i]));

  dlogpsi2      = 0.0;
  dhpsioverpsi2 = 0.0;
  j2->evaluateDerivatives(elec_, opt_vars2, dlogpsi2, dhpsioverpsi2);
  for (int i = 0; i < nopt2; i++)
    CHECK(std::real(dlogpsi2[i]) == Approx(refvals3[i]));
  for (int i = 0; i < nopt2; i++)
    CHECK(std::real(dhpsioverpsi2[i]) == Approx(refvals4[i]));
}
} // namespace qmcplusplus
