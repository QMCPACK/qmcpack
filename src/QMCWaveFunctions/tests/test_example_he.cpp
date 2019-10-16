//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Configuration.h"
#include "Message/Communicate.h"
#include "OhmmsData/Libxml2Doc.h"
#include "QMCWaveFunctions/WaveFunctionFactory.h"
#include "QMCWaveFunctions/ExampleHeComponent.h"

namespace qmcplusplus
{
typedef WaveFunctionComponent::RealType RealType;
typedef WaveFunctionComponent::ValueType ValueType;

TEST_CASE("ExampleHe", "[wavefunction]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate* c = OHMMS::Controller;

  ParticleSet* elec = new ParticleSet;
  std::vector<int> agroup(1);
  int nelec = 2;
  agroup[0] = nelec;
  elec->setName("e");
  elec->create(agroup);
  elec->R[0][0] = 1.0;
  elec->R[0][1] = 2.0;
  elec->R[0][2] = 3.0;
  elec->R[1][0] = 0.0;
  elec->R[1][1] = 1.1;
  elec->R[1][2] = 2.2;

  SpeciesSet& tspecies       = elec->getSpeciesSet();
  int upIdx                  = tspecies.addSpecies("u");
  int downIdx                = tspecies.addSpecies("d");
  int massIdx                = tspecies.addAttribute("mass");
  tspecies(massIdx, upIdx)   = 1.0;
  tspecies(massIdx, downIdx) = 1.0;


  WaveFunctionFactory::PtclPoolType particle_set_map;
  particle_set_map["e"] = elec;

  ParticleSet* ions = new ParticleSet;
  ions->setName("ion0");
  ions->create(1);
  ions->R[0][0] = 0.0;
  ions->R[0][1] = 0.0;
  ions->R[0][2] = 0.0;

  SpeciesSet& he_species      = ions->getSpeciesSet();
  int He_Idx                  = he_species.addSpecies("He");
  int chargeIdx               = he_species.addAttribute("charge");
  tspecies(chargeIdx, He_Idx) = 2.0;
  tspecies(massIdx, upIdx)    = 1.0;
  particle_set_map["ion0"]    = ions;


  WaveFunctionFactory wff(elec, particle_set_map, c);

  const char* wavefunction_xml = "<wavefunction> \
  <example_he name=\"mine\" source=\"ion0\"> \
        <var id=\"B\" name=\"B\">0.8</var> \
  </example_he> \
</wavefunction>";
  Libxml2Document doc;
  bool okay = doc.parseFromString(wavefunction_xml);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();
  wff.put(root);

  REQUIRE(wff.targetPsi != NULL);
  REQUIRE(wff.targetPsi->size() == 1);

  WaveFunctionComponent* base_example_he = wff.targetPsi->getOrbitals()[0];
  REQUIRE(base_example_he != NULL);

  ExampleHeComponent* example_he = dynamic_cast<ExampleHeComponent*>(base_example_he);
  REQUIRE(example_he != NULL);

  ions->update();
  elec->update();

  ParticleSet::ParticleGradient_t all_grad;
  ParticleSet::ParticleLaplacian_t all_lap;
  all_grad.resize(nelec);
  all_lap.resize(nelec);

  // Set the base expectations for wavefunction value and derivatives
  RealType logpsi = example_he->evaluateLog(*elec, all_grad, all_lap);


  // Comparisons are performed at a single set of electron coordinates.  This should be expanded.

  // Compare with evalGrad

  ParticleSet::GradType grad0;
  int iat = 0;
  grad0   = example_he->evalGrad(*elec, iat);

  REQUIRE(grad0[0] == ValueApprox(all_grad[0][0]));
  REQUIRE(grad0[1] == ValueApprox(all_grad[0][1]));
  REQUIRE(grad0[2] == ValueApprox(all_grad[0][2]));

  ParticleSet::GradType grad1;
  iat   = 1;
  grad1 = example_he->evalGrad(*elec, iat);

  REQUIRE(grad1[0] == ValueApprox(all_grad[1][0]));
  REQUIRE(grad1[1] == ValueApprox(all_grad[1][1]));
  REQUIRE(grad1[2] == ValueApprox(all_grad[1][2]));


  // Compare ratio and ratioGrad with a zero displacement
  ParticleSet::SingleParticlePos_t zero_displ(0.0, 0.0, 0.0);
  iat = 0;
  elec->makeMove(iat, zero_displ);


  ValueType ratio = example_he->ratio(*elec, iat);
  REQUIRE(std::real(ratio) == Approx(1.0));

  ratio = example_he->ratioGrad(*elec, iat, grad0);

  REQUIRE(std::real(ratio) == Approx(1.0));

  REQUIRE(grad0[0] == ValueApprox(all_grad[0][0]));
  REQUIRE(grad0[1] == ValueApprox(all_grad[0][1]));
  REQUIRE(grad0[2] == ValueApprox(all_grad[0][2]));

  iat = 1;
  elec->makeMove(iat, zero_displ);
  ratio = example_he->ratio(*elec, iat);
  REQUIRE(std::real(ratio) == Approx(1.0));


  ratio = example_he->ratioGrad(*elec, iat, grad1);

  REQUIRE(std::real(ratio) == Approx(1.0));
  REQUIRE(grad1[0] == ValueApprox(all_grad[1][0]));
  REQUIRE(grad1[1] == ValueApprox(all_grad[1][1]));
  REQUIRE(grad1[2] == ValueApprox(all_grad[1][2]));

  // Compare ratio and ratioGrad with a non-zero displacement
  // Should compare more displacements
  ParticleSet::SingleParticlePos_t oldpos = elec->R[0];
  ParticleSet::SingleParticlePos_t displ(0.15, 0.10, 0.21);
  elec->R[0] = oldpos + displ;
  elec->update();
  ParticleSet::ParticleGradient_t new_grad;
  ParticleSet::ParticleLaplacian_t new_lap;
  new_grad.resize(nelec);
  new_lap.resize(nelec);

  // wavefunction value and derivatives at new position
  RealType new_logpsi = example_he->evaluateLog(*elec, new_grad, new_lap);
  elec->R[0]          = oldpos;
  elec->update();

  iat = 0;
  elec->makeMove(iat, displ);

  ratio = example_he->ratio(*elec, iat);
  REQUIRE(std::real(ratio) == Approx(std::exp(new_logpsi - logpsi)));

  ratio = example_he->ratioGrad(*elec, iat, grad0);

  REQUIRE(std::real(ratio) == Approx(std::exp(new_logpsi - logpsi)));

  REQUIRE(grad0[0] == ValueApprox(new_grad[0][0]));
  REQUIRE(grad0[1] == ValueApprox(new_grad[0][1]));
  REQUIRE(grad0[2] == ValueApprox(new_grad[0][2]));

  // Compare parameter derivatives

  const int nparam = 1;
  optimize::VariableSet var_param;
  example_he->checkInVariables(var_param);
  REQUIRE(var_param.size_of_active() == nparam);

  example_he->checkOutVariables(var_param);

  RealType old_B = example_he->B;
  // Interval size for finite-difference approximation to parameter derivative
  RealType h     = 0.01;
  RealType new_B = old_B + h;

  var_param["B"] = new_B;
  example_he->resetParameters(var_param);
  REQUIRE(example_he->B == Approx(new_B));

  ParticleSet::ParticleGradient_t grad_plus_h;
  ParticleSet::ParticleLaplacian_t lap_plus_h;
  grad_plus_h.resize(nelec);
  lap_plus_h.resize(nelec);

  RealType logpsi_plus_h = example_he->evaluateLog(*elec, grad_plus_h, lap_plus_h);


  // Finite difference derivative approximation
  RealType fd_logpsi = (logpsi_plus_h - logpsi) / h;

  std::vector<ValueType> dlogpsi(nparam);
  std::vector<ValueType> dhpsioverpsi(nparam);
  example_he->evaluateDerivatives(*elec, var_param, dlogpsi, dhpsioverpsi);

  REQUIRE(dlogpsi[0] == ValueApprox(fd_logpsi).epsilon(h));

  ValueType eloc   = -0.5 * (Sum(all_lap) + Dot(all_grad, all_grad));
  ValueType eloc_h = -0.5 * (Sum(lap_plus_h) + Dot(grad_plus_h, grad_plus_h));

  ValueType fd_eloc = (eloc_h - eloc) / h;

  REQUIRE(dhpsioverpsi[0] == ValueApprox(fd_eloc).epsilon(h));
}
} // namespace qmcplusplus
