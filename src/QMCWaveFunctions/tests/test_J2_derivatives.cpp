//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include <memory>
#include "Jastrow/J2OrbitalSoA.h"
#include "Jastrow/FakeFunctor.h"

namespace qmcplusplus
{

using FakeJasFunctor = FakeFunctor<OHMMS_PRECISION>;

TEST_CASE("J2OrbitalSoA simple", "[wavefunction]")
{
  const SimulationCell simulation_cell;
  ParticleSet elec(simulation_cell);
  elec.setName("e");
  elec.create({1,1});
  J2OrbitalSoA<FakeJasFunctor> jorb("J2_fake", elec);

  opt_variables_type active;
  jorb.checkOutVariables(active);
}

TEST_CASE("J2OrbitalSoA one species and two variables", "[wavefunction]")
{
  const SimulationCell simulation_cell;
  ParticleSet elec(simulation_cell);
  elec.setName("e");
  elec.create({1,1});
  J2OrbitalSoA<FakeJasFunctor> jorb("J2_fake", elec);

  auto j2_uptr = std::make_unique<FakeJasFunctor>();
  auto& j2     = *j2_uptr;
  j2.myVars.insert("opt1", 1.0);
  j2.myVars.insert("opt2", 2.0);
  // update num_active_vars
  j2.myVars.resetIndex();
  jorb.addFunc(0, 0, std::move(j2_uptr));

  opt_variables_type global_active;
  global_active.insertFrom(j2.myVars);

  jorb.checkOutVariables(global_active);

  CHECK(global_active.size_of_active() == 2);
}

ParticleSet get_two_species_particleset(const SimulationCell& simulation_cell)
{
  ParticleSet elec(simulation_cell);
  std::vector<int> ud{2, 2};
  elec.setName("e");
  elec.create(ud);

  elec.R[0] = {1.0, 0.0, 0.0};
  elec.R[1] = {1.1, 1.0, 0.1};
  elec.R[2] = {0.9, 0.8, 1.0};
  elec.R[3] = {0.9, 0.5, 1.1};

  SpeciesSet& tspecies = elec.getSpeciesSet();
  int upIdx            = tspecies.addSpecies("u");
  int downIdx          = tspecies.addSpecies("d");
  elec.resetGroups();

  elec.addTable(elec);
  elec.update();

  return elec;
}

// Two variables, both active
TEST_CASE("J2OrbitalSoA two variables", "[wavefunction]")
{
  const SimulationCell simulation_cell;
  ParticleSet elec = get_two_species_particleset(simulation_cell);

  J2OrbitalSoA<FakeJasFunctor> jorb("J2_fake", elec);

  auto j2a_uptr = std::make_unique<FakeJasFunctor>();
  auto& j2a     = *j2a_uptr;
  j2a.myVars.insert("opt1", 1.0);
  // update num_active_vars
  j2a.myVars.resetIndex();
  jorb.addFunc(0, 0, std::move(j2a_uptr));

  auto j2b_uptr = std::make_unique<FakeJasFunctor>();
  auto& j2b     = *j2b_uptr;
  j2b.myVars.insert("opt2", 2.0);
  // update num_active_vars
  j2b.myVars.resetIndex();
  jorb.addFunc(0, 1, std::move(j2b_uptr));

  opt_variables_type global_active;
  global_active.insertFrom(j2a.myVars);
  global_active.insertFrom(j2b.myVars);
  global_active.resetIndex();

  jorb.checkOutVariables(global_active);

  // Formatted output of variables involved
  //j2a.myVars.print(std::cout,0,true);
  //j2b.myVars.print(std::cout,0,true);
  //global_active.print(std::cout,0,true);
  //jorb.getComponentVars().print(std::cout,0,true);

  CHECK(global_active.size_of_active() == 2);

  // Order is based on the function list F
  // For two species (ia - index of first species, ib - index of second species)
  // F[0] is (0,0)
  // F[1] and F[2] are (0,1),(1,0) - use the same functions
  // F[3] is (1,1) b-b (by default uses the same function as a-a)

  // Index into global_active
  auto o1 = jorb.getComponentOffset(0);
  CHECK(o1.first == 0);
  CHECK(o1.second == 1);

  auto o2 = jorb.getComponentOffset(1);
  CHECK(o2.first == 1);
  CHECK(o2.second == 2);

  auto o3 = jorb.getComponentOffset(2);
  CHECK(o3.first == 1);
  CHECK(o3.second == 2);

  auto o4 = jorb.getComponentOffset(3);
  CHECK(o4.first == 0);
  CHECK(o4.second == 1);
}

// Reproduce 2814.  If the variational parameter for the first Jastrow factor
// is not active, the offsets are not correct.
// "First" means the first in the function list F, which has species indices 0,0
TEST_CASE("J2OrbitalSoA variables fail", "[wavefunction]")
{
  const SimulationCell simulation_cell;
  ParticleSet elec = get_two_species_particleset(simulation_cell);

  J2OrbitalSoA<FakeJasFunctor> jorb("J2_fake", elec);

  auto j2a_uptr = std::make_unique<FakeJasFunctor>();
  auto& j2a     = *j2a_uptr;
  j2a.myVars.insert("opt1", 1.0);
  // update num_active_vars
  j2a.myVars.resetIndex();
  jorb.addFunc(0, 0, std::move(j2a_uptr));

  auto j2b_uptr = std::make_unique<FakeJasFunctor>();
  auto& j2b     = *j2b_uptr;
  j2b.myVars.insert("opt2", 2.0);
  // update num_active_vars
  j2b.myVars.resetIndex();
  jorb.addFunc(0, 1, std::move(j2b_uptr));

  opt_variables_type global_active;
  global_active.insertFrom(j2b.myVars);
  global_active.resetIndex();


  jorb.checkOutVariables(global_active);

  CHECK(global_active.size_of_active() == 1);
  // Not optimizing the parameter in this Jastrow factor, indicated by first index is -1
  auto o1 = jorb.getComponentOffset(0);
  CHECK(o1.first == -1);

  // Offset into set of active variables (global_active)
  auto o2 = jorb.getComponentOffset(1);
  CHECK(o2.first == 0);
  CHECK(o2.second == 1);

  auto o3 = jorb.getComponentOffset(2);
  CHECK(o3.first == 0);
  CHECK(o3.second == 1);

  auto o4 = jorb.getComponentOffset(3);
  CHECK(o4.first == -1);

  using ValueType = QMCTraits::ValueType;

  // Check derivative indexing
  int num_vars = 1;
  j2b.derivs_.resize(num_vars);
  // Elements are d/dp_i u(r), d/dp_i du/dr,  d/dp_i d2u/dr2
  j2b.derivs_[0] = {0.5, 1.3, 2.4};
  std::vector<ValueType> dlogpsi(num_vars);
  jorb.evaluateDerivativesWF(elec, global_active, dlogpsi);

  CHECK(dlogpsi[0] == ValueApprox(-2.0)); // 4 * derivs_[0][0]

  std::vector<ValueType> dlogpsi2(num_vars);
  std::vector<ValueType> dhpsioverpsi(num_vars);

  jorb.evaluateDerivatives(elec, global_active, dlogpsi2, dhpsioverpsi);
  CHECK(dlogpsi2[0] == ValueApprox(-2.0));
}

// Other variational parameters in the wavefunction (e.g. one-body Jastrow)

TEST_CASE("J2OrbitalSoA other variables", "[wavefunction]")
{
  const SimulationCell simulation_cell;
  ParticleSet elec = get_two_species_particleset(simulation_cell);

  J2OrbitalSoA<FakeJasFunctor> jorb("J2_fake", elec);

  auto j2a_uptr = std::make_unique<FakeJasFunctor>();
  auto j2a      = *j2a_uptr;
  j2a.myVars.insert("opt1", 1.0);
  // update num_active_vars
  j2a.myVars.resetIndex();
  jorb.addFunc(0, 0, std::move(j2a_uptr));

  auto j2b_uptr = std::make_unique<FakeJasFunctor>();
  auto& j2b     = *j2b_uptr;
  j2b.myVars.insert("opt2", 2.0);
  // update num_active_vars
  j2b.myVars.resetIndex();
  jorb.addFunc(0, 1, std::move(j2b_uptr));

  opt_variables_type global_active;
  // This is a parameter from another part of the wavefunction
  global_active.insert("other_opt", 1.0);
  global_active.insertFrom(j2b.myVars);
  global_active.resetIndex();


  jorb.checkOutVariables(global_active);

  //global_active.print(std::cout,0,true);
  //jorb.getComponentVars().print(std::cout,0,true);

  CHECK(global_active.size_of_active() == 2);
  // Not optimizing the parameter in this Jastrow factor, indicated by first index is -1
  auto o1 = jorb.getComponentOffset(0);
  CHECK(o1.first < 0);

  // Offset into set of active variables (global_active)
  auto o2 = jorb.getComponentOffset(1);
  CHECK(o2.first == 0);
  CHECK(o2.second == 1);

  auto o3 = jorb.getComponentOffset(2);
  CHECK(o3.first == 0);
  CHECK(o3.second == 1);

  auto o4 = jorb.getComponentOffset(3);
  CHECK(o4.first < 0);


  using ValueType = QMCTraits::ValueType;

  // Check derivative indexing
  int num_vars = 2;
  j2b.derivs_.resize(num_vars);
  // Elements are d/dp_i u(r), d/dp_i du/dr,  d/dp_i d2u/dr2
  j2b.derivs_[0] = {0.5, 1.3, 2.4};
  std::vector<ValueType> dlogpsi(num_vars);
  jorb.evaluateDerivativesWF(elec, global_active, dlogpsi);

  CHECK(dlogpsi[1] == ValueApprox(-2.0)); // 4 * derivs_[0][0]

  std::vector<ValueType> dlogpsi2(num_vars);
  std::vector<ValueType> dhpsioverpsi(num_vars);

  jorb.evaluateDerivatives(elec, global_active, dlogpsi2, dhpsioverpsi);
  CHECK(dlogpsi2[1] == ValueApprox(-2.0));
}

// Reproduce 3137.  If the number of particle types equals the number of particles
// the two body jastrow is not constructed correctly (except in the case of two
// particles).
TEST_CASE("J2OrbitalSoA Jastrow three particles of three types", "[wavefunction]")
{
  const SimulationCell simulation_cell;
  ParticleSet ions(simulation_cell);
  ParticleSet elec(simulation_cell);

  ions.setName("ion");
  ions.create(1);
  ions.R[0][0] = 0.0;
  ions.R[0][1] = 0.0;
  ions.R[0][2] = 0.0;

  elec.setName("elec");
  std::vector<int> udp(3);
  udp[0] = 1;
  udp[1] = 1;
  udp[2] = 1;
  elec.create(udp);
  elec.R[0][0] = -0.28;
  elec.R[0][1] = 0.0225;
  elec.R[0][2] = -2.709;
  elec.R[1][0] = -1.08389;
  elec.R[1][1] = 1.9679;
  elec.R[1][2] = -0.0128914;
  elec.R[2][0] = -2.08389;
  elec.R[2][1] = 0.9679;
  elec.R[2][2] = 0.0128914;


  J2OrbitalSoA<FakeJasFunctor> jorb("J2_fake", elec);

  // 0 uu  (0,0)
  // 1 ud  (0,1)
  // 2 up  (0,2)
  // 3 du  (1,0)
  // 4 dd  (1,1)
  // 5 dp  (1,2)
  // 6 pu  (2,0)
  // 7 pd  (2,1)
  // 8 pp  (2,2)

  auto j2a_uptr = std::make_unique<FakeJasFunctor>();
  auto& j2a     = *j2a_uptr;
  j2a.myVars.insert("opt1", 1.0);
  // update num_active_vars
  j2a.myVars.resetIndex();
  jorb.addFunc(0, 1, std::move(j2a_uptr));

  auto j2b_uptr = std::make_unique<FakeJasFunctor>();
  auto j2b      = *j2b_uptr;
  j2b.myVars.insert("opt2", 2.0);
  // update num_active_vars
  j2b.myVars.resetIndex();
  jorb.addFunc(0, 2, std::make_unique<FakeJasFunctor>(*j2b_uptr));

  // currently opposite spins won't be set to be equivalent
  // setting u,p doesn't set d,p
  jorb.addFunc(1, 2, std::move(j2b_uptr));

  auto& F = jorb.getPairFunctions();
  for (size_t i = 0; i < F.size(); ++i)
    CHECK(F[i] != nullptr);
}

} // namespace qmcplusplus
