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

#include "QMCWaveFunctions/Jastrow/DiffTwoBodyJastrowOrbital.h"

namespace qmcplusplus
{
// Mock radial function to use in DiffTwoBodyJastrowOrbital
class FakeJastrow
{
public:
  std::string name;

  opt_variables_type myVars;

  using RealType = QMCTraits::RealType;

  bool evaluateDerivatives(RealType r, std::vector<TinyVector<RealType, 3>>& derivs)
  {
    derivs = derivs_;
    return true;
  }

  void resetParameters(const opt_variables_type& active) {}

  std::vector<TinyVector<RealType, 3>> derivs_;
};

TEST_CASE("DiffTwoBodyJastrowOrbital simple", "[wavefunction]")
{
  ParticleSet elec;
  elec.setName("e");
  DiffTwoBodyJastrowOrbital<FakeJastrow> jorb(elec);

  opt_variables_type active;
  jorb.checkOutVariables(active);
}

TEST_CASE("DiffTwoBodyJastrowOrbital one species and two variables", "[wavefunction]")
{
  ParticleSet elec;
  elec.setName("e");
  DiffTwoBodyJastrowOrbital<FakeJastrow> jorb(elec);

  FakeJastrow j2;
  j2.myVars.insert("opt1", 1.0);
  j2.myVars.insert("opt2", 2.0);
  // update num_active_vars
  j2.myVars.resetIndex();
  jorb.addFunc(0, 0, &j2);

  opt_variables_type global_active;
  global_active.insertFrom(j2.myVars);

  jorb.checkOutVariables(global_active);

  CHECK(global_active.size_of_active() == 2);
}

ParticleSet get_two_species_particleset()
{
  ParticleSet elec;
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
TEST_CASE("DiffTwoBodyJastrowOrbital two variables", "[wavefunction]")
{
  ParticleSet elec = get_two_species_particleset();

  DiffTwoBodyJastrowOrbital<FakeJastrow> jorb(elec);

  FakeJastrow j2a;
  j2a.myVars.insert("opt1", 1.0);
  j2a.name = "j2a";
  // update num_active_vars
  j2a.myVars.resetIndex();
  jorb.addFunc(0, 0, &j2a);

  FakeJastrow j2b;
  j2b.myVars.insert("opt2", 2.0);
  j2b.name = "j2b";
  // update num_active_vars
  j2b.myVars.resetIndex();
  jorb.addFunc(0, 1, &j2b);

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
TEST_CASE("DiffTwoBodyJastrowOrbital variables fail", "[wavefunction]")
{
  ParticleSet elec = get_two_species_particleset();

  DiffTwoBodyJastrowOrbital<FakeJastrow> jorb(elec);

  FakeJastrow j2a;
  j2a.myVars.insert("opt1", 1.0);
  j2a.name = "j2a";
  // update num_active_vars
  j2a.myVars.resetIndex();
  jorb.addFunc(0, 0, &j2a);

  FakeJastrow j2b;
  j2b.myVars.insert("opt2", 2.0);
  j2b.name = "j2b";
  // update num_active_vars
  j2b.myVars.resetIndex();
  jorb.addFunc(0, 1, &j2b);

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

TEST_CASE("DiffTwoBodyJastrowOrbital other variables", "[wavefunction]")
{
  ParticleSet elec = get_two_species_particleset();

  DiffTwoBodyJastrowOrbital<FakeJastrow> jorb(elec);

  FakeJastrow j2a;
  j2a.myVars.insert("opt1", 1.0);
  j2a.name = "j2a";
  // update num_active_vars
  j2a.myVars.resetIndex();
  jorb.addFunc(0, 0, &j2a);

  FakeJastrow j2b;
  j2b.myVars.insert("opt2", 2.0);
  j2b.name = "j2b";
  // update num_active_vars
  j2b.myVars.resetIndex();
  jorb.addFunc(0, 1, &j2b);

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

} // namespace qmcplusplus
