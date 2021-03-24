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

  bool evaluateDerivatives(RealType r, std::vector<TinyVector<RealType, 3>>& derivs) { return true; }

  void resetParameters(const opt_variables_type& active) {}
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
  SpeciesSet& tspecies = elec.getSpeciesSet();
  int upIdx            = tspecies.addSpecies("u");
  int downIdx          = tspecies.addSpecies("d");
  elec.resetGroups();
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
  //jorb.getVars().print(std::cout,0,true);

  CHECK(global_active.size_of_active() == 2);

  // Order is based on the function list F
  // For two species (ia - index of first species, ib - index of second species)
  // F[0] is (0,0)
  // F[1] and F[2] are (0,1),(1,0) - use the same functions
  // F[3] is (1,1) b-b (by default uses the same function as a-a)

  // Index into global_active
  auto o1 = jorb.getOffset(0);
  CHECK(o1.first == 0);
  CHECK(o1.second == 1);

  auto o2 = jorb.getOffset(1);
  CHECK(o2.first == 1);
  CHECK(o2.second == 2);

  auto o3 = jorb.getOffset(2);
  CHECK(o3.first == 1);
  CHECK(o3.second == 2);

  auto o4 = jorb.getOffset(3);
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
  auto o1 = jorb.getOffset(0);
  CHECK(o1.first == -1);

  // Offset into set of active variables (global_active)
  auto o2 = jorb.getOffset(1);
  CHECK(o2.first == 0);
  CHECK(o2.second == 1);

  auto o3 = jorb.getOffset(2);
  CHECK(o3.first == 0);
  CHECK(o3.second == 1);

  auto o4 = jorb.getOffset(3);
  CHECK(o4.first == -1);
}

} // namespace qmcplusplus
