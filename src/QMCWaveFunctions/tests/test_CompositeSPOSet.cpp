//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"
#include "QMCWaveFunctions/CompositeSPOSet.h"
#include <exception>
#include "Particle/tests/MinimalParticlePool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"

namespace qmcplusplus
{

TEST_CASE("CompositeSPO::diamond_1x1x1", "[wavefunction")
{
  Libxml2Document doc;

  Communicate* comm;
  comm = OHMMS::Controller;
  outputManager.pause();

  auto particle_pool     = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto wavefunction_pool = MinimalWaveFunctionPool::make_diamondC_1x1x1(comm, particle_pool);
  auto& pset             = *particle_pool.getParticleSet("e");
  auto& twf              = *wavefunction_pool.getWaveFunction("wavefunction");

  CompositeSPOSet comp_sposet;

  std::vector<std::string> sposets{"spo_ud", "spo_dm"};
  for (auto sposet_str : sposets)
  {
    auto& sposet = twf.getSPOSet(sposet_str);
    comp_sposet.add(sposet.makeClone());
  }
  CHECK(comp_sposet.size() == 8);

  SPOSet::ValueMatrix psiM(pset.R.size(), comp_sposet.getOrbitalSetSize());
  SPOSet::GradMatrix dpsiM(pset.R.size(), comp_sposet.getOrbitalSetSize());
  SPOSet::ValueMatrix d2psiM(pset.R.size(), comp_sposet.getOrbitalSetSize());
  comp_sposet.evaluate_notranspose(pset, 0, pset.R.size(), psiM, dpsiM, d2psiM);
}
} // namespace qmcplusplus
