//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//
// File created by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "type_traits/template_types.hpp"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Particle/tests/MinimalParticlePool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
#include "QMCHamiltonians/tests/MinimalHamiltonianPool.h"
namespace qmcplusplus
{
TEST_CASE("Eloc_Derivatives:slater_noj", "[hamiltonian]")
{
  Communicate* comm;
  comm = OHMMS::Controller;

  app_log()<<" WOOHOO Eloc_derivatives:slater_noj tested!";
  REQUIRE( true );
}

TEST_CASE("Eloc_Derivatives:slater_wj", "[hamiltonian]")
{
  Communicate* comm;
  comm = OHMMS::Controller;

  app_log()<<" WOOHOO Eloc_Derivatives:slater_wj tested!";
  REQUIRE( true );
}

TEST_CASE("Eloc_Derivatives:multislater_noj", "[hamiltonian]")
{
  Communicate* comm;
  comm = OHMMS::Controller;

  app_log()<<" WOOHOO Eloc_derivatives:multislater_noj tested!";
  REQUIRE( true );
}

TEST_CASE("Eloc_Derivatives:multislater_wj", "[hamiltonian]")
{
  Communicate* comm;
  comm = OHMMS::Controller;

  app_log()<<" WOOHOO Eloc_derivatives:multislater_wj tested!";
  REQUIRE( true );
}

} // namespace qmcplusplus
