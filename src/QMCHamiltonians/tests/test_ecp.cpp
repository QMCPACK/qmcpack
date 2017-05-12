//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "Utilities/OhmmsInfo.h"
#include "QMCHamiltonians/ECPComponentBuilder.h"


#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{

TEST_CASE("CheckSphericalIntegration", "[hamiltonian]")
{
  Communicate *c;
  OHMMS::Controller->initialize(0, NULL);
  c = OHMMS::Controller;
  OhmmsInfo("testlogfile");

  ECPComponentBuilder ecp("test_ecp",c);

  // Use the built-in quadrature rule check
  for (int quadrature_index = 1; quadrature_index < 8; quadrature_index++)
  {
    ecp.pp_nonloc = new NonLocalECPComponent;
    REQUIRE(ecp.SetQuadratureRule(quadrature_index));
  }
}

}

