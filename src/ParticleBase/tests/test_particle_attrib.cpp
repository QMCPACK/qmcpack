//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign 
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#include "Message/catch_mpi_main.hpp"

#include "Message/Communicate.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OhmmsPETE/TinyVector.h"
#include "Lattice/CrystalLattice.h"
#include "ParticleBase/ParticleAttrib.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{

TEST_CASE("particle_attrib_scalar", "[particle_base]")
{
  OHMMS::Controller->initialize(0, NULL);

  ParticleAttrib<double> PA1;
  REQUIRE(PA1.size() == 0);

  PA1.resize(4);
  REQUIRE(PA1.size() == 4);

  ParticleAttrib<double> PA3(3);
  REQUIRE(PA3.size() == 3);

  REQUIRE(PA3[0] == 0.0);
  // Whole array operation
  PA3 = 1.0;
  REQUIRE(PA3[0] == 1.0);
  REQUIRE(PA3[1] == 1.0);
  REQUIRE(PA3[2] == 1.0);

  // Single element
  PA3[0] = 2.0;
  REQUIRE(PA3[0] == 2.0);
  REQUIRE(PA3[1] == 1.0);
}

TEST_CASE("particle_attrib_vector", "[particle_base]")
{
  OHMMS::Controller->initialize(0, NULL);

  ParticleAttrib<TinyVector<double, 2> > PA1;
  REQUIRE(PA1.size() == 0);

  PA1.resize(3);
  REQUIRE(PA1.size() == 3);

  PA1[0] = 1.0;
  PA1[1] = 0.0;
  REQUIRE(PA1[0][0] == 1.0);
  REQUIRE(PA1[0][1] == 1.0);
  REQUIRE(PA1[1][0] == 0.0);
}

}
