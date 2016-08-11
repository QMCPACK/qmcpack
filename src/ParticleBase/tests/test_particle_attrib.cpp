
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "Message/Communicate.h"
#include "Utilities/OhmmsInfo.h"
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
  OhmmsInfo("testlogfile");

  ParticleAttrib<double> PA1;
  REQUIRE(PA1.size() == 0);

  PA1.create(4);
  REQUIRE(PA1.size() == 4);

  ParticleAttrib<double> PA2("type_name", "object_name");
  REQUIRE(PA2.size() == 0);

  ParticleAttrib<double> PA3("type_name", "object_name",3);
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
  OhmmsInfo("testlogfile");

  ParticleAttrib<TinyVector<double, 2> > PA1;
  REQUIRE(PA1.size() == 0);

  PA1.create(3);
  REQUIRE(PA1.size() == 3);

  PA1[0] = 1.0;
  REQUIRE(PA1[0][0] == 1.0);
  REQUIRE(PA1[0][1] == 1.0);
  REQUIRE(PA1[1][0] == 0.0);
}

}
