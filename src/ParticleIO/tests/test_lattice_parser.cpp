
#include "catch.hpp"


#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/Tensor.h"
#include "Utilities/OhmmsInfo.h"
#include "Lattice/Uniform3DGridLayout.h"
#include "Particle/ParticleSet.h"
#include "ParticleIO/ParticleLayoutIO.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
TEST_CASE("read_lattice_xml", "[particle_io][xml]")
{

  OHMMS::Controller->initialize(0, NULL);
  Communicate *c = OHMMS::Controller;
  OhmmsInfo("testlogfile");


const char *particles = \
"<tmp> \
 <parameter name=\"lattice\" units=\"bohr\"> \
                 3.80000000       0.00000000       0.00000000 \
                 0.00000000       3.80000000       0.00000000 \
                 0.00000000       0.00000000       3.80000000 \
         </parameter> \
         <parameter name=\"bconds\"> \
            p p p \
         </parameter> \
         <parameter name=\"LR_dim_cutoff\">20</parameter> \
</tmp> \
";
  Libxml2Document doc;
  bool okay = doc.parseFromString(particles);
  REQUIRE(okay);

  xmlNodePtr root = doc.getRoot();

  Uniform3DGridLayout ugrid;
  LatticeParser lp(ugrid);
  lp.put(root);

 
  REQUIRE(ugrid.R[0] == Approx(3.8));
  REQUIRE(ugrid.Volume == Approx(3.8*3.8*3.8));

  REQUIRE(ugrid.LR_dim_cutoff == Approx(20));

}
}
