//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"


#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/Tensor.h"
#include "Particle/ParticleSet.h"
#include "ParticleIO/LatticeIO.h"

#include <stdio.h>
#include <string>

using std::string;

namespace qmcplusplus
{
TEST_CASE("read_lattice_xml", "[particle_io][xml]")
{
  const char* particles = "<tmp> \
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

  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> uLattice;
  LatticeParser lp(uLattice);
  lp.put(root);


  REQUIRE(uLattice.R[0] == Approx(3.8));
  REQUIRE(uLattice.Volume == Approx(3.8 * 3.8 * 3.8));

  REQUIRE(uLattice.LR_dim_cutoff == Approx(20));
}
} // namespace qmcplusplus
