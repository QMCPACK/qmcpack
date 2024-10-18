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
  SECTION("valid p p p input")
  {
    const char* const particles = R"(
      <tmp>
        <parameter name="lattice" units="bohr">
          3.80000000       0.00000000       0.00000000
          0.00000000       3.80000000       0.00000000
          0.00000000       0.00000000       3.80000000
        </parameter>
        <parameter name="bconds">
          p p p
        </parameter>
        <parameter name="LR_dim_cutoff">20</parameter>
      </tmp>
    )";

    Libxml2Document doc;
    bool okay = doc.parseFromString(particles);
    REQUIRE(okay);

    xmlNodePtr root = doc.getRoot();

    CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> uLattice;
    LatticeParser lp(uLattice);
    REQUIRE_NOTHROW(lp.put(root));

    CHECK(uLattice.R[0] == Approx(3.8));
    CHECK(uLattice.Volume == Approx(3.8 * 3.8 * 3.8));

    CHECK(uLattice.LR_dim_cutoff == Approx(20));
  }

  SECTION("invalid n p p input")
  {
    const char* const particles = R"(
      <tmp>
        <parameter name="lattice" units="bohr">
          3.80000000       0.00000000       0.00000000
          0.00000000       3.80000000       0.00000000
          0.00000000       0.00000000       3.80000000
        </parameter>
        <parameter name="bconds">
          n p p
        </parameter>
        <parameter name="LR_dim_cutoff">20</parameter>
      </tmp>
    )";

    Libxml2Document doc;
    bool okay = doc.parseFromString(particles);
    REQUIRE(okay);

    xmlNodePtr root = doc.getRoot();

    CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> uLattice;
    LatticeParser lp(uLattice);
    REQUIRE_THROWS_WITH(lp.put(root),
                        "LatticeParser::put. In \"bconds\", non periodic directions must be placed after the periodic "
                        "ones.");
  }

  SECTION("invalid p n p input")
  {
    const char* const particles = R"(
      <tmp>
        <parameter name="lattice" units="bohr">
          3.80000000       0.00000000       0.00000000
          0.00000000       3.80000000       0.00000000
          0.00000000       0.00000000       3.80000000
        </parameter>
        <parameter name="bconds">
          p n p
        </parameter>
        <parameter name="LR_dim_cutoff">20</parameter>
      </tmp>
    )";

    Libxml2Document doc;
    bool okay = doc.parseFromString(particles);
    REQUIRE(okay);

    xmlNodePtr root = doc.getRoot();

    CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> uLattice;
    LatticeParser lp(uLattice);
    REQUIRE_THROWS_WITH(lp.put(root),
                        "LatticeParser::put. In \"bconds\", non periodic directions must be placed after the periodic "
                        "ones.");
  }
}

TEST_CASE("read_lattice_xml_lrhandle", "[particle_io][xml]")
{
  SECTION("valid p p n input")
  {
    const char* const particles = R"(
      <tmp>
        <parameter name="lattice" units="bohr">
          3.80000000       0.00000000       0.00000000
          0.00000000       3.80000000       0.00000000
          0.00000000       0.00000000      10.00000000
        </parameter>
        <parameter name="bconds">
          p p n
        </parameter>
        <parameter name="LR_dim_cutoff">30</parameter>
        <parameter name="LR_handler"> opt_breakup </parameter>
      </tmp>
    )";

    Libxml2Document doc;
    bool okay = doc.parseFromString(particles);
    REQUIRE(okay);

    xmlNodePtr root = doc.getRoot();

    CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> uLattice;
    LatticeParser lp(uLattice);
    REQUIRE_NOTHROW(lp.put(root));

    CHECK(uLattice.R[8] == Approx(10.0));
    CHECK(uLattice.LR_dim_cutoff == Approx(30));
  }

  SECTION("valid p p n ewald_quasi2d input")
  {
    const char* const particles = R"(
      <tmp>
        <parameter name="lattice" units="bohr">
          3.80000000       0.00000000       0.00000000
          0.00000000       3.80000000       0.00000000
          0.00000000       0.00000000      10.00000000
        </parameter>
        <parameter name="bconds">
          p p n
        </parameter>
        <parameter name="LR_dim_cutoff">30</parameter>
        <parameter name="LR_handler"> ewald_quasi2d </parameter>
      </tmp>
    )";

    Libxml2Document doc;
    bool okay = doc.parseFromString(particles);
    REQUIRE(okay);

    xmlNodePtr root = doc.getRoot();

    CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> uLattice;
    LatticeParser lp(uLattice);
    REQUIRE_NOTHROW(lp.put(root));
  }

  SECTION("invalid p p p ewald_quasi2d input")
  {
    const char* const particles = R"(
      <tmp>
        <parameter name="lattice" units="bohr">
          3.80000000       0.00000000       0.00000000
          0.00000000       3.80000000       0.00000000
          0.00000000       0.00000000      10.00000000
        </parameter>
        <parameter name="bconds">
          p p p
        </parameter>
        <parameter name="LR_dim_cutoff">30</parameter>
        <parameter name="LR_handler"> ewald_quasi2d </parameter>
      </tmp>
    )";

    Libxml2Document doc;
    bool okay = doc.parseFromString(particles);
    REQUIRE(okay);

    xmlNodePtr root = doc.getRoot();

    CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> uLattice;
    LatticeParser lp(uLattice);
    REQUIRE_THROWS_WITH(lp.put(root),
                        "LatticeParser::put. Quasi 2D Ewald only works with boundary condition 'p p n'!");
  }
}
} // namespace qmcplusplus
