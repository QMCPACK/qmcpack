//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "Configuration.h"
#include "Numerics/Quadrature.h"
#include "QMCHamiltonians/ECPComponentBuilder.h"

namespace qmcplusplus
{

TEST_CASE("CheckSphericalIntegration", "[hamiltonian]")
{
  // Use the built-in quadrature rule check
  for (int quadrature_index = 1; quadrature_index < 8; quadrature_index++)
  {
    Quadrature3D<QMCTraits::RealType> myRule(quadrature_index,false);
    REQUIRE(myRule.quad_ok);
  }
}

TEST_CASE("ReadFileBuffer_no_file","[hamiltonian]")
{
  ReadFileBuffer buf(NULL);
  bool open_okay = buf.open_file("does_not_exist");
  REQUIRE(open_okay == false);
}

TEST_CASE("ReadFileBuffer_simple_serial","[hamiltonian]")
{
  // Initializing with no Communicate pointer under MPI,
  //   this will read the file on every node.  Should be okay
  //   for testing purposes.
  ReadFileBuffer buf(NULL);
  bool open_okay = buf.open_file("simple.txt");
  REQUIRE(open_okay == true);

  bool read_okay = buf.read_contents();
  REQUIRE(read_okay);
  REQUIRE(buf.length == 14);
  REQUIRE(std::string("File contents\n") == buf.contents());
}

TEST_CASE("ReadFileBuffer_simple_mpi","[hamiltonian]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate *c = OHMMS::Controller;

  ReadFileBuffer buf(c);
  bool open_okay = buf.open_file("simple.txt");
  REQUIRE(open_okay == true);

  bool read_okay = buf.read_contents();
  REQUIRE(read_okay);
  REQUIRE(buf.length == 14);
  REQUIRE(std::string("File contents\n") == buf.contents());
}

TEST_CASE("ReadFileBuffer_ecp","[hamiltonian]")
{
  OHMMS::Controller->initialize(0, NULL);
  Communicate *c = OHMMS::Controller;

  ECPComponentBuilder ecp("test_read_ecp",c);

  bool okay = ecp.read_pp_file("C.BFD.xml");
  REQUIRE(okay);

  REQUIRE(ecp.Zeff == 4);

  // TODO: add more checks that pseudopotential file was read correctly
}

TEST_CASE("ReadFileBuffer_reopen","[hamiltonian]")
{
  // Initializing with no Communicate pointer under MPI,
  //   this will read the file on every node.  Should be okay
  //   for testing purposes.
  ReadFileBuffer buf(NULL);
  bool open_okay = buf.open_file("simple.txt");
  REQUIRE(open_okay == true);

  bool read_okay = buf.read_contents();
  REQUIRE(read_okay);
  REQUIRE(buf.length == 14);

  open_okay = buf.open_file("C.BFD.xml");
  REQUIRE(open_okay == true);

  read_okay = buf.read_contents();
  REQUIRE(read_okay);
  REQUIRE(buf.length > 14);
}

}
