
//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include <vector>
#include "hdf/hdf_archive.h"
#include "Message/Communicate.h"

using namespace qmcplusplus;

// So far these tests only instantiate the types based on the type of MPI wrapper passed to the
//  constructor.
// Todo: Add tests that check parallel IO in the various configurations (ENABLE_PDHF5, etc)

TEST_CASE("hdf_archive_parallel", "[hdf]")
{
  Communicate* c = OHMMS::Controller;
  hdf_archive hd(c);
  hd.create("test_communicate.hdf");
  int i     = 23;
  bool okay = hd.writeEntry(i, "int");
  REQUIRE(okay);
}

#ifdef HAVE_MPI
TEST_CASE("hdf_archive_parallel_mpi3", "[hdf]")
{
  Communicate* c = OHMMS::Controller;
  hdf_archive hd(c->comm);
  hd.create("test_mpi3_communicator.hdf");
  int i     = 23;
  bool okay = hd.writeEntry(i, "int");
  REQUIRE(okay);
}
#endif
