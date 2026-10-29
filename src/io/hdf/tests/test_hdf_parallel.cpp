
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
#include <catch2/catch_test_macros.hpp>

#include <vector>
#include "hdf/hdf_archive.h"
#include "Message/Communicate.h"

using namespace qmcplusplus;

// So far these tests only instantiate the types based on the type of MPI wrapper passed to the
//  constructor.
// Todo: Add tests that check parallel IO in the various configurations (ENABLE_PDHF5, etc)

TEST_CASE("hdf_archive_parallel", "[hdf]")
{
  Communicate& c(*OHMMS::Controller);
#if defined(ENABLE_PHDF5)
  {
    hdf_archive hd_write(c, true);
    REQUIRE(hd_write.create("test_communicate.hdf"));
    int i     = 23;
    bool okay = hd_write.writeEntry(i, "int");
    REQUIRE(okay);
  }
#else
  REQUIRE_THROWS_AS(hdf_archive(c, true), std::runtime_error);
  if (c.rank() == 0)
  {
    hdf_archive hd_write;
    REQUIRE(hd_write.create("test_communicate.hdf"));
    int i     = 23;
    bool okay = hd_write.writeEntry(i, "int");
    REQUIRE(okay);
  }
#endif
  c.barrier();

  hdf_archive hd_read(c, false);
  REQUIRE_THROWS_AS(hd_read.create("test_communicate_fail.hdf"), std::runtime_error);
  REQUIRE_THROWS_AS(hd_read.open("test_communicate.hdf", H5F_ACC_RDWR), std::runtime_error);
  REQUIRE(hd_read.open("test_communicate.hdf", H5F_ACC_RDONLY));
  int i_read = 0;
  REQUIRE(hd_read.readEntry(i_read, "int"));
  REQUIRE(i_read == 23);
  REQUIRE_THROWS_AS(hd_read.writeEntry(i_read, "int"), std::runtime_error);
}

#ifdef HAVE_MPI
TEST_CASE("hdf_archive_parallel_mpi3", "[hdf]")
{
  Communicate& c(*OHMMS::Controller);
  auto& comm = c.getCommMPI3();
#if defined(ENABLE_PHDF5)
  {
    hdf_archive hd_write(comm, true);
    REQUIRE(hd_write.create("test_mpi3_communicator.hdf"));
    int i     = 23;
    bool okay = hd_write.writeEntry(i, "int");
    REQUIRE(okay);
  }
#else
  REQUIRE_THROWS_AS(hdf_archive(comm, true), std::runtime_error);
  if (comm.rank() == 0)
  {
    hdf_archive hd_write;
    REQUIRE(hd_write.create("test_mpi3_communicator.hdf"));
    int i     = 23;
    bool okay = hd_write.writeEntry(i, "int");
    REQUIRE(okay);
  }
#endif
  comm.barrier();

  hdf_archive hd_read(comm, false);
  REQUIRE_THROWS_AS(hd_read.create("test_mpi3_communicator_fail.hdf"), std::runtime_error);
  REQUIRE_THROWS_AS(hd_read.open("test_mpi3_communicator.hdf", H5F_ACC_RDWR), std::runtime_error);
  REQUIRE(hd_read.open("test_mpi3_communicator.hdf", H5F_ACC_RDONLY));
  int i_read = 0;
  REQUIRE(hd_read.readEntry(i_read, "int"));
  REQUIRE(i_read == 23);
  REQUIRE_THROWS_AS(hd_read.writeEntry(i_read, "int"), std::runtime_error);
}
#endif
