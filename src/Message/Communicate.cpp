//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Anouar Benali, benali@anl.gov, Argonne National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//                    Alfredo A. Correa, correaa@llnl.gov, Lawrence Livermore National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "Communicate.h"
#include <iostream>
#include <sstream>
#include <string>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <numeric>
#include <cassert>
#include "config.h"
#include "Utilities/FairDivide.h"

#ifdef ENABLE_GCOV
#ifdef __GNUC__
extern "C" void __gcov_dump();
#endif
#else
#endif

#ifdef HAVE_MPI
#include "mpi3/shared_communicator.hpp"
#endif

//Global Communicator is created without initialization
std::unique_ptr<Communicate> OHMMS::Controller = std::make_unique<Communicate>();

//default constructor: ready for a serial execution
Communicate::Communicate() : myMPI(MPI_COMM_NULL), d_mycontext(0), d_ncontexts(1), d_groupid(0), d_ngroups(1) {}

Communicate::~Communicate() = default;
Communicate::Communicate(Communicate&&) = default;

//exclusive:  MPI or Serial
#ifdef HAVE_MPI


// Takes ownership of a communicator (e.g. from split or split_shared) via move semantics,
// avoiding the collective overhead and extra context allocation of MPI_Comm_dup.
Communicate::Communicate(mpi3::communicator&& in_comm) : d_groupid(0), d_ngroups(1), comm{std::move(in_comm)}
{
  myMPI       = comm.get();
  d_mycontext = comm.rank();
  d_ncontexts = comm.size();
}

Communicate Communicate::splitByStride(int nparts, int stride) const
{
  assert(stride >= 1);
  assert(size() % stride == 0);

  // group index
  const int gid = rank() / stride % nparts;
  Communicate comm_split(comm.split(gid, rank()));
  comm_split.d_groupid = gid;
  comm_split.d_ngroups = nparts;

  // create an inter group communicator
  comm_split.inter_group_comm_ = std::make_unique<Communicate>(comm.split(comm_split.comm.rank(), rank()));
  return comm_split;
}

Communicate Communicate::split(int nparts, int multiple) const
{
  assert(multiple >= 1);
  assert(size() % multiple == 0);
  std::vector<int> nplist(nparts + 1);
  const int gid = FairDivideLow(rank() / multiple, size() / multiple, nparts, nplist);
  Communicate comm_split(comm.split(gid, rank()));
  comm_split.d_groupid = gid;
  comm_split.d_ngroups = nparts;

  comm_split.inter_group_comm_ = std::make_unique<Communicate>(comm.split(comm_split.comm.rank(), rank()));
  return comm_split;
}

/** provide a node/shared-memory communicator from current (parent) communicator
 *
 *  Note: the size of the node communicator may differ on different nodes.
 */
Communicate Communicate::createNodeComm() const
{
  return Communicate{comm.split_shared()};
}

int Communicate::getUniformNodeRankCount() const
{
  Communicate node_comm = createNodeComm();
  const int local_node_size = node_comm.size();
  std::vector<int> all_node_sizes(size());
  MPI_Allgather(&local_node_size, 1, MPI_INT, all_node_sizes.data(), 1, MPI_INT, myMPI);
  int gcd_size = all_node_sizes.empty() ? 1 : all_node_sizes.front();
  for (int s : all_node_sizes)
    gcd_size = std::gcd(gcd_size, s);

  return gcd_size;
}


void Communicate::abort() const { comm.abort(1); }

void Communicate::barrier() const { comm.barrier(); }
#else

/** provide a node/shared-memory communicator from current (parent) communicator
 *
 *  Note: in non-MPI (serial) builds, this returns a single-rank communicator.
 */
Communicate Communicate::createNodeComm() const
{
  return Communicate{};
}

int Communicate::getUniformNodeRankCount() const
{
  return 1;
}


void Communicate::abort() const { std::_Exit(EXIT_FAILURE); }

void Communicate::barrier() const {}


Communicate Communicate::splitByStride(int nparts, int stride) const
{
  assert(stride >= 1);
  assert(size() % stride == 0);

  Communicate comm_split;
  comm_split.d_ngroups = nparts;
  comm_split.inter_group_comm_ = std::make_unique<Communicate>();
  return comm_split;
}

Communicate Communicate::split(int nparts, int multiple) const
{
  assert(multiple >= 1);
  assert(size() % multiple == 0);
  Communicate comm_split;
  comm_split.d_ngroups = nparts;
  comm_split.inter_group_comm_ = std::make_unique<Communicate>();
  return comm_split;
}
#endif // !HAVE_MPI

void Communicate::barrier_and_abort(const std::string& msg) const
{
  if (!rank())
    std::cerr << "Fatal Error. Aborting at " << msg << std::endl;
  Communicate::barrier();
#ifdef ENABLE_GCOV
#ifdef __GNUC__
  __gcov_dump();
#endif
#endif
  Communicate::abort();
}

Communicate& Communicate::getInterGroupComm() const
{
  if (!inter_group_comm_)
    throw std::runtime_error("Communicate::getInterGroupComm(): inter_group_comm_ is not initialized.");
  return *inter_group_comm_;
}
