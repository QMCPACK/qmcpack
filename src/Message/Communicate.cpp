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

Communicate::Communicate(const Communicate& in_comm, int nparts, int stripe)
{
  std::vector<int> nplist(nparts + 1);
  // group index
  const int gid =
      stripe == 0 ? FairDivideLow(in_comm.rank(), in_comm.size(), nparts, nplist) : in_comm.rank() / stripe % nparts;
  // comm is mutable member
  comm  = in_comm.comm.split(gid, in_comm.rank());
  myMPI = comm.get();
  // TODO: mpi3 needs to define comm
  d_mycontext = comm.rank();
  d_ncontexts = comm.size();
  d_groupid   = gid;
  d_ngroups   = nparts;

  // create an inter group communicator
  inter_group_comm_ = std::make_unique<Communicate>(in_comm.comm.split(comm.rank(), in_comm.rank()));
}

/** provide a node/shared-memory communicator from current (parent) communicator
 *
 *  The inter_group_comm_ is created across nodes for ranks sharing the same intra-node rank.
 *  Note: the size of the node communicator and its inter_group_comm_ may differ on different nodes.
 */
Communicate Communicate::NodeComm() const
{
  Communicate node_comm{comm.split_shared()};
  node_comm.inter_group_comm_ = std::make_unique<Communicate>(comm.split(node_comm.rank(), comm.rank()));
  return node_comm;
}

void Communicate::abort() const { comm.abort(1); }

void Communicate::barrier() const { comm.barrier(); }
#else

/** provide a node/shared-memory communicator from current (parent) communicator
 *
 *  Note: in non-MPI (serial) builds, this returns a single-rank communicator.
 */
Communicate Communicate::NodeComm() const
{
  Communicate node_comm;
  node_comm.inter_group_comm_ = std::make_unique<Communicate>();
  return node_comm;
}

void Communicate::abort() const { std::_Exit(EXIT_FAILURE); }

void Communicate::barrier() const {}


Communicate::Communicate(const Communicate& in_comm, int nparts, int stripe)
    : myMPI(MPI_COMM_NULL), d_mycontext(0), d_ncontexts(1), d_groupid(0), d_ngroups(nparts)
{ inter_group_comm_ = std::make_unique<Communicate>(); }
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
