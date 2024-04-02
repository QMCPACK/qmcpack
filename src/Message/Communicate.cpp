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

#ifdef HAVE_MPI
#include "mpi3/shared_communicator.hpp"
#endif


//Global Communicator is created without initialization
Communicate* OHMMS::Controller = new Communicate;

//default constructor: ready for a serial execution
Communicate::Communicate() : myMPI(MPI_COMM_NULL), d_mycontext(0), d_ncontexts(1), d_groupid(0), d_ngroups(1) {}

Communicate::~Communicate() = default;

//exclusive:  MPI or Serial
#ifdef HAVE_MPI

  // in_comm needs to be mutable to be duplicated
Communicate::Communicate(mpi3::communicator& in_comm) : d_groupid(0), d_ngroups(1), comm{in_comm.duplicate()}
{
  myMPI       = comm.get();
  d_mycontext = comm.rank();
  d_ncontexts = comm.size();
}

Communicate::Communicate(mpi3::communicator&& in_comm) : d_groupid(0), d_ngroups(1), comm{std::move(in_comm)}
{
  myMPI       = comm.get();
  d_mycontext = comm.rank();
  d_ncontexts = comm.size();
}

Communicate::Communicate(const Communicate& in_comm, int nparts)
{
  std::vector<int> nplist(nparts + 1);
  int p = FairDivideLow(in_comm.rank(), in_comm.size(), nparts, nplist); //group
  // comm is mutable member
  comm  = in_comm.comm.split(p, in_comm.rank());
  myMPI = comm.get();
  // TODO: mpi3 needs to define comm
  d_mycontext = comm.rank();
  d_ncontexts = comm.size();
  d_groupid   = p;
  d_ngroups   = nparts;
  // create a communicator among group leaders.

  nplist.pop_back();
  mpi3::communicator leader_comm = in_comm.comm.subcomm(nplist);
  if (isGroupLeader())
    GroupLeaderComm = std::make_unique<Communicate>(leader_comm);
  else
    GroupLeaderComm.reset();
}

// For unit tests until they can be changed and this will be removed.
void Communicate::initialize(int argc, char** argv) {}

Communicate Communicate::NodeComm() const {
  // comm is mutable member
  return Communicate{comm.split_shared()};
}

void Communicate::finalize()
{
  static bool has_finalized = false;

  if (!has_finalized)
  {
    has_finalized = true;
  }
}

void Communicate::cleanupMessage(void*) {}

void Communicate::abort() const { comm.abort(1); }

void Communicate::barrier() const { comm.barrier(); }
#else

void Communicate::initialize(int argc, char** argv) {}

Communicate Communicate::NodeComm() const {return Communicate{};}

void Communicate::finalize() {}

void Communicate::abort() const { std::_Exit(EXIT_FAILURE); }

void Communicate::barrier() const {}

void Communicate::cleanupMessage(void*) {}

Communicate::Communicate(const Communicate& in_comm, int nparts)
    : myMPI(MPI_COMM_NULL), d_mycontext(0), d_ncontexts(1), d_groupid(0)
{
  GroupLeaderComm = std::make_unique<Communicate>();
}
#endif // !HAVE_MPI

void Communicate::barrier_and_abort(const std::string& msg) const
{
  if (!rank())
    std::cerr << "Fatal Error. Aborting at " << msg << std::endl;
  Communicate::barrier();
  Communicate::abort();
}
