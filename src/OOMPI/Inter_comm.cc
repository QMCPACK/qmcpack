// -*- c++ -*-
//
// Copyright (c) 2002-2003 Indiana University.  All rights reserved.
// Copyright (c) 1996, 1997, 1998, 2000 University of Notre Dame.
//                         All rights reserved.
// 
// This file is part of the OOMPI software package.  For license
// information, see the LICENSE file in the top level directory of the
// OOMPI source distribution.
//
// $Id$
//
// OOMPI Class library
// Inter_communicators
//

#include "Inter_comm.h"


// Creates an OOMPI_Inter_comm from an mpi_comm
void OOMPI_Inter_comm::do_full_init(MPI_Comm mpi_comm, bool needs_to_be_freed)
{
  if (mpi_comm == MPI_COMM_NULL) 
    return;

  int flag(0);
  MPI_Comm_test_inter(mpi_comm, &flag);
  
  if (flag != 0)
    MPI_Constructor(mpi_comm, needs_to_be_freed);
}


OOMPI_Inter_comm::OOMPI_Inter_comm(const OOMPI_Inter_comm &a)
:OOMPI_Comm(a)
{
}

OOMPI_Inter_comm::~OOMPI_Inter_comm()
{
}


OOMPI_Inter_comm &
OOMPI_Inter_comm::operator=(const OOMPI_Inter_comm &a)
{
  (OOMPI_Comm &) *this = (OOMPI_Comm &) a; 
  return *this;
}


// Does an MPI_Intercomm_create 
OOMPI_Inter_comm::OOMPI_Inter_comm(OOMPI_Intra_comm &local_comm, 
				   int local_leader, 
				   OOMPI_Intra_comm &peer_comm, 
				   int remote_leader, int tag)
{
  MPI_Comm mpi_comm;

  MPI_Intercomm_create(local_comm.comm_wrapper->Get(), local_leader, 
		       peer_comm.comm_wrapper->Get(), remote_leader, 
		       tag, &mpi_comm);


  MPI_Constructor(mpi_comm, true);

}


// Does an MPI_Comm_dup and returns the new inter comm
OOMPI_Inter_comm 
OOMPI_Inter_comm::Dup(void)
{
  MPI_Comm mpi_comm;
  if (MPI_Comm_dup(Get_mpi(), &mpi_comm) != MPI_SUCCESS)
    return OOMPI_Inter_comm(MPI_COMM_NULL);
 
  return OOMPI_Inter_comm(mpi_comm, true);
}



// ----------- Type - virtual function from comm
// returns flag  

bool
OOMPI_Inter_comm::Test_inter()
{
  int flag(0);

  MPI_Comm_test_inter(comm_wrapper->Get(), &flag);

  return (bool) flag;
}


// returns an OOMPI_Group that holds the mpi_group
OOMPI_Group
OOMPI_Inter_comm::Remote_group()
{
  MPI_Group mpi_group;

  if (MPI_Comm_remote_group(comm_wrapper->Get(), &mpi_group) != MPI_SUCCESS)
    return OOMPI_Group(MPI_GROUP_NULL);
  
  return OOMPI_Group(mpi_group);
}
 

// returns size
int 
OOMPI_Inter_comm::Remote_size()
{
  int size;
  
  if (MPI_Comm_remote_size(comm_wrapper->Get(), &size) != MPI_SUCCESS)
    return -1;

  return size;
}


// returns mpi_comm
OOMPI_Intra_comm
OOMPI_Inter_comm::Merge(bool high)
{
  MPI_Comm mpi_comm;
  if (MPI_Intercomm_merge(comm_wrapper->Get(), (int) high, &mpi_comm) != 
      MPI_SUCCESS)
    return OOMPI_Intra_comm(MPI_COMM_NULL);

  return OOMPI_Intra_comm(mpi_comm, true);
}

