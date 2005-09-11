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
// Port class
//

#include "Error.h"
#include "Port.h"
#include "Packed.h"
#include "Intra_comm.h"
#include "Inter_comm.h"


//
// Global instances
//

OOMPI_Port OOMPI_PORT_NULL(MPI_COMM_WORLD, MPI_PROC_NULL);


//
// Default consutrctor
// Makes an invalid port
//
OOMPI_Port::OOMPI_Port(void)
: my_rank(OOMPI_PROC_NULL), writable(true)
{
  // This is bogus.  Should be able to set it to MPI_COMM_NULL, but the
  // IBM SP-2 implementation won't allow it.
  comm_wrapper->Set(MPI_COMM_WORLD, 0);
}


//
// Main constructor
//
OOMPI_Port::OOMPI_Port(MPI_Comm c, int rank) 
: my_rank(rank), writable(false)
{
  comm_wrapper->Set(c, 0);
}


//
// Copy constructor
//
OOMPI_Port::OOMPI_Port(const OOMPI_Port &a)
: comm_wrapper(a.comm_wrapper), my_rank(a.my_rank), writable(true)
{
}


//
// Assignment operator
//
OOMPI_Port &
OOMPI_Port::operator=(const OOMPI_Port &a)
{
  if (this != &a && writable) {
    comm_wrapper = a.comm_wrapper;
    my_rank = a.my_rank;
  }
  
  return *this;
}


//
// These functions aren't inlined because:
// - they take too long anyway
// - they create #include loops
//

//
// Making intercommunicators
//
OOMPI_Inter_comm
OOMPI_Port::Intercomm_create(OOMPI_Intra_comm &peer_comm, 
			     int remote_leader, int tag)
{
  MPI_Comm comm = comm_wrapper->Get();
  MPI_Comm mpi_peer_comm = peer_comm.Get_mpi();
  MPI_Comm mpi_newintercomm;
  
  MPI_Intercomm_create(comm, my_rank, mpi_peer_comm,
		       remote_leader, tag, &mpi_newintercomm);

  return OOMPI_Inter_comm(mpi_newintercomm, true);
}

OOMPI_Inter_comm
OOMPI_Port::Intercomm_create(OOMPI_Port &peer_port, int tag)
{
  MPI_Comm comm = comm_wrapper->Get();
  MPI_Comm mpi_peer_comm = peer_port.comm_wrapper->Get();
  MPI_Comm mpi_newintercomm;
  
  MPI_Intercomm_create(comm, my_rank, mpi_peer_comm,
		       peer_port.Rank(), tag, &mpi_newintercomm);

  return OOMPI_Inter_comm(mpi_newintercomm, true);
}


