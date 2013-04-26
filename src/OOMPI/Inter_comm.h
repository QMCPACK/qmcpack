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

#ifndef _OOMPI_INTER_COMM_H_
#define _OOMPI_INTER_COMM_H_

#include "Intra_comm.h"


class OOMPI_Inter_comm : public OOMPI_Comm
{
  friend class OOMPI_Intra_comm;
  friend class OOMPI_Port;

public:
  // default constructor.  sets comm to MPI_COMM_NULL
  inline OOMPI_Inter_comm(MPI_Comm mpi_comm = MPI_COMM_NULL);

  // Does a shallow copy
  OOMPI_Inter_comm(const OOMPI_Inter_comm &a);
  OOMPI_Inter_comm &operator=(const OOMPI_Inter_comm &a);

  // Does an MPI_Intercomm_create
  OOMPI_Inter_comm(OOMPI_Intra_comm &local_comm, int local_leader,
                   OOMPI_Intra_comm &peer_comm, int remote_leader,
                   int tag = OOMPI_INTERCOMM_CREATE_TAG);

  virtual ~OOMPI_Inter_comm();

  // Does an MPI_Comm_dup and returns the new inter comm
  OOMPI_Inter_comm Dup(void);

  // ----------- Type - virtual function from comm
  // returns flag
  virtual bool Test_inter(void);

  // returns an OOMPI_Group that holds the mpi_group
  OOMPI_Group Remote_group(void);

  // returns size
  int Remote_size(void);

  // returns newintracomm
  OOMPI_Intra_comm Merge(bool high = true);

protected:
private:

  // Internal constructor (like default, but adds needs_to_be_freed)
  inline OOMPI_Inter_comm(MPI_Comm mpi_comm, bool needs_to_be_freed);
  void do_full_init(MPI_Comm mpi_comm, bool needs_to_be_freed);
  // This exists since constructors can't be called directly.
};

// MPI and default constructor
inline OOMPI_Inter_comm::OOMPI_Inter_comm(MPI_Comm mpi_comm)
{
  do_full_init(mpi_comm, false);
}

// Full constructor
inline OOMPI_Inter_comm::OOMPI_Inter_comm(MPI_Comm mpi_comm, bool needs_to_be_freed)
{
  do_full_init(mpi_comm, needs_to_be_freed);
}

#endif

