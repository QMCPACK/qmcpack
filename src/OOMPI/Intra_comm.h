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
// Intra_communicators
//

#ifndef _OOMPI_INTRA_COMM_H_
#define _OOMPI_INTRA_COMM_H_

#include "Group.h"
#include "Comm.h"
#include "Inter_comm.h"
#include "Message.h"


class OOMPI_Cart_comm;
class OOMPI_Graph_comm;


class OOMPI_Intra_comm : public OOMPI_Comm
{

  friend class OOMPI_Inter_comm;
  friend class OOMPI_Graph_comm;
  friend class OOMPI_Cart_comm;

public:
  // default constructor.
  inline OOMPI_Intra_comm(MPI_Comm mpi_comm = MPI_COMM_NULL);

  // Shallow copy
  OOMPI_Intra_comm(const OOMPI_Intra_comm &a);

  // Shallow copy
  OOMPI_Intra_comm &operator=(const OOMPI_Intra_comm &a);

  // destructor
  virtual ~OOMPI_Intra_comm(void);


  //
  // Create
  //

  OOMPI_Intra_comm Create(OOMPI_Group &group);

  //
  // Dup
  //

  OOMPI_Intra_comm Dup(void);

  //
  // Split
  //

  OOMPI_Intra_comm Split(int color, int key = 0);

  //
  // Intercomm_create
  //

  OOMPI_Inter_comm Intercomm_create(int local_leader,
                                    OOMPI_Intra_comm &peer_comm,
                                    int remote_leader,
                                    int tag = OOMPI_INTERCOMM_CREATE_TAG);

  OOMPI_Inter_comm Intercomm_create(int local_leader,
                                    OOMPI_Port &peer_port,
                                    int tag = OOMPI_INTERCOMM_CREATE_TAG);

  //
  // Test_inter
  //

  virtual bool Test_inter(void);

  //
  // Collective communication
  // Allgather
  //

  void Allgather(OOMPI_Message sendbuf, OOMPI_Message recvbuf);
  void Allgather(OOMPI_Message sendbuf,
                 OOMPI_Array_message recvbuf, int recvcount);
  void Allgather(OOMPI_Array_message sendbuf, int sendcount,
                 OOMPI_Message recvbuf);
  void Allgather(OOMPI_Array_message sendbuf, int sendcount,
                 OOMPI_Array_message recvbuf, int recvcount);

  //
  // Allgatherv
  //

  void Allgatherv(OOMPI_Message sendbuf,
                  OOMPI_Array_message recvbuf, int recvcounts[],
                  int displs[]);
  void Allgatherv(OOMPI_Array_message sendbuf, int sendcount,
                  OOMPI_Array_message recvbuf, int recvcounts[],
                  int displs[]);

  //
  // Allreduce
  //

  void Allreduce(OOMPI_Message sendbuf,
                 OOMPI_Message recvbuf, const OOMPI_Op& op);
  void Allreduce(OOMPI_Array_message sendbuf, int sendcount,
                 OOMPI_Array_message recvbuf, const OOMPI_Op& op);

  //
  // Alltoall
  //

  void Alltoall(OOMPI_Message sendbuf, OOMPI_Message recvbuf);
  void Alltoall(OOMPI_Message sendbuf,
                OOMPI_Array_message recvbuf, int recvcount);
  void Alltoall(OOMPI_Array_message sendbuf, int sendcount,
                OOMPI_Message recvbuf);
  void Alltoall(OOMPI_Array_message sendbuf, int sendcount,
                OOMPI_Array_message recvbuf, int recvcount);

  //
  // Alltoallv
  //

  void Alltoallv(OOMPI_Array_message sendbuf, int sendcounts[],
                 int sdispls[], OOMPI_Array_message recvbuf,
                 int recvcounts[], int rdispls[]);

  //
  // Barrier
  //

  void Barrier(void);

  //
  // Reduce_scatter
  //

  void Reduce_scatter(OOMPI_Message sendbuf,
                      OOMPI_Array_message recvbuf, int recvcounts[],
                      const OOMPI_Op& op);
  void Reduce_scatter(OOMPI_Array_message sendbuf,
                      OOMPI_Array_message recvbuf, int recvcounts[],
                      const OOMPI_Op& op);

  //
  // Scan
  //

  void Scan(OOMPI_Message sendbuf,
            OOMPI_Array_message recvbuf, const OOMPI_Op& op);
  void Scan(OOMPI_Array_message sendbuf,
            OOMPI_Array_message recvbuf, int count, const OOMPI_Op& op);

  //
  // Receives on any port
  // operator>>
  //

  OOMPI_Intra_comm &operator>>(OOMPI_Message buf);

  //
  // Recv
  //

  OOMPI_Status Recv(OOMPI_Message buf, int tag = OOMPI_NO_TAG);
  OOMPI_Status Recv(OOMPI_Array_message buf, int count,
                    int tag = OOMPI_NO_TAG);

  //
  // Irecv
  //

  OOMPI_Request Irecv(OOMPI_Message buf, int tag = OOMPI_NO_TAG);
  OOMPI_Request Irecv(OOMPI_Array_message buf, int count,
                      int tag = OOMPI_NO_TAG);

  //
  // Recv_init
  //

  OOMPI_Request Recv_init(OOMPI_Message buf, int tag = OOMPI_NO_TAG);
  OOMPI_Request Recv_init(OOMPI_Array_message buf, int count,
                          int tag = OOMPI_NO_TAG);

  //
  // Probe
  //

  OOMPI_Status Probe(int tag);
  OOMPI_Status Iprobe(int tag, bool &flag);
  bool Iprobe(int tag);


  //
  // Topology
  // MPI_Topo_test
  //

  OOMPI_Topology Topo_test(void);

  // MPI_Graph_map

  int Graph_map(int nnodes, int index[], int edges[]);
  int Cart_map(int ndims, int dims[], bool periods[]);
  int Cart_map(int ndims, int dims[], bool periods);


  // MPI_Dims_create .  If the default is used then nnodes will be
  // set to the number of ports in the intra communicator
  // returns dims
  int *Dims_create(int ndims, int dims[], int nnodes = 0);

protected:

private:

  // normal constructor.  sets comm to MPI_COMM_NULL
  inline OOMPI_Intra_comm(MPI_Comm mpi_comm, bool needs_to_be_freed);
  void do_full_init(MPI_Comm mpi_comm, bool needs_to_be_freed);
  // This exists since constructors can't be called directly.

};


//
// OOMPI_COMM_SELF
//

extern OOMPI_Intra_comm OOMPI_COMM_SELF;

// MPI and default constructor
inline OOMPI_Intra_comm::OOMPI_Intra_comm(MPI_Comm mpi_comm)
  : OOMPI_Comm(MPI_COMM_NULL)
{
  do_full_init(mpi_comm, false);
}

// Standard constructor
inline OOMPI_Intra_comm::OOMPI_Intra_comm(MPI_Comm mpi_comm, bool needs_to_be_freed)
  : OOMPI_Comm(MPI_COMM_NULL)
{
  do_full_init(mpi_comm, needs_to_be_freed);
}

#endif
