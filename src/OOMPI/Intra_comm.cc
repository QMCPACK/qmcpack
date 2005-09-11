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

#include <mpi.h>
#include "Comm.h"
#include "Intra_comm.h"
#include "Inter_comm.h"
#include "oompi-config.h"
#include "Cart_comm.h"
#include "Graph_comm.h"
#include "Op.h"
#include "Constants.h"


//
// Instiantiate OOMPI_COMM_SELF
//

OOMPI_Intra_comm OOMPI_COMM_SELF(MPI_COMM_SELF);


void OOMPI_Intra_comm::do_full_init(MPI_Comm mpi_comm, bool needs_to_be_freed)
{
  if (mpi_comm == MPI_COMM_NULL)
    return;

  int flag(0);
  if (mpi_comm != MPI_COMM_SELF)
    MPI_Comm_test_inter(mpi_comm, &flag);

  if (flag == 0)
    MPI_Constructor(mpi_comm, needs_to_be_freed);
}


OOMPI_Intra_comm::OOMPI_Intra_comm(const OOMPI_Intra_comm &a)
: OOMPI_Comm(a)
{
}


//
// operator=
//
OOMPI_Intra_comm &
OOMPI_Intra_comm::operator=(const OOMPI_Intra_comm &a)
{
  (OOMPI_Comm &) *this = (OOMPI_Comm &) a; 
  return *this;
}


OOMPI_Intra_comm::~OOMPI_Intra_comm()
{
}


//
// Create
//
OOMPI_Intra_comm
OOMPI_Intra_comm::Create(OOMPI_Group &ingroup)
{
  MPI_Comm mpi_comm;
  if (MPI_Comm_create(comm_wrapper->Get(), ingroup.Get_mpi(), &mpi_comm) !=
      MPI_SUCCESS)
    return OOMPI_Intra_comm(MPI_COMM_NULL);

  return OOMPI_Intra_comm(mpi_comm, true);
}


//
// Dup
//
OOMPI_Intra_comm
OOMPI_Intra_comm::Dup()
{
  MPI_Comm mpi_comm;
  if (MPI_Comm_dup(Get_mpi(), &mpi_comm) != MPI_SUCCESS)
    return OOMPI_Intra_comm(MPI_COMM_NULL);

  return OOMPI_Intra_comm(mpi_comm, true);
}


//
// Split
//
OOMPI_Intra_comm
OOMPI_Intra_comm::Split(int color, int key)
{
  MPI_Comm mpi_comm;
  if (MPI_Comm_split(comm_wrapper->Get(), color, key, &mpi_comm) !=
      MPI_SUCCESS)
    return OOMPI_Intra_comm(MPI_COMM_NULL);

  return OOMPI_Intra_comm(mpi_comm, true);
}


//
// Intercomm_create
//
OOMPI_Inter_comm
OOMPI_Intra_comm::Intercomm_create(int local_leader, 
				   OOMPI_Intra_comm &peer_comm,
				   int remote_leader, 
				   int tag)
{
  MPI_Comm mpi_comm;
  if (MPI_Intercomm_create(comm_wrapper->Get(), local_leader, 
			   peer_comm.Get_mpi(), remote_leader, tag, 
			   &mpi_comm) != MPI_SUCCESS)
    return OOMPI_Inter_comm(MPI_COMM_NULL);

  return OOMPI_Inter_comm(mpi_comm, true);
}

			       
OOMPI_Inter_comm
OOMPI_Intra_comm::Intercomm_create(int local_leader, 
				   OOMPI_Port &peer_port,
				   int tag)
{
  MPI_Comm mpi_comm;
  
  if (MPI_Intercomm_create(comm_wrapper->Get(), local_leader, 
			   peer_port.comm_wrapper->Get(),
			   peer_port.Rank(), tag, &mpi_comm) != MPI_SUCCESS)
    return OOMPI_Inter_comm(MPI_COMM_NULL);

  return OOMPI_Inter_comm(mpi_comm, true);
} 


//
// Test_inter
//
bool
OOMPI_Intra_comm::Test_inter()
{
  int flag(0);
  
  MPI_Comm_test_inter(comm_wrapper->Get(), &flag);

  return (bool) flag;
}


//
// Collective communication
// Allgather
//
void 
OOMPI_Intra_comm::Allgather(OOMPI_Message sendbuf, OOMPI_Message recvbuf)
{
  MPI_Allgather(sendbuf.Get_top(), sendbuf.Get_count(), sendbuf.Get_type(),
		recvbuf.Get_top(), recvbuf.Get_count(), recvbuf.Get_type(),
		comm_wrapper->Get());
}


void 
OOMPI_Intra_comm::Allgather(OOMPI_Array_message sendbuf, int sendcount, 
			    OOMPI_Message recvbuf)
{
  MPI_Allgather(sendbuf.Get_top(), sendcount, sendbuf.Get_type(),
		recvbuf.Get_top(), recvbuf.Get_count(), recvbuf.Get_type(),
		comm_wrapper->Get());
}


void 
OOMPI_Intra_comm::Allgather(OOMPI_Message sendbuf, 
			    OOMPI_Array_message recvbuf, int recvcount)
{
  MPI_Allgather(sendbuf.Get_top(), sendbuf.Get_count(), sendbuf.Get_type(),
		recvbuf.Get_top(), recvcount, recvbuf.Get_type(),
		comm_wrapper->Get());
}


void 
OOMPI_Intra_comm::Allgather(OOMPI_Array_message sendbuf, int sendcount,
			    OOMPI_Array_message recvbuf, int recvcount)
{
  MPI_Allgather(sendbuf.Get_top(), sendcount, sendbuf.Get_type(),
		recvbuf.Get_top(), recvcount, recvbuf.Get_type(),
		comm_wrapper->Get());
}


//
// Allgatherv
//
void 
OOMPI_Intra_comm::Allgatherv(OOMPI_Message sendbuf,
			     OOMPI_Array_message recvbuf, int recvcounts[],
			     int displs[])
{
  MPI_Allgatherv(sendbuf.Get_top(), sendbuf.Get_count(), sendbuf.Get_type(),
		 recvbuf.Get_top(), recvcounts, displs,
		 recvbuf.Get_type(), comm_wrapper->Get());
}

void 
OOMPI_Intra_comm::Allgatherv(OOMPI_Array_message sendbuf, int sendcount,
			     OOMPI_Array_message recvbuf, int recvcounts[],
			     int displs[])
{
  MPI_Allgatherv(sendbuf.Get_top(), sendcount, sendbuf.Get_type(),
		 recvbuf.Get_top(), recvcounts, displs,
		 recvbuf.Get_type(), comm_wrapper->Get());
}


//
// Allreduce
//
void 
OOMPI_Intra_comm::Allreduce(OOMPI_Array_message sendbuf, int sendcount,
			    OOMPI_Array_message recvbuf, const OOMPI_Op& op)
{
  MPI_Allreduce(sendbuf.Get_top(), recvbuf.Get_top(), sendcount,
		sendbuf.Get_type(), op.op_wrapper->Get(), 
		comm_wrapper->Get());
}


void 
OOMPI_Intra_comm::Allreduce(OOMPI_Message sendbuf, 
			    OOMPI_Message recvbuf, const OOMPI_Op& op)
			    
{
  MPI_Allreduce(sendbuf.Get_top(), recvbuf.Get_top(), sendbuf.Get_count(),
		sendbuf.Get_type(), op.op_wrapper->Get(), 
		comm_wrapper->Get());
}


//
// Alltoall
//
void 
OOMPI_Intra_comm::Alltoall(OOMPI_Message sendbuf, OOMPI_Message recvbuf)
{
  MPI_Alltoall(sendbuf.Get_top(), sendbuf.Get_count(), sendbuf.Get_type(),
	       recvbuf.Get_top(), recvbuf.Get_count(), recvbuf.Get_type(),
	       comm_wrapper->Get());
}


void 
OOMPI_Intra_comm::Alltoall(OOMPI_Message sendbuf, 
			   OOMPI_Array_message recvbuf, int recvcount)
{
  MPI_Alltoall(sendbuf.Get_top(), sendbuf.Get_count(), sendbuf.Get_type(),
	       recvbuf.Get_top(), recvcount, recvbuf.Get_type(),
	       comm_wrapper->Get());
}


void 
OOMPI_Intra_comm::Alltoall(OOMPI_Array_message sendbuf, int sendcount,
			   OOMPI_Message recvbuf)
{
  MPI_Alltoall(sendbuf.Get_top(), sendcount, sendbuf.Get_type(),
	       recvbuf.Get_top(), recvbuf.Get_count(), recvbuf.Get_type(),
	       comm_wrapper->Get());
}


void 
OOMPI_Intra_comm::Alltoall(OOMPI_Array_message sendbuf, int sendcount,
			   OOMPI_Array_message recvbuf, int recvcount)
{
  MPI_Alltoall(sendbuf.Get_top(), sendcount, sendbuf.Get_type(),
	       recvbuf.Get_top(), recvcount, recvbuf.Get_type(),
	       comm_wrapper->Get());
}


//
// Alltoallv
//
void 
OOMPI_Intra_comm::Alltoallv(OOMPI_Array_message sendbuf, int sendcounts[], 
			    int sdispls[], OOMPI_Array_message recvbuf, 
			    int recvcounts[], int rdispls[])
{
  MPI_Alltoallv(sendbuf.Get_top(), sendcounts, sdispls,
		sendbuf.Get_type(), recvbuf.Get_top(),
		recvcounts, rdispls, recvbuf.Get_type(),
		comm_wrapper->Get());
}


//
// Barrier
//
void 
OOMPI_Intra_comm::Barrier()
{
  MPI_Barrier(comm_wrapper->Get());
}

 
//
// Reduce scatter
//
void 
OOMPI_Intra_comm::Reduce_scatter(OOMPI_Message sendbuf, 
				 OOMPI_Array_message recvbuf, int recvcounts[],
				 const OOMPI_Op& op)
{
  MPI_Reduce_scatter(sendbuf.Get_top(), recvbuf.Get_top(),
		     recvcounts, sendbuf.Get_type(),
		     op.op_wrapper->Get(), comm_wrapper->Get());
}


void 
OOMPI_Intra_comm::Reduce_scatter(OOMPI_Array_message sendbuf,
				 OOMPI_Array_message recvbuf, int recvcounts[],
				 const OOMPI_Op& op)
{
  MPI_Reduce_scatter(sendbuf.Get_top(), recvbuf.Get_top(),
		     recvcounts, sendbuf.Get_type(),
		     op.op_wrapper->Get(), comm_wrapper->Get());
}


//
// Scan
//
void 
OOMPI_Intra_comm::Scan(OOMPI_Message sendbuf, 
		       OOMPI_Array_message recvbuf, const OOMPI_Op& op)
{
  MPI_Scan(sendbuf.Get_top(), recvbuf.Get_top(), sendbuf.Get_count(), 
	   sendbuf.Get_type(), op.op_wrapper->Get(), comm_wrapper->Get());
}


void 
OOMPI_Intra_comm::Scan(OOMPI_Array_message sendbuf, 
		       OOMPI_Array_message recvbuf, int count, 
		       const OOMPI_Op& op)
{
  MPI_Scan(sendbuf.Get_top(), recvbuf.Get_top(), count, 
	   sendbuf.Get_type(), op.op_wrapper->Get(), comm_wrapper->Get());
}

//
// Receives on any port
// operator>>
//

OOMPI_Intra_comm &
OOMPI_Intra_comm::operator>>(OOMPI_Message buf)
{
  MPI_Status mpi_status;
  MPI_Recv(buf.Get_top(), buf.Get_count(), buf.Get_type(), MPI_ANY_SOURCE,
	   buf.Get_tag(), comm_wrapper->Get(), &mpi_status);

  return *this;
}


//
// Recv
//
OOMPI_Status 
OOMPI_Intra_comm::Recv(OOMPI_Message buf, int tag)
{
  MPI_Status mpi_status;
  int my_tag= (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Recv(buf.Get_top(), buf.Get_count(), buf.Get_type(), MPI_ANY_SOURCE, 
	   my_tag, comm_wrapper->Get(), &mpi_status);

  return OOMPI_Status(mpi_status);
}


OOMPI_Status
OOMPI_Intra_comm::Recv(OOMPI_Array_message buf, int count, int tag)
{
  MPI_Status mpi_status;
  int my_tag= (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Recv(buf.Get_top(), count, buf.Get_type(), MPI_ANY_SOURCE, 
	   my_tag, comm_wrapper->Get(), &mpi_status);

  return OOMPI_Status(mpi_status);
}


//
// Irecv
//
OOMPI_Request
OOMPI_Intra_comm::Irecv(OOMPI_Message buf, int tag)
{
  MPI_Request mpi_request;
  int my_tag= (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Irecv(buf.Get_top(), buf.Get_count(), buf.Get_type(),
	    MPI_ANY_SOURCE, my_tag, comm_wrapper->Get(), &mpi_request);

  return OOMPI_Request(mpi_request);
}


OOMPI_Request
OOMPI_Intra_comm::Irecv(OOMPI_Array_message buf, int count, int tag)
{
  MPI_Request mpi_request;
  int my_tag= (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Irecv(buf.Get_top(), count, buf.Get_type(),
	    MPI_ANY_SOURCE, my_tag, comm_wrapper->Get(), &mpi_request);

  return OOMPI_Request(mpi_request);
}


//
// Recv_init
//
OOMPI_Request
OOMPI_Intra_comm::Recv_init(OOMPI_Message buf, int tag)
{
  MPI_Request mpi_request;
  int my_tag= (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Recv_init(buf.Get_top(), buf.Get_count(), buf.Get_type(),
		MPI_ANY_SOURCE, my_tag, comm_wrapper->Get(), &mpi_request);

  return OOMPI_Request(mpi_request);
}


OOMPI_Request
OOMPI_Intra_comm::Recv_init(OOMPI_Array_message buf, int count, int tag)
{
  MPI_Request mpi_request;
  int my_tag= (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Recv_init(buf.Get_top(), count, buf.Get_type(),
		MPI_ANY_SOURCE, my_tag, comm_wrapper->Get(), &mpi_request);

  return OOMPI_Request(mpi_request);
}


//
// Probe on any port
//

OOMPI_Status
OOMPI_Intra_comm::Probe(int tag)
{
  MPI_Status mpi_status;
  MPI_Probe(MPI_ANY_SOURCE, tag, comm_wrapper->Get(), &mpi_status);

  return OOMPI_Status(mpi_status);
}


OOMPI_Status
OOMPI_Intra_comm::Iprobe(int tag, bool &flag)
{
  MPI_Status mpi_status;
#if OOMPI_BOOL_NE_INT
  int i;
  MPI_Iprobe(MPI_ANY_SOURCE, tag, comm_wrapper->Get(), &i, &mpi_status);
  flag = (bool) i;
#else
  MPI_Iprobe(MPI_ANY_SOURCE, tag, comm_wrapper->Get(), (int *) &flag, 
	     &mpi_status);
#endif

  return OOMPI_Status(mpi_status);
}


bool
OOMPI_Intra_comm::Iprobe(int tag)
{
  MPI_Status mpi_status;
  int flag;

  MPI_Iprobe(MPI_ANY_SOURCE, tag, comm_wrapper->Get(), &flag, &mpi_status);

  return (bool) flag;
}


//
// Topology
//

// returns status
OOMPI_Topology
OOMPI_Intra_comm::Topo_test()
{
  int status;
  MPI_Topo_test(comm_wrapper->Get(), &status);
  return (OOMPI_Topology) status;
}


// MPI_Dims_create .  If the default is used then nnodes will be
// set to the number of ports in the intra communicator
int *
OOMPI_Intra_comm::Dims_create(int ndims, int dims[], int nnodes)
{
  if (nnodes == 0)
    nnodes = num_ports;
  
  MPI_Dims_create(nnodes, ndims, dims);
  return dims;
}


// MPI_Cart_map
int 
OOMPI_Intra_comm::Cart_map(int ndims, int dims[], bool periods[])
{
  int newrank(0);
#if OOMPI_BOOL_NE_INT
  int *int_ptr = new int[ndims];
  for (int i = 0; i < ndims; i++)
    int_ptr[i] = (int) periods[i];
  if (MPI_Cart_map(comm_wrapper->Get(), ndims, dims, int_ptr, 
		   &newrank) != MPI_SUCCESS) {
    delete[] int_ptr;
    return OOMPI_UNDEFINED;
  }

  delete[] int_ptr;
#else
  if (MPI_Cart_map(comm_wrapper->Get(), ndims, dims, (int *) periods, 
		   &newrank) != MPI_SUCCESS)
    return OOMPI_UNDEFINED;
#endif
  return newrank;
}

// MPI_Cart_map
int 
OOMPI_Intra_comm::Cart_map(int ndims, int dims[], bool periods)
{
  int newrank(0);
  int *period_array = new int [ndims];
  for (int i = 0; i < ndims; i++)
    period_array[i] = (int) periods;

  if (MPI_Cart_map(comm_wrapper->Get(), ndims, dims, (int *) period_array, 
	       &newrank) != MPI_SUCCESS)
    return OOMPI_UNDEFINED;
  delete[] period_array;
  return newrank;
}

// MPI_Graph_map
int 
OOMPI_Intra_comm::Graph_map(int nnodes, int index[], int edges[])
{
  int newrank;
  if (MPI_Graph_map(comm_wrapper->Get(), nnodes, index, edges, &newrank) !=
      MPI_SUCCESS)
    return OOMPI_UNDEFINED;

  return newrank;
}
