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

#ifndef _OOMPI_PORT_H_
#define _OOMPI_PORT_H_

//
// Forward reference
//

class OOMPI_Port;
class OOMPI_Inter_comm;
class OOMPI_Intra_comm;
class OOMPI_Status;
class OOMPI_Comm;


#include <mpi.h>
#include "Op.h"
#include "Status.h"
#include "Request.h"
#include "Wrapper_ptr.cct"
#include "Message.h"


class OOMPI_Port
{
  friend class OOMPI_Comm;
  friend class OOMPI_Intra_comm;

public:

  //
  // Default constructor makes an invalid port
  //

  OOMPI_Port(void);
  OOMPI_Port(MPI_Comm a, int rank);
  OOMPI_Port(const OOMPI_Port &a);
  OOMPI_Port &operator=(const OOMPI_Port &a);

  //
  // Utility functions
  //

  inline int Rank(void);

  //
  // Receives
  //

  inline OOMPI_Port &operator>>(OOMPI_Message buf);

  inline OOMPI_Status Recv(OOMPI_Message buf, int tag = OOMPI_NO_TAG);
  inline OOMPI_Status Recv(OOMPI_Array_message buf, int count,
                           int tag = OOMPI_NO_TAG);

  inline OOMPI_Request Irecv(OOMPI_Message buf, int tag = OOMPI_NO_TAG);
  inline OOMPI_Request Irecv(OOMPI_Array_message buf, int count,
                             int tag = OOMPI_NO_TAG);

  inline OOMPI_Request Recv_init(OOMPI_Message buf, int tag = OOMPI_NO_TAG);
  inline OOMPI_Request Recv_init(OOMPI_Array_message buf, int count,
                                 int tag = OOMPI_NO_TAG);

  //
  // Probe
  //

  inline OOMPI_Status Probe(int tag);
  inline OOMPI_Status Iprobe(int tag, int &flag);

  //
  // Making intercommunicators
  //

  OOMPI_Inter_comm Intercomm_create(OOMPI_Intra_comm &peer_comm,
                                    int remote_leader,
                                    int tag= OOMPI_INTERCOMM_CREATE_TAG);
  OOMPI_Inter_comm Intercomm_create(OOMPI_Port &peer_port,
                                    int tag= OOMPI_INTERCOMM_CREATE_TAG);

  //
  // Rooted collective functions
  // Bcast
  //

  inline void Bcast(OOMPI_Message buf);
  inline void Bcast(OOMPI_Array_message buf, int count);

  //
  // Gather
  // Root functions
  inline void Gather(OOMPI_Message sendbuf, OOMPI_Message recvbuf);
  inline void Gather(OOMPI_Message sendbuf,
                     OOMPI_Array_message recvbuf, int recvcount = 1);
  inline void Gather(OOMPI_Array_message sendbuf, int sendcount,
                     OOMPI_Message recvbuf);
  inline void Gather(OOMPI_Array_message sendbuf, int sendcount,
                     OOMPI_Array_message recvbuf, int recvcount);
  // Non-root functions
  inline void Gather(OOMPI_Message sendbuf);
  inline void Gather(OOMPI_Array_message sendbuf, int sendcount);

  //
  // Gatherv
  // Root functions
  inline void Gatherv(OOMPI_Message sendbuf,
                      OOMPI_Array_message recvbuf, int recvcounts[], int displs[]);
  inline void Gatherv(OOMPI_Array_message sendbuf, int sendcount,
                      OOMPI_Array_message recvbuf, int recvcounts[], int displs[]);
  // Non-root functions
  inline void Gatherv(OOMPI_Message sendbuf);
  inline void Gatherv(OOMPI_Array_message sendbuf, int sendcount);

  //
  // Reduce
  // Root functions
  inline void Reduce(OOMPI_Message sendbuf,
                     OOMPI_Message recvbuf, const OOMPI_Op& op);
  inline void Reduce(OOMPI_Message sendbuf,
                     OOMPI_Array_message recvbuf, const OOMPI_Op& op);
  inline void Reduce(OOMPI_Array_message sendbuf,
                     OOMPI_Message recvbuf, const OOMPI_Op& op);
  inline void Reduce(OOMPI_Array_message sendbuf,
                     OOMPI_Array_message recvbuf, int count, const OOMPI_Op& op);
  // Non-root functions
  inline void Reduce(OOMPI_Message sendbuf, const OOMPI_Op& op);
  inline void Reduce(OOMPI_Array_message sendbuf, int sendcount,
                     const OOMPI_Op& op);

  //
  // Scatter
  // Root functions
  inline void Scatter(OOMPI_Message sendbuf, OOMPI_Message recvbuf);
  inline void Scatter(OOMPI_Message sendbuf,
                      OOMPI_Array_message recvbuf, int recvcount);
  inline void Scatter(OOMPI_Array_message sendbuf, int sendcount,
                      OOMPI_Message recvbuf);
  inline void Scatter(OOMPI_Array_message sendbuf, int sendcount,
                      OOMPI_Array_message recvbuf, int recvcount);
  // Non-root functions
  inline void Scatter(OOMPI_Message recvbuf);
  inline void Scatter(OOMPI_Array_message sendbuf, int recvcount);

  //
  // Scatterv
  // Root functions
  inline void Scatterv(OOMPI_Array_message sendbuf, int sendcounts[],
                       int displs[], OOMPI_Message recvbuf);
  inline void Scatterv(OOMPI_Array_message sendbuf, int sendcounts[],
                       int displs[], OOMPI_Array_message recvbuf,
                       int recvcount);
  // Non-root functions
  inline void Scatterv(OOMPI_Message recvbuf);
  inline void Scatterv(OOMPI_Array_message recvbuf, int recvcount);

  //
  // Sends
  //

  inline OOMPI_Port &operator<<(OOMPI_Message buf);

  inline void Bsend(OOMPI_Message buf, int tag = OOMPI_NO_TAG);
  inline void Bsend(OOMPI_Array_message buf, int count,
                    int tag = OOMPI_NO_TAG);

  inline void Rsend(OOMPI_Message buf, int tag = OOMPI_NO_TAG);
  inline void Rsend(OOMPI_Array_message buf, int count,
                    int tag = OOMPI_NO_TAG);

  inline void Send(OOMPI_Message buf, int tag = OOMPI_NO_TAG);
  inline void Send(OOMPI_Array_message buf, int count, int tag = OOMPI_NO_TAG);

  inline void Ssend(OOMPI_Message buf, int tag = OOMPI_NO_TAG);
  inline void Ssend(OOMPI_Array_message buf, int count,
                    int tag = OOMPI_NO_TAG);

  inline OOMPI_Request Bsend_init(OOMPI_Message buf,
                                  int tag = OOMPI_NO_TAG);
  inline OOMPI_Request Bsend_init(OOMPI_Array_message buf, int count,
                                  int tag = OOMPI_NO_TAG);

  inline OOMPI_Request Rsend_init(OOMPI_Message buf,
                                  int tag = OOMPI_NO_TAG);
  inline OOMPI_Request Rsend_init(OOMPI_Array_message buf, int count,
                                  int tag = OOMPI_NO_TAG);

  inline OOMPI_Request Send_init(OOMPI_Message buf,
                                 int tag = OOMPI_NO_TAG);
  inline OOMPI_Request Send_init(OOMPI_Array_message buf, int count,
                                 int tag = OOMPI_NO_TAG);

  inline OOMPI_Request Ssend_init(OOMPI_Message buf,
                                  int tag = OOMPI_NO_TAG);
  inline OOMPI_Request Ssend_init(OOMPI_Array_message buf, int count,
                                  int tag = OOMPI_NO_TAG);


  inline OOMPI_Request Ibsend(OOMPI_Message buf,
                              int tag = OOMPI_NO_TAG);
  inline OOMPI_Request Ibsend(OOMPI_Array_message buf, int count,
                              int tag = OOMPI_NO_TAG);

  inline OOMPI_Request Irsend(OOMPI_Message buf,
                              int tag = OOMPI_NO_TAG);
  inline OOMPI_Request Irsend(OOMPI_Array_message buf, int count,
                              int tag = OOMPI_NO_TAG);

  inline OOMPI_Request Isend(OOMPI_Message buf,
                             int tag = OOMPI_NO_TAG);
  inline OOMPI_Request Isend(OOMPI_Array_message buf, int count,
                             int tag = OOMPI_NO_TAG);

  inline OOMPI_Request Issend(OOMPI_Message buf,
                              int tag = OOMPI_NO_TAG);
  inline OOMPI_Request Issend(OOMPI_Array_message buf, int count,
                              int tag = OOMPI_NO_TAG);

protected:

  //
  // Protected attributes
  //

  OOMPI_Wrapper_ptr<MPI_Comm> comm_wrapper;
  int my_rank;
  bool writable;

private:
};


//
// Global instances
//

extern OOMPI_Port OOMPI_PORT_NULL;




//
// Utility functions
//
int
OOMPI_Port::Rank()
{
  return my_rank;
}


//
// Receives
//
OOMPI_Port &
OOMPI_Port::operator>>(OOMPI_Message buf)
{
  MPI_Status mpi_status;
  MPI_Recv(buf.Get_top(), buf.Get_count(), buf.Get_type(), my_rank,
           buf.Get_tag(), comm_wrapper->Get(), &mpi_status);
  return *this;
}


OOMPI_Status
OOMPI_Port::Recv(OOMPI_Message buf, int tag)
{
  MPI_Status mpi_status;
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Recv(buf.Get_top(), buf.Get_count(), buf.Get_type(), my_rank,
           my_tag, comm_wrapper->Get(), &mpi_status);
  return OOMPI_Status(mpi_status);
}


OOMPI_Status
OOMPI_Port::Recv(OOMPI_Array_message buf, int count, int tag)
{
  MPI_Status mpi_status;
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Recv(buf.Get_top(), count, buf.Get_type(), my_rank,
           my_tag, comm_wrapper->Get(), &mpi_status);
  return OOMPI_Status(mpi_status);
}


OOMPI_Request
OOMPI_Port::Irecv(OOMPI_Message buf, int tag)
{
  MPI_Request mpi_request;
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Irecv(buf.Get_top(), buf.Get_count(), buf.Get_type(),
            my_rank, my_tag, comm_wrapper->Get(), &mpi_request);
  return OOMPI_Request(mpi_request);
}


OOMPI_Request
OOMPI_Port::Irecv(OOMPI_Array_message buf, int count, int tag)
{
  MPI_Request mpi_request;
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Irecv(buf.Get_top(), count, buf.Get_type(),
            my_rank, my_tag, comm_wrapper->Get(), &mpi_request);
  return OOMPI_Request(mpi_request);
}


OOMPI_Request
OOMPI_Port::Recv_init(OOMPI_Message buf, int tag)
{
  MPI_Request mpi_request;
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Recv_init(buf.Get_top(), buf.Get_count(), buf.Get_type(),
                my_rank, my_tag, comm_wrapper->Get(), &mpi_request);
  return OOMPI_Request(mpi_request);
}


OOMPI_Request
OOMPI_Port::Recv_init(OOMPI_Array_message buf, int count, int tag)
{
  MPI_Request mpi_request;
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Recv_init(buf.Get_top(), count, buf.Get_type(),
                my_rank, my_tag, comm_wrapper->Get(), &mpi_request);
  return OOMPI_Request(mpi_request);
}


//
// Rooted collective functions
// Bcast
//
void
OOMPI_Port::Bcast(OOMPI_Message buf)
{
  MPI_Bcast(buf.Get_top(), buf.Get_count(), buf.Get_type(),
            my_rank, comm_wrapper->Get());
}


void
OOMPI_Port::Bcast(OOMPI_Array_message buf, int count)
{
  MPI_Bcast(buf.Get_top(), count, buf.Get_type(),
            my_rank, comm_wrapper->Get());
}


//
// Gather
// Root functions
//
void
OOMPI_Port::Gather(OOMPI_Message sendbuf, OOMPI_Message recvbuf)
{
  MPI_Gather(sendbuf.Get_top(), sendbuf.Get_count(),
             sendbuf.Get_type(), recvbuf.Get_top(),
             recvbuf.Get_count(), recvbuf.Get_type(),
             my_rank, comm_wrapper->Get());
}


void
OOMPI_Port::Gather(OOMPI_Message sendbuf,
                   OOMPI_Array_message recvbuf, int recvcount)
{
  MPI_Gather(sendbuf.Get_top(), sendbuf.Get_count(),
             sendbuf.Get_type(), recvbuf.Get_top(),
             recvcount, recvbuf.Get_type(),
             my_rank, comm_wrapper->Get());
}


void
OOMPI_Port::Gather(OOMPI_Array_message sendbuf, int sendcount,
                   OOMPI_Message recvbuf)
{
  MPI_Gather(sendbuf.Get_top(), sendcount,
             sendbuf.Get_type(), recvbuf.Get_top(),
             recvbuf.Get_count(), recvbuf.Get_type(),
             my_rank, comm_wrapper->Get());
}


void
OOMPI_Port::Gather(OOMPI_Array_message sendbuf, int sendcount,
                   OOMPI_Array_message recvbuf, int recvcount)
{
  MPI_Gather(sendbuf.Get_top(), sendcount,
             sendbuf.Get_type(), recvbuf.Get_top(),
             recvcount, recvbuf.Get_type(),
             my_rank, comm_wrapper->Get());
}


//
// Non-root functions
//
void
OOMPI_Port::Gather(OOMPI_Message sendbuf)
{
  // JMS LAM 6.0 is not right -- should be able to use MPI_DATATYPE_NULL
  MPI_Gather(sendbuf.Get_top(), sendbuf.Get_count(),
             sendbuf.Get_type(), 0, 0, MPI_INT,
             my_rank, comm_wrapper->Get());
}


void
OOMPI_Port::Gather(OOMPI_Array_message sendbuf, int sendcount)
{
  // JMS LAM 6.0 is not right -- should be able to use MPI_DATATYPE_NULL
  MPI_Gather(sendbuf.Get_top(), sendcount,
             sendbuf.Get_type(), 0, 0, MPI_INT,
             my_rank, comm_wrapper->Get());
}

//
// Gatherv
// Root functions
//
void
OOMPI_Port::Gatherv(OOMPI_Message sendbuf,
                    OOMPI_Array_message recvbuf, int recvcounts[],
                    int displs[])
{
  MPI_Gatherv(sendbuf.Get_top(), sendbuf.Get_count(),
              sendbuf.Get_type(), recvbuf.Get_top(), recvcounts, displs,
              recvbuf.Get_type(), my_rank, comm_wrapper->Get());
}


void
OOMPI_Port::Gatherv(OOMPI_Array_message sendbuf, int sendcount,
                    OOMPI_Array_message recvbuf, int recvcounts[],
                    int displs[])
{
  MPI_Gatherv(sendbuf.Get_top(), sendcount, sendbuf.Get_type(),
              recvbuf.Get_top(), recvcounts, displs,
              recvbuf.Get_type(), my_rank, comm_wrapper->Get());
}


//
// Non-root functions
//
void
OOMPI_Port::Gatherv(OOMPI_Message sendbuf)
{
  // JMS LAM 6.0 is not right -- should be able to use MPI_DATATYPE_NULL
  MPI_Gatherv(sendbuf.Get_top(), sendbuf.Get_count(),
              sendbuf.Get_type(), 0, 0, 0, MPI_INT,
              my_rank, comm_wrapper->Get());
}


void
OOMPI_Port::Gatherv(OOMPI_Array_message sendbuf, int sendcount)
{
  // JMS LAM 6.0 is not right -- should be able to use MPI_DATATYPE_NULL
  MPI_Gatherv(sendbuf.Get_top(), sendcount, sendbuf.Get_type(),
              0, 0, 0, MPI_INT, my_rank, comm_wrapper->Get());
}


//
// Reduce
// Root functions
//
void
OOMPI_Port::Reduce(OOMPI_Message sendbuf,
                   OOMPI_Message recvbuf, const OOMPI_Op& op)
{
  MPI_Reduce(sendbuf.Get_top(), recvbuf.Get_top(), sendbuf.Get_count(),
             sendbuf.Get_type(), op.op_wrapper->Get(),
             my_rank, comm_wrapper->Get());
}


void
OOMPI_Port::Reduce(OOMPI_Message sendbuf,
                   OOMPI_Array_message recvbuf, const OOMPI_Op& op)
{
  MPI_Reduce(sendbuf.Get_top(), recvbuf.Get_top(), sendbuf.Get_count(),
             sendbuf.Get_type(), op.op_wrapper->Get(),
             my_rank, comm_wrapper->Get());
}


void
OOMPI_Port::Reduce(OOMPI_Array_message sendbuf,
                   OOMPI_Message recvbuf, const OOMPI_Op& op)
{
  MPI_Reduce(sendbuf.Get_top(), recvbuf.Get_top(), recvbuf.Get_count(),
             sendbuf.Get_type(), op.op_wrapper->Get(),
             my_rank, comm_wrapper->Get());
}


void
OOMPI_Port::Reduce(OOMPI_Array_message sendbuf,
                   OOMPI_Array_message recvbuf, int count,
                   const OOMPI_Op& op)
{
  MPI_Reduce(sendbuf.Get_top(), recvbuf.Get_top(), count,
             sendbuf.Get_type(), op.op_wrapper->Get(),
             my_rank, comm_wrapper->Get());
}


//
// Non-root functions
//
void
OOMPI_Port::Reduce(OOMPI_Message sendbuf, const OOMPI_Op& op)
{
  // We shouldn't need to put in a "valid" recvbuf arg here.  Yuk.
  MPI_Reduce(sendbuf.Get_top(), (void *) 1, sendbuf.Get_count(),
             sendbuf.Get_type(), op.op_wrapper->Get(),
             my_rank, comm_wrapper->Get());
}


void
OOMPI_Port::Reduce(OOMPI_Array_message sendbuf, int sendcount,
                   const OOMPI_Op& op)
{
  // We shouldn't need to put in a "valid" recvbuf arg here.  Yuk.
  MPI_Reduce(sendbuf.Get_top(), (void *) 1, sendcount,
             sendbuf.Get_type(), op.op_wrapper->Get(),
             my_rank, comm_wrapper->Get());
}


//
// Scatter
// Root functions
//
void
OOMPI_Port::Scatter(OOMPI_Message sendbuf, OOMPI_Message recvbuf)
{
  MPI_Scatter(sendbuf.Get_top(), sendbuf.Get_count(),
              sendbuf.Get_type(), recvbuf.Get_top(), recvbuf.Get_count(),
              recvbuf.Get_type(), my_rank, comm_wrapper->Get());
}


void
OOMPI_Port::Scatter(OOMPI_Array_message sendbuf, int sendcount,
                    OOMPI_Message recvbuf)
{
  MPI_Scatter(sendbuf.Get_top(), sendcount,
              sendbuf.Get_type(), recvbuf.Get_top(), recvbuf.Get_count(),
              recvbuf.Get_type(), my_rank, comm_wrapper->Get());
}


void
OOMPI_Port::Scatter(OOMPI_Message sendbuf,
                    OOMPI_Array_message recvbuf, int recvcount)
{
  MPI_Scatter(sendbuf.Get_top(), sendbuf.Get_count(),
              sendbuf.Get_type(), recvbuf.Get_top(), recvcount,
              recvbuf.Get_type(), my_rank, comm_wrapper->Get());
}


void
OOMPI_Port::Scatter(OOMPI_Array_message sendbuf, int sendcount,
                    OOMPI_Array_message recvbuf, int recvcount)
{
  MPI_Scatter(sendbuf.Get_top(), sendcount,
              sendbuf.Get_type(), recvbuf.Get_top(), recvcount,
              recvbuf.Get_type(), my_rank, comm_wrapper->Get());
}


//
// Non-root functions
//
void
OOMPI_Port::Scatter(OOMPI_Message recvbuf)
{
  // JMS LAM 6.0 is not right -- should be able to use MPI_DATATYPE_NULL
  MPI_Scatter(0, 0, MPI_INT,
              recvbuf.Get_top(), recvbuf.Get_count(),
              recvbuf.Get_type(), my_rank, comm_wrapper->Get());
}


void
OOMPI_Port::Scatter(OOMPI_Array_message recvbuf, int recvcount)
{
  // JMS LAM 6.0 is not right -- should be able to use MPI_DATATYPE_NULL
  MPI_Scatter(0, 0, MPI_INT, recvbuf.Get_top(), recvcount,
              recvbuf.Get_type(), my_rank, comm_wrapper->Get());
}


//
// Scatterv
// Root functions
//
void
OOMPI_Port::Scatterv(OOMPI_Array_message sendbuf, int sendcounts[],
                     int displs[], OOMPI_Message recvbuf)
{
  MPI_Scatterv(sendbuf.Get_top(), sendcounts, displs,
               sendbuf.Get_type(), recvbuf.Get_top(), recvbuf.Get_count(),
               recvbuf.Get_type(), my_rank, comm_wrapper->Get());
}


void
OOMPI_Port::Scatterv(OOMPI_Array_message sendbuf, int sendcounts[],
                     int displs[], OOMPI_Array_message recvbuf, int recvcount)
{
  MPI_Scatterv(sendbuf.Get_top(), sendcounts, displs,
               sendbuf.Get_type(), recvbuf.Get_top(), recvcount,
               recvbuf.Get_type(), my_rank, comm_wrapper->Get());
}


//
// Non-root functions
//
void
OOMPI_Port::Scatterv(OOMPI_Message recvbuf)
{
  // JMS LAM 6.0 is not right -- should be able to use MPI_DATATYPE_NULL
  MPI_Scatterv(0, 0, 0, MPI_INT, recvbuf.Get_top(), recvbuf.Get_count(),
               recvbuf.Get_type(), my_rank, comm_wrapper->Get());
}


void
OOMPI_Port::Scatterv(OOMPI_Array_message recvbuf, int recvcount)
{
  // JMS LAM 6.0 is not right -- should be able to use MPI_DATATYPE_NULL
  MPI_Scatterv(0, 0, 0, MPI_INT, recvbuf.Get_top(), recvcount,
               recvbuf.Get_type(), my_rank, comm_wrapper->Get());
}


//
// Sends
//
OOMPI_Port &
OOMPI_Port::operator<<(OOMPI_Message buf)
{
  MPI_Send(buf.Get_top(), buf.Get_count(), buf.Get_type(),
           my_rank, buf.Get_tag(), comm_wrapper->Get());
  return *this;
}

void
OOMPI_Port::Bsend(OOMPI_Message buf, int tag)
{
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Bsend(buf.Get_top(), buf.Get_count(), buf.Get_type(),
            my_rank, my_tag, comm_wrapper->Get());
}


void
OOMPI_Port::Bsend(OOMPI_Array_message buf, int count, int tag)
{
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Bsend(buf.Get_top(), count, buf.Get_type(),
            my_rank, my_tag, comm_wrapper->Get());
}


void
OOMPI_Port::Rsend(OOMPI_Message buf, int tag)
{
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Rsend(buf.Get_top(), buf.Get_count(), buf.Get_type(),
            my_rank, my_tag, comm_wrapper->Get());
}


void
OOMPI_Port::Rsend(OOMPI_Array_message buf, int count, int tag)
{
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Rsend(buf.Get_top(), count, buf.Get_type(),
            my_rank, my_tag, comm_wrapper->Get());
}


void
OOMPI_Port::Send(OOMPI_Message buf, int tag)
{
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Send(buf.Get_top(), buf.Get_count(), buf.Get_type(),
           my_rank, my_tag, comm_wrapper->Get());
}


void
OOMPI_Port::Send(OOMPI_Array_message buf, int count, int tag)
{
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Send(buf.Get_top(), count, buf.Get_type(),
           my_rank, my_tag, comm_wrapper->Get());
}


void
OOMPI_Port::Ssend(OOMPI_Message buf, int tag)
{
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Ssend(buf.Get_top(), buf.Get_count(), buf.Get_type(),
            my_rank, my_tag, comm_wrapper->Get());
}


void
OOMPI_Port::Ssend(OOMPI_Array_message buf, int count, int tag)
{
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Ssend(buf.Get_top(), count, buf.Get_type(),
            my_rank, my_tag, comm_wrapper->Get());
}


OOMPI_Request
OOMPI_Port::Bsend_init(OOMPI_Message buf, int tag)
{
  MPI_Request mpi_request;
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Bsend_init(buf.Get_top(), buf.Get_count(), buf.Get_type(),
                 my_rank, my_tag, comm_wrapper->Get(), &mpi_request);
  return OOMPI_Request(mpi_request);
}


OOMPI_Request
OOMPI_Port::Bsend_init(OOMPI_Array_message buf, int count, int tag)
{
  MPI_Request mpi_request;
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Bsend_init(buf.Get_top(), count, buf.Get_type(),
                 my_rank, my_tag, comm_wrapper->Get(), &mpi_request);
  return OOMPI_Request(mpi_request);
}


OOMPI_Request
OOMPI_Port::Rsend_init(OOMPI_Message buf, int tag)
{
  MPI_Request mpi_request;
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Rsend_init(buf.Get_top(), buf.Get_count(), buf.Get_type(),
                 my_rank, my_tag, comm_wrapper->Get(), &mpi_request);
  return OOMPI_Request(mpi_request);
}


OOMPI_Request
OOMPI_Port::Rsend_init(OOMPI_Array_message buf, int count, int tag)
{
  MPI_Request mpi_request;
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Rsend_init(buf.Get_top(), count, buf.Get_type(),
                 my_rank, my_tag, comm_wrapper->Get(), &mpi_request);
  return OOMPI_Request(mpi_request);
}


OOMPI_Request
OOMPI_Port::Send_init(OOMPI_Message buf, int tag)
{
  MPI_Request mpi_request;
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Send_init(buf.Get_top(), buf.Get_count(), buf.Get_type(),
                my_rank, my_tag, comm_wrapper->Get(), &mpi_request);
  return OOMPI_Request(mpi_request);
}


OOMPI_Request
OOMPI_Port::Send_init(OOMPI_Array_message buf, int count, int tag)
{
  MPI_Request mpi_request;
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Send_init(buf.Get_top(), count, buf.Get_type(),
                my_rank, my_tag, comm_wrapper->Get(), &mpi_request);
  return OOMPI_Request(mpi_request);
}


OOMPI_Request
OOMPI_Port::Ssend_init(OOMPI_Message buf, int tag)
{
  MPI_Request mpi_request;
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Ssend_init(buf.Get_top(), buf.Get_count(), buf.Get_type(),
                 my_rank, my_tag, comm_wrapper->Get(), &mpi_request);
  return OOMPI_Request(mpi_request);
}


OOMPI_Request
OOMPI_Port::Ssend_init(OOMPI_Array_message buf, int count, int tag)
{
  MPI_Request mpi_request;
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Ssend_init(buf.Get_top(), count, buf.Get_type(),
                 my_rank, my_tag, comm_wrapper->Get(), &mpi_request);
  return OOMPI_Request(mpi_request);
}


OOMPI_Request
OOMPI_Port::Ibsend(OOMPI_Message buf, int tag)
{
  MPI_Request mpi_request;
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Ibsend(buf.Get_top(), buf.Get_count(), buf.Get_type(),
             my_rank, my_tag, comm_wrapper->Get(), &mpi_request);
  return OOMPI_Request(mpi_request);
}


OOMPI_Request
OOMPI_Port::Ibsend(OOMPI_Array_message buf, int count, int tag)
{
  MPI_Request mpi_request;
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Ibsend(buf.Get_top(), count, buf.Get_type(),
             my_rank, my_tag, comm_wrapper->Get(), &mpi_request);
  return OOMPI_Request(mpi_request);
}


OOMPI_Request
OOMPI_Port::Irsend(OOMPI_Message buf, int tag)
{
  MPI_Request mpi_request;
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Irsend(buf.Get_top(), buf.Get_count(), buf.Get_type(),
             my_rank, my_tag, comm_wrapper->Get(), &mpi_request);
  return OOMPI_Request(mpi_request);
}


OOMPI_Request
OOMPI_Port::Irsend(OOMPI_Array_message buf, int count, int tag)
{
  MPI_Request mpi_request;
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Irsend(buf.Get_top(), count, buf.Get_type(),
             my_rank, my_tag, comm_wrapper->Get(), &mpi_request);
  return OOMPI_Request(mpi_request);
}


OOMPI_Request
OOMPI_Port::Isend(OOMPI_Message buf, int tag)
{
  MPI_Request mpi_request;
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Isend(buf.Get_top(), buf.Get_count(), buf.Get_type(),
            my_rank, my_tag, comm_wrapper->Get(), &mpi_request);
  return OOMPI_Request(mpi_request);
}


OOMPI_Request
OOMPI_Port::Isend(OOMPI_Array_message buf, int count, int tag)
{
  MPI_Request mpi_request;
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Isend(buf.Get_top(), count, buf.Get_type(),
            my_rank, my_tag, comm_wrapper->Get(), &mpi_request);
  return OOMPI_Request(mpi_request);
}


OOMPI_Request
OOMPI_Port::Issend(OOMPI_Message buf, int tag)
{
  MPI_Request mpi_request;
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Issend(buf.Get_top(), buf.Get_count(), buf.Get_type(),
             my_rank, my_tag, comm_wrapper->Get(), &mpi_request);
  return OOMPI_Request(mpi_request);
}


OOMPI_Request
OOMPI_Port::Issend(OOMPI_Array_message buf, int count, int tag)
{
  MPI_Request mpi_request;
  int my_tag = (tag == OOMPI_NO_TAG) ? buf.Get_tag() : tag;
  MPI_Issend(buf.Get_top(), count, buf.Get_type(),
             my_rank, my_tag, comm_wrapper->Get(), &mpi_request);
  return OOMPI_Request(mpi_request);
}


//
// Probe
//

//
// Non-blocking probe on a port
//
OOMPI_Status
OOMPI_Port::Iprobe(int tag, int &flag)
{
  MPI_Status mpi_status;
  MPI_Iprobe(my_rank, tag, comm_wrapper->Get(), &flag, &mpi_status);
  return OOMPI_Status(mpi_status);
}


//
// Blocking probe on a port
//
OOMPI_Status
OOMPI_Port::Probe(int tag)
{
  MPI_Status mpi_status;
  MPI_Probe(my_rank, tag, comm_wrapper->Get(), &mpi_status);
  return OOMPI_Status(mpi_status);
}

#endif
