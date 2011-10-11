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
// Communicator base class
//


#include <iostream>

#include <mpi.h>

#include "Comm.h"
#include "Port.h"
#include "Group.h"
#include "Error.h"
#include "Message.h"
#include "Comm_world.h"
#include "Util.h"
#include "Constants.h"


//
// Wrapper delete function
//
static int 
Free_mpi_comm(MPI_Comm *ptr)
{
  if (ptr != 0 && *ptr != MPI_COMM_NULL && !OOMPI_COMM_WORLD.Finalized() &&
      *ptr != MPI_COMM_WORLD && *ptr != MPI_COMM_SELF) {
    OOMPI_ERROR.Delete(*ptr);
    if (OOMPI_COMM_WORLD.Finalized())
      std::cerr << "OOMPI error: Attempt to free communicator after finalize (ignored, probably resulting in memory leak)." << std::endl;
    else
      MPI_Comm_free(ptr);
  }

  return MPI_SUCCESS;
}


// Free only the OOMPI error handler for the communicator.  This is
// mostly copied from Free_mpi_comm, so they need to be updated
// together.
static int
Free_mpi_comm_errhandler_only(MPI_Comm *ptr)
{
  if (ptr != 0 && *ptr != MPI_COMM_NULL && !OOMPI_COMM_WORLD.Finalized() &&
      *ptr != MPI_COMM_WORLD && *ptr != MPI_COMM_SELF)
    OOMPI_ERROR.Delete(*ptr);

  return MPI_SUCCESS;
}
  

OOMPI_Comm::OOMPI_Comm(const OOMPI_Comm &a)
  : comm_wrapper(a.comm_wrapper), num_ports(a.num_ports), 
    err_handler(a.err_handler)
{
  create_ports();
}


//
// Assignment operator
//
OOMPI_Comm &
OOMPI_Comm::operator=(const OOMPI_Comm &a)
{
  if (&a == this)
    return *this;

  if (num_ports > 0)
    delete_ports();
  
  // Copy MPI_COMM_NULL

  if (a.num_ports <= 0) {
    comm_wrapper->Set(MPI_COMM_NULL, 0);
    num_ports = 0;
  }
  
  // Shallow copy

  else {
    comm_wrapper = a.comm_wrapper;
    num_ports = a.num_ports;
  }
 
  create_ports();
  return *this;
}


//
// Destructor
//
OOMPI_Comm::~OOMPI_Comm()
{ 
  delete_ports(); 
};


//			  
// Communicator Management
// Return true MPI_Comm_compare
//

int 
OOMPI_Comm::Compare(OOMPI_Comm &a)
{
  int result;
  if (MPI_Comm_compare(comm_wrapper->Get(), a.comm_wrapper->Get(), 
		       &result) != MPI_SUCCESS)
    return MPI_UNEQUAL;

  return result;
}


//
// operator==
// Return true only if MPI_Comm_compare returns MPI_IDENT
//

bool
operator==(OOMPI_Comm &a, OOMPI_Comm &b)
{
  int result(0);
  // No use checking for MPI_SUCCESS here -- only have bool to return
  MPI_Comm_compare(a.comm_wrapper->Get(), b.comm_wrapper->Get(), 
		   &result);

  return (bool) (result != MPI_UNEQUAL);
}


//
// operator!=
// Return true only if MPI_Comm_compare returns MPI_IDENT
//

bool
operator!=(OOMPI_Comm &a, OOMPI_Comm &b)
{
  return (bool) !(a == b);
}


OOMPI_Group
OOMPI_Comm::Group()
{
  MPI_Group mpi_group(0);
  if (MPI_Comm_group(comm_wrapper->Get(), &mpi_group) != MPI_SUCCESS)
    return OOMPI_Group(MPI_GROUP_NULL);

  return OOMPI_Group(mpi_group);
}


int 
OOMPI_Comm::Rank()
{
  int rank(0);
  if (MPI_Comm_rank(comm_wrapper->Get(), &rank) != MPI_SUCCESS)
    return OOMPI_UNDEFINED;

  return rank;
}


int 
OOMPI_Comm::Size() 
{
  int size(0);
  if (MPI_Comm_size(comm_wrapper->Get(), &size) != MPI_SUCCESS)
    return OOMPI_UNDEFINED;

  return size;
}


//
// Packed
//
int 
OOMPI_Comm::Pack_size(OOMPI_Datatype type, int incount)
{
  int size(0);
  if (MPI_Pack_size(incount, type.Get_mpi(), comm_wrapper->Get(), 
		    &size) != MPI_SUCCESS)
    return OOMPI_UNDEFINED;

  return size;
}


int 
OOMPI_Comm::Pack_size(OOMPI_Message data)
{
  int size(0);
  if (MPI_Pack_size(data.Get_count(), data.Get_type(), comm_wrapper->Get(), 
		    &size) != MPI_SUCCESS)
    return OOMPI_UNDEFINED;

  return size;
}


int 
OOMPI_Comm::Pack_size(OOMPI_Array_message data, int incount)
{
  int size(0);
  if (MPI_Pack_size(incount, data.Get_type(), comm_wrapper->Get(), 
		    &size) != MPI_SUCCESS)
    return OOMPI_UNDEFINED;

  return size;
}

//
// Testing
//

bool
OOMPI_Comm::Initialized()
{
  int flag(0);
  MPI_Initialized(&flag);
  return (bool) flag;
}


void 
OOMPI_Comm::Abort(int errorcode)
{
  MPI_Abort(comm_wrapper->Get(), errorcode);
}


//
// Port access
//
OOMPI_Port
OOMPI_Comm::operator[](int i)
{
  if ((i < 0) || (i >= num_ports))
    return OOMPI_PORT_NULL;

  return *ports[i];
}


//
// Sendrecv
//
OOMPI_Status 
OOMPI_Comm::Sendrecv(OOMPI_Message sendbuf, int dest, int sendtag,
		     OOMPI_Message recvbuf, int source, int recvtag)
{
  MPI_Status mpi_status;
 
  if (MPI_Sendrecv(sendbuf.Get_top(), sendbuf.Get_count(),
		   sendbuf.Get_type(), dest, sendtag,
		   recvbuf.Get_top(), recvbuf.Get_count(), 
		   recvbuf.Get_type(), source, recvtag,
		   Get_mpi(), &mpi_status) != MPI_SUCCESS)
    return OOMPI_Status();
  
  return OOMPI_Status(mpi_status);
}
 
 
OOMPI_Status 
OOMPI_Comm::Sendrecv(OOMPI_Array_message sendbuf, int sendcount, 
                     int dest, int sendtag,
                     OOMPI_Message recvbuf, int source, int recvtag)
{
  MPI_Status mpi_status;
  
  if (MPI_Sendrecv(sendbuf.Get_top(), sendcount, 
		   sendbuf.Get_type(), dest, sendtag,
		   recvbuf.Get_top(), recvbuf.Get_count(), 
		   recvbuf.Get_type(), source, recvtag,
		   Get_mpi(), &mpi_status) != MPI_SUCCESS)
    return OOMPI_Status();
 
  return OOMPI_Status(mpi_status);
}
 
 
OOMPI_Status 
OOMPI_Comm::Sendrecv(OOMPI_Message sendbuf, int dest, int sendtag, 
                     OOMPI_Array_message recvbuf, int recvcount, int source, 
                     int recvtag)
{
  MPI_Status mpi_status;

  if (MPI_Sendrecv(sendbuf.Get_top(), sendbuf.Get_count(), 
		   sendbuf.Get_type(), dest, sendtag,
		   recvbuf.Get_top(), recvcount, 
		   recvbuf.Get_type(), source, recvtag,
		   Get_mpi(), &mpi_status) != MPI_SUCCESS)
    return OOMPI_Status();

  return OOMPI_Status(mpi_status);
}


OOMPI_Status 
OOMPI_Comm::Sendrecv(OOMPI_Array_message sendbuf, int sendcount, 
		     int dest, int sendtag,
                     OOMPI_Array_message recvbuf, int recvcount, 
                     int source, int recvtag)
{
  MPI_Status mpi_status;
  
  if (MPI_Sendrecv(sendbuf.Get_top(), sendcount, 
		   sendbuf.Get_type(), dest, sendtag,
		   recvbuf.Get_top(), recvcount, 
		   recvbuf.Get_type(), source, recvtag,
		   Get_mpi(), &mpi_status) != MPI_SUCCESS)
    return OOMPI_Status();
 
  return OOMPI_Status(mpi_status);
}
 
 
//
// Sendrecv_replace
//
OOMPI_Status 
OOMPI_Comm::Sendrecv_replace(OOMPI_Message buf, int dest, int sendtag,
                             int source, int recvtag)
{
  MPI_Status mpi_status;

  if (MPI_Sendrecv_replace(buf.Get_top(), buf.Get_count(), 
			   buf.Get_type(), dest, sendtag,
			   source, recvtag, Get_mpi(), 
			   &mpi_status) != MPI_SUCCESS)
    return OOMPI_Status();

  return OOMPI_Status(mpi_status);
}
 
 
OOMPI_Status 
OOMPI_Comm::Sendrecv_replace(OOMPI_Array_message buf, int count, 
                             int dest, int sendtag,
                             int source, int recvtag)
{
  MPI_Status mpi_status;
 
  if (MPI_Sendrecv_replace(buf.Get_top(), count,
			   buf.Get_type(), dest, sendtag,
			   source, recvtag, Get_mpi(), 
			   &mpi_status) != MPI_SUCCESS)
    return OOMPI_Status();
 
  return OOMPI_Status(mpi_status);
}


//
// Exceptions
//
void
OOMPI_Comm::Set_error_action(OOMPI_Error_action action, bool have_sem)
{
  OOMPI_Util util;

  if (!have_sem)
    util.Get_sem(SEM_ERR_HANDLER);

  OOMPI_ERROR.Change(comm_wrapper->Get(), action);

  if (!have_sem)
    util.Release_sem(SEM_ERR_HANDLER);
}


//
// Find current value of error handling status
//
OOMPI_Error_action
OOMPI_Comm::Get_error_action()
{
  OOMPI_Util util;
  OOMPI_Error_action action;

  util.Get_sem(SEM_ERR_HANDLER);

  action = OOMPI_ERROR.Get_action(comm_wrapper->Get());

  util.Release_sem(SEM_ERR_HANDLER);
  return action;
}


//
// Protected helper function
//
void
OOMPI_Comm::MPI_Constructor(MPI_Comm mpi_comm, bool free_upon_dereference)
{
  if (free_upon_dereference)
    comm_wrapper->Set(mpi_comm, Free_mpi_comm);
  else
    comm_wrapper->Set(mpi_comm, Free_mpi_comm_errhandler_only);

  if (mpi_comm == MPI_COMM_SELF)
    num_ports = 1;
  else if (mpi_comm != MPI_COMM_NULL) {
    if (MPI_Comm_size(mpi_comm, &num_ports) != MPI_SUCCESS) {
      num_ports = 0;
      return;
    }
  }
  else
    num_ports = 0;

  create_ports();

  if (mpi_comm != MPI_COMM_NULL && mpi_comm != MPI_COMM_SELF)
    install_handler();
}


//
// Delete ports
//
void
OOMPI_Comm::delete_ports()
{
  int i;

  if (ports != 0) {
    for (i = 0; i < num_ports; i++) 
      delete ports[i];
    delete [] ports;
    ports=0;
  }
}


//
// Create ports
//
void
OOMPI_Comm::create_ports()
{
  int i;

  if (num_ports == 0) {
    ports = 0;
    return;
  }

  MPI_Comm comm = comm_wrapper->Get();
  ports = new OOMPI_Port *[num_ports];
  for (i = 0; i < num_ports; i++)
    ports[i] = new OOMPI_Port(comm, i);
}


//
// This routine only called when OOMPI_Comm's are made either by
// default or via MPI consturctor
//
void
OOMPI_Comm::install_handler()
{
  OOMPI_Util util;
  util.Get_sem(SEM_ERR_HANDLER);

  MPI_Comm comm = comm_wrapper->Get();
  if (comm != MPI_COMM_NULL) {
    OOMPI_ERROR.Add(comm, this, OOMPI_ERRORS_ARE_FATAL);
    MPI_Errhandler_create((MPI_Handler_function *) OOMPI_Error_handler, 
			  &err_handler);
    MPI_Errhandler_set(comm, err_handler);
    MPI_Errhandler_free(&err_handler);
    // The MPI 1.1 spec claims that it is OK to free the error handler right away, as the implementation will handle deleting it at the appropriate time.
  }
  util.Release_sem(SEM_ERR_HANDLER);
}
