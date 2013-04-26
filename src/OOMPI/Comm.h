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

#ifndef _OOMPI_COMM_H_
#define _OOMPI_COMM_H_

//
// Forward references
//

class OOMPI_Comm;
class OOMPI_Packed;
class OOMPI_Group;
class OOMPI_Any_port;
class OOMPI_Port;


//
// Other definitions
//

#include <mpi.h>
#include "Wrapper_ptr.cct"
#include "Port.h"
#include "Group.h"
#include "Error.h"


//
// Class declaraction
//

class OOMPI_Comm
{
  friend class OOMPI_Group;
  friend class OOMPI_Packed;
  friend class OOMPI_Port;

public:
  inline OOMPI_Comm(MPI_Comm mpi_comm = MPI_COMM_NULL)
    : err_handler(MPI_ERRHANDLER_NULL)
  {
    MPI_Constructor(mpi_comm, false);
  };
  OOMPI_Comm(const OOMPI_Comm& a);
  OOMPI_Comm& operator=(const OOMPI_Comm& a);
  virtual ~OOMPI_Comm();

  // returns the MPI_Comm
  inline MPI_Comm& Get_mpi()
  {
    return comm_wrapper->Get();
  };

  //
  // Communicator Management
  //

  // returns result
  OOMPI_Compare Compare(OOMPI_Comm& a);
  friend bool operator==(OOMPI_Comm& a, OOMPI_Comm& b);
  friend bool operator!=(OOMPI_Comm& a, OOMPI_Comm& b);

  // returns group
  OOMPI_Group Group(void);

  // returns rank
  int Rank();

  // returns size
  int Size();

  //
  // Packed
  // returns size
  int Pack_size(OOMPI_Datatype type, int count);
  int Pack_size(OOMPI_Message data);
  int Pack_size(OOMPI_Array_message data, int count);

  //
  // Testing
  //

  inline bool Is_null(void)
  {
    return (bool) (comm_wrapper->Get() == MPI_COMM_NULL);
  };

  // abstract function
  virtual bool Test_inter(void) = 0;

  // returns flag
  bool Initialized(void);

  // returns void
  void Abort(int errorcode = 1);

  //
  // Port access
  //

  OOMPI_Port operator[](int i);

  //
  // SendRecv
  //

  OOMPI_Status Sendrecv(OOMPI_Message sendbuf, int dest, int sendtag,
                        OOMPI_Message recvbuf, int source, int recvtag);
  OOMPI_Status Sendrecv(OOMPI_Array_message sendbuf, int sendcount,
                        int dest, int sendtag,
                        OOMPI_Message recvbuf, int source, int recvtag);
  OOMPI_Status Sendrecv(OOMPI_Message sendbuf, int dest, int sendtag,
                        OOMPI_Array_message recvbuf, int recvcount,
                        int source, int recvtag);
  OOMPI_Status Sendrecv(OOMPI_Array_message sendbuf, int sendcount,
                        int dest, int sendtag,
                        OOMPI_Array_message recvbuf, int recvcount,
                        int source, int recvtag);

  OOMPI_Status Sendrecv_replace(OOMPI_Message buf, int dest, int sendtag,
                                int source, int recvtag);
  OOMPI_Status Sendrecv_replace(OOMPI_Array_message buf, int count,
                                int dest, int sendtag,
                                int source, int recvtag);

  //
  // Exceptions
  //

  inline void Set_error_action(OOMPI_Error_action action)
  {
    Set_error_action(action, false);
  };
  OOMPI_Error_action Get_error_action(void);

protected:
  OOMPI_Wrapper_ptr<MPI_Comm> comm_wrapper;

  OOMPI_Any_port *any_port;
  OOMPI_Port **ports;
  int num_ports;

  // Error stuff

  MPI_Errhandler err_handler;

  // Protected helper functions

  void MPI_Constructor(MPI_Comm a, bool free_upon_dereference);
  void Set_error_action(OOMPI_Error_action action, bool have_sem);

private:

  // Private helper functions

  void delete_ports();
  void create_ports();

  void install_handler();
};

#endif
