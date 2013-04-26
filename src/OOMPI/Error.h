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
// OOMPI Error handling
// Handling MPI errors in OOMPI
//

#ifndef _OOMPI_ERROR_H_
#define _OOMPI_ERROR_H_

//
// Other definitions
//

#include <mpi.h>
#include "Constants.h"
#include "Error_table.h"


//
// Use a forward reference instead of including "Comm.h" because:
//   a) We are only using OOMPI_Comm pointers here, and
//   b) it would create a circular #include problem
//
class OOMPI_Comm;


//
// OOMPI error class
// Used for throwing
//
class OOMPI_Error
{
public:
  inline OOMPI_Error()
    : oompi_comm(0), code(0)
  {};
  inline OOMPI_Error(MPI_Comm *a, int error)
  // JMS Calling OOMPI_Error_table without semaphore may not
  // be thread safe
    : oompi_comm(OOMPI_ERROR.Get_oompi(*a)), code(error)
  {};
  inline OOMPI_Error(OOMPI_Comm *a, int error)
    : oompi_comm(a), code(error)
  {};

  //
  // Access functions
  //

  inline OOMPI_Comm &Get_comm(void)
  {
    return *oompi_comm;
  };
  inline int Get_code(void)
  {
    return code;
  };
  inline int Get_class(void)
  {
    int c;
    MPI_Error_class(code, &c);
    return c;
  };
  inline char *Get_string(void)
  {
    char *msg = new char[MPI_MAX_ERROR_STRING];
    int len;
    MPI_Error_string(code, msg, &len);
    return msg;
  };
  inline void Get_string(char msg[], int &len)
  {
    MPI_Error_string(code, msg, &len);
  };

protected:
  OOMPI_Comm *oompi_comm;
  int code;

private:
};


#endif
