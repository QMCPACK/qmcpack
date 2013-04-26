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
// Reduction operator class
//

#ifndef _OOMPI_OP_H_
#define _OOMPI_OP_H_

#include <mpi.h>
#include "oompi-config.h"
#include "Wrapper_ptr.cct"


class OOMPI_Hidden;


class OOMPI_Op
{
  friend class OOMPI_Port;
  friend class OOMPI_Intra_comm;


public:

  OOMPI_Op(MPI_Op a = MPI_OP_NULL);
  OOMPI_Op(const OOMPI_Op &a);
  OOMPI_Op &operator=(const OOMPI_Op &a);

  //
  // Main constructor
  //
  OOMPI_Op(MPI_User_function *f, bool commutative = true);

  //
  // Hidden constructor for making const intrinsic operators
  //
  OOMPI_Op(MPI_Op a, OOMPI_Hidden &hidden);

  //
  // Access functions
  //
  inline bool Is_null(void)
  {
    return (bool) (op_wrapper->Get() == MPI_OP_NULL);
  };
  inline MPI_Op &Get_mpi(void)
  {
    return op_wrapper->Get();
  };

private:
  OOMPI_Wrapper_ptr<MPI_Op> op_wrapper;
};

#endif
