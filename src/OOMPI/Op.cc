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

#include <iostream>

#include "Op.h"
#include "Comm_world.h"
#include "Hidden.h"


//
// Wrapper delete function
//
static int 
Free_mpi_op(MPI_Op *op_wrapper)
{
  if (op_wrapper != 0 && *op_wrapper != MPI_OP_NULL ) {
    if (OOMPI_COMM_WORLD.Finalized()) {
      std::cerr << "Attempt to free operator after finalize (ignored, probably resulting in memory leak)." 
		<< std::endl;
    } else {
      MPI_Op_free(op_wrapper);
    }
  }

  return MPI_SUCCESS;
}


//
// Constructor
//
OOMPI_Op::OOMPI_Op(MPI_Op a)
  : op_wrapper(a, 0)
{
}


//
// Copy constructor
//
OOMPI_Op::OOMPI_Op(const OOMPI_Op &a)
  : op_wrapper(a.op_wrapper)
{
}


//
// Assignment operator
//
OOMPI_Op &
OOMPI_Op::operator=(const OOMPI_Op &a)
{
  if (this != &a)
    op_wrapper = a.op_wrapper;

  return *this;
}


//
// Main Constructor
//
OOMPI_Op::OOMPI_Op(MPI_User_function *f, bool c)
{
  MPI_Op *op = new MPI_Op;

  if (MPI_Op_create(f, (int) c, op) != MPI_SUCCESS) {
    op_wrapper->Set(MPI_OP_NULL, Free_mpi_op);
    return;
  }

  op_wrapper->Set(*op, Free_mpi_op);
  delete op;
}


//
// Hidden constructor for making const intrinsic operators
//
OOMPI_Op::OOMPI_Op(MPI_Op a, OOMPI_Hidden& hidden)
{
  op_wrapper->Set(a, 0);
  hidden.do_nothing(); // For picky compilers
}
