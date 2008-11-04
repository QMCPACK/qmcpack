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

#include <iostream>
#include <cstdlib>
#include <stdarg.h>

#include "Constants.h"
#include "Error.h"
#include "Comm.h"


//
// Global instance
//

OOMPI_Error_table OOMPI_ERROR;


//
// Class variable
//

bool OOMPI_Error_table::init = false;


//
// Callback MPI error function
//
void
OOMPI_Error_handler(MPI_Comm *mpi_comm, int *err, ...)
{
  OOMPI_Error_action action;
  char string[MPI_MAX_ERROR_STRING];
  int len;
  
  MPI_Error_string(*err, string, &len);

  OOMPI_errno = *err;

  // Find the data associated with the mpi_comm

  action = OOMPI_ERROR.Get_action(*mpi_comm);

  // Perform the action

  switch(action) {
  case OOMPI_ERRORS_ARE_FATAL:
    std::cerr << std::endl 
	      << "OOMPI detects MPI error " << *err << ": " 
	      << string << std::endl;
    std::cerr << "OOMPI is aborting (OOMPI_ERRORS_ARE_FATAL)" << std::endl;
    MPI_Abort(*mpi_comm, *err);
    // Should never return
    //OOMPI_Exit(1);
    exit(1);
    break;
  case OOMPI_ERRORS_RETURN:
    break;
  case OOMPI_ERRORS_EXCEPTION:
#if OOMPI_HAVE_EXCEPTIONS
    throw OOMPI_Error(mpi_comm, *err);
#else
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::cerr << std::endl << "OOMPI detects MPI error " << *err << ": " 
	      << string << std::endl
	      << "MPI error in MPI_COMM_WORLD rank " << rank << " of " 
	      << size << "." << std::endl
	      << "OOMPI set for exception handling (OOMPI_ERROR_EXCEPTIONS), " 
	      << std::endl
	      << "but OOMPI was compiled with no exception handling." 
	      << std::endl
	      << "Defaulting to OOMPI_ERRORS_RETURN." << std::endl;
    return;
#endif
  }
}


//
// Default constructor
//
OOMPI_Error_table::OOMPI_Error_table()
  : root(0)
{
  if (!init)
    valid = init = true;
  else
    valid = false;
}


//
// Destructor
//
OOMPI_Error_table::~OOMPI_Error_table()
{
  if (valid) {
    OOMPI_Error_entry *cur, *next = 0;
    
    for (cur = root; cur != 0; cur = next) {
      next = cur->next;
      delete cur;
    }
  }
}


//
// Convenience function
// Cannot be inline in .h file because of include file loop
//
void
OOMPI_Error_table::Handler(OOMPI_Comm *comm, int err)
{
  OOMPI_Error_handler(&(comm->Get_mpi()), &err);
}


//
// Table maintenance
//
void
OOMPI_Error_table::Add(MPI_Comm comm, OOMPI_Comm *oompi, 
		       OOMPI_Error_action action)
{
  if (comm == MPI_COMM_NULL || oompi == 0)
    return;

  OOMPI_Error_entry *cur, *prev, *n = new OOMPI_Error_entry;
  int result;

  n->mpi_comm = comm;
  n->oompi_comm = oompi;
  n->action = action;
  n->next = 0;

  if (root == (OOMPI_Error_entry *) 0) {
    root = n;
  }
  else {
    for (prev = root, cur = root; cur != 0; prev = cur, cur = cur->next) {
      // Check to see if some other instance of OOMPI_Comm is in the table
      // that has the same MPI_Comm
 
      MPI_Comm_compare(comm, cur->mpi_comm, &result);
      if (result == MPI_IDENT) {
	delete n;
	return;
      }
    }

    prev->next = n;
  }
}


void
OOMPI_Error_table::Change(MPI_Comm comm, OOMPI_Comm *oompi)
{
  OOMPI_Error_entry *cur = get_entry(comm);
  if (cur == 0) {
    Add(comm, oompi, OOMPI_ERRORS_ARE_FATAL);
    return;
  }

  cur->oompi_comm = oompi;
}


void
OOMPI_Error_table::Change(MPI_Comm comm, OOMPI_Error_action action)
{
  OOMPI_Error_entry *cur = get_entry(comm);
  if (cur == 0)
    return;

  cur->action = action;
}


void
OOMPI_Error_table::Delete(MPI_Comm comm)
{
  if (comm == MPI_COMM_NULL)
    return;

  OOMPI_Error_entry *cur, *prev;
  int result;

  for (prev = 0, cur = root; cur != 0; prev = cur, cur = cur->next) {
    MPI_Comm_compare(comm, cur->mpi_comm, &result);
    if (result == MPI_IDENT) {
      if (prev == 0) {
	root = cur->next;
	delete cur;
      }
      else {
	prev->next = cur->next;
	delete cur;
      }

      return;
    }
  }
}


//
// Given MPI_Comm, get OOMPI_Comm *
//
OOMPI_Comm *
OOMPI_Error_table::Get_oompi(MPI_Comm comm)
{
  OOMPI_Error_entry *cur = get_entry(comm);
  if (cur != 0)
    return cur->oompi_comm;

  return 0;
}


//
// Given MPI_Comm, get action
//
OOMPI_Error_action
OOMPI_Error_table::Get_action(MPI_Comm comm)
{
  OOMPI_Error_entry *cur = get_entry(comm);
  if (cur != 0)
    return cur->action;

  return OOMPI_ERRORS_ARE_FATAL;
}


//
// Do a lookup on the table
//
OOMPI_Error_entry *
OOMPI_Error_table::get_entry(MPI_Comm comm)
{
  if (comm == MPI_COMM_NULL)
    return 0;

  OOMPI_Error_entry *cur;
  int result(0);

  for (cur = root; cur != 0; cur = cur->next) {
    MPI_Comm_compare(comm, cur->mpi_comm, &result);
    if (result == MPI_IDENT)
      return cur;
  }

  return 0;
}

