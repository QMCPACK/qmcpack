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

#ifndef _OOMPI_ERROR_TABLE_H_
#define _OOMPI_ERROR_TABLE_H_


//
// Other definitions
//

#include <mpi.h>
#include "oompi-config.h"
#include "Constants.h"


//
// Use a forward reference instead of including "Comm.h" because:
//   a) We are only using OOMPI_Comm pointers here, and
//   b) it would create a circular #include problem
//
class OOMPI_Comm;

//
// Forward declaration for function used as friend in OOMPI_Error_table
// class
//
extern void OOMPI_Error_handler(MPI_Comm *mpi_comm, int *err, ...);

//
// Structure used for entries in the error lookup table
//

struct OOMPI_Error_entry
{
  MPI_Comm mpi_comm;
  OOMPI_Comm *oompi_comm;
  OOMPI_Error_action action;
  OOMPI_Error_entry *next;
};


//
// OOMPI error table class
// Maintains communications and actions to be performed when something
// Bad happens
//
class OOMPI_Error_table
{
  friend void OOMPI_Error_handler(MPI_Comm *mpi_comm, int *err, ...);

public:
  OOMPI_Error_table();
  inline OOMPI_Error_table(const OOMPI_Error_table &t)
    : valid(false), root(0)
  {
    t.do_nothing();
  };
  inline OOMPI_Error_table &operator=(const OOMPI_Error_table &t)
  {
    valid = false;
    root = 0;
    t.do_nothing();
    return *this;
  };
  ~OOMPI_Error_table();

  // Convenience routines

  inline void Handler(MPI_Comm comm, int err)
  {
    OOMPI_Error_handler(&comm, &err);
  };
  // Unfortunately cannot be inline because of .h loop
  void Handler(OOMPI_Comm *comm, int err);

  // Table maintenance

  void Add(MPI_Comm comm, OOMPI_Comm *oompi, OOMPI_Error_action action);
  void Change(MPI_Comm comm, OOMPI_Comm *oompi);
  void Change(MPI_Comm comm, OOMPI_Error_action action);
  void Delete(MPI_Comm comm);
  OOMPI_Comm *Get_oompi(MPI_Comm comm);
  OOMPI_Error_action Get_action(MPI_Comm comm);

protected:
  static bool init;
  bool valid;

  OOMPI_Error_entry *root;

private:
  // Private helper function

  OOMPI_Error_entry *get_entry(MPI_Comm comm);

  // Stupid function to avoid compiler warnings

  inline void do_nothing(void) const {};
};

#endif
