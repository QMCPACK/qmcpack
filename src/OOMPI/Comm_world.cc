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
// MPI_COMM_WORLD class
//

#include <Configuration.h>
#include "oompi-config.h"

#include "Comm_world.h"
#include "Util.h"
#include "Error.h"
#include "Error_table.h"
#include "Message.h"
#include "Datatype.h"
#include <Message/OpenMP.h>

//
// Instantiate OOMPI_COMM_WORLD
//

OOMPI_Comm_world OOMPI_COMM_WORLD;


//
// Instantiate and initialize the static variables
//

bool OOMPI_Comm_world::created = false;
bool OOMPI_Comm_world::oompi_finalized = false;


//
// Default constructor
//
OOMPI_Comm_world::OOMPI_Comm_world(void) 
: OOMPI_Intra_comm(MPI_COMM_NULL), valid(false)
{
  if (created)
    OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_COMM);
  else
    valid = true;
  
  created = true;
}


//
// Copy constructor
//
OOMPI_Comm_world::OOMPI_Comm_world(const OOMPI_Comm_world &)
: OOMPI_Intra_comm(MPI_COMM_NULL), valid(false)
{
  OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_COMM);
}


//
// Assignment operator
//
OOMPI_Comm_world &
OOMPI_Comm_world::operator=(const OOMPI_Comm_world &)
{
  OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_COMM);
  return *this;
}


//
// Destructor 
//
OOMPI_Comm_world::~OOMPI_Comm_world(void)
{
  // Do not call MPI_Finalize here.  It is the user's responsibility to
  // call OOMPI_COMM_WORLD.Finalize()

  // Calling it here could result in calling MPI_Finalize() being called
  // after MPI has been shut down! (i.e. this destructor may be called 
  // after other various onexit()/atexit() functions, or other destructors)
}


//
// MPI_COMM_world functions
// MPI_Init
//
void
OOMPI_Comm_world::Init(int& argc, char**& argv)
{
  Init(argc, argv, true);
}


//
// Initializes OOMPI.  This assumes that MPI_Init has already been called!
//
void
OOMPI_Comm_world::Init(void)
{
  int argc;
  char **argv;

  Init(argc, argv, false);
}


//
// The "real" Init function
//
void
OOMPI_Comm_world::Init(int& argc, char**& argv, bool call_init)
{
  int flag(0);

#if 0
  // Not supported at present.  Some MPI's deadlock here (SGI).
  // This will require re-thinking to make thread safe...
  OOMPI_Util util;
  util.Get_sem(SEM_INIT);
#endif

  if (valid) {

    // Make it ok to call the OOMPI Init even deep into user code, long
    // after MPI_Init has been called.

    if (call_init) {
      MPI_Initialized(&flag);
      if (!flag) 
      {
#if defined(ENABLE_OPENMP)
    int provided, claimed;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Query_thread(&claimed);
    if (claimed != provided) 
    {
      std::ostringstream o;
  o << "OOMPI_Comm_world::init"
	<< "\n  MPI_Query_thread thread level " << claimed
    << "\n  MPI_Init_threadthread level " << provided;
      APP_ABORT(o.str());
    }
#else
MPI_Init(&argc, &argv);
#endif
      }
    }

    MPI_Constructor(MPI_COMM_WORLD, false);
    Set_error_action(OOMPI_ERRORS_ARE_FATAL);

#if OOMPI_HAVE_ANSI_COMPLEX || OOMPI_HAVE_LONG_DOUBLE || OOMPI_HAVE_LONG_LONG_INT
    Message a(MPI_DATATYPE_NULL, 0);
    a.Init();
#endif

    int flag2, *addr;

    // Initialize other global instances from built in attributes

    MPI_Attr_get(MPI_COMM_WORLD, MPI_HOST, &addr, &flag2);
    OOMPI_HOST = *addr;
    if (!flag2)
      OOMPI_HOST = OOMPI_PROC_NULL;

    MPI_Attr_get(MPI_COMM_WORLD, MPI_IO, &addr, &flag2);
    OOMPI_IO = *addr;
    if (!flag2)
      OOMPI_IO = OOMPI_PROC_NULL;

    MPI_Attr_get(MPI_COMM_WORLD, MPI_WTIME_IS_GLOBAL, &addr, &flag2);
    OOMPI_WTIME_IS_GLOBAL = (bool) *addr;
    if (!flag2)
      OOMPI_WTIME_IS_GLOBAL = false;
  }

#if 0
  // See above note on deadlocking/thread safe
  util.Release_sem(SEM_INIT);
#endif
}


//
// MPI_Finalize
//
void
OOMPI_Comm_world::Finalize(void)
{
  int flag;
  OOMPI_Util a;

  a.Get_sem(SEM_FINAL);
  MPI_Initialized(&flag);

  if (valid && flag && !oompi_finalized) {
    // Free all remaining MPI datatypes.
    MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
    OOMPI_Datatype::Free_all_datatypes();

    oompi_finalized = true;

    MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
    if (err_handler != MPI_ERRHANDLER_NULL)
      MPI_Errhandler_free(&err_handler);

    MPI_Finalize();
  }
  a.Release_sem(SEM_FINAL);
}


//
// Virtual function from abstract base class
//
bool
OOMPI_Comm_world::Test_inter(void)
{
  return false;
}

