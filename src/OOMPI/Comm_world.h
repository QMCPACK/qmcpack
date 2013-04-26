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


#ifndef _OOMPI_COMM_WORLD_H_
#define _OOMPI_COMM_WORLD_H_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "oompi-config.h"
#include "Intra_comm.h"


class OOMPI_Comm_world : public OOMPI_Intra_comm
{
public:
  //
  // Big 4
  //

  OOMPI_Comm_world(void);
  OOMPI_Comm_world(const OOMPI_Comm_world& a);
  OOMPI_Comm_world& operator=(const OOMPI_Comm_world& a);
  ~OOMPI_Comm_world(void);

  //
  // MPI_COMM_WORLD functions
  //

  // Calls MPI_Init, then initializes OOMPI
  void Init(int& argc, char**& argv);
  // Initializes OOMPI -- assumes that MPI_Init has been called already
  void Init(void);

  // returns the return of MPI_Finalize
  void Finalize(void);

  // Return if we are an intercommunicator or not
  virtual bool Test_inter(void);

  // Return whether we have been finalized or not
  inline bool Finalized(void)
  {
    return oompi_finalized;
  };

protected:
  static bool created;
  static bool oompi_finalized;

  bool valid;

private:
  // The "real" Init function
  void Init(int& argc, char**& argv, bool call_init);

};


//
// OOMPI_COMM_WORLD
//

extern OOMPI_Comm_world OOMPI_COMM_WORLD;


#endif
