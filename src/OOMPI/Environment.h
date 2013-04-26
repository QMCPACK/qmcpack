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

#ifndef _OOMPI_ENVIRONMENT_H_
#define _OOMPI_ENVIRONMENT_H_

#include <mpi.h>
#include "oompi-config.h"
#include "Error.h"
#include "Comm_world.h"


//
// Class declaraction
//

class OOMPI_Environment
{
public:
  OOMPI_Environment(void);
  OOMPI_Environment(const OOMPI_Environment &a)
    : valid(false), alloc(false), buffer(0)
  {
    a.do_nothing();
    OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_OTHER);
  };
  OOMPI_Environment &operator=(const OOMPI_Environment &a)
  {
    buffer = 0;
    valid = false;
    alloc = false;
    a.do_nothing();
    OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_OTHER);
    return *this;
  };
  inline ~OOMPI_Environment() {};

  // Access functions
  // Buffer

  void Buffer_attach(int size);
  void Buffer_attach(void *buf, int size);
  int Buffer_detach(void);

  // Processor dependant stuff

  char *Get_processor_name(void);
  char *Get_processor_name(int &len);
  void Get_processor_name(char name[], int &len);
  inline double Wtick(void)
  {
    return MPI_Wtick();
  };
  inline double Wtime(void)
  {
    return MPI_Wtime();
  };

  // Profiling

  int Pcontrol(int level, ...);

protected:
  static bool initialized;
  bool valid;
  bool alloc;

  char *buffer;
  int buf_size;

private:
  // Stupid function to avoid compiler warnings
  inline void do_nothing(void) const {} ;
};

#endif
