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
// Environment class
//

#include <stdarg.h>
#include "Environment.h"
#include "Constants.h"
#include "Comm_world.h"


//
// Global instance
//

OOMPI_Environment OOMPI_ENV;


//
// Class variable
//

bool OOMPI_Environment::initialized = false;


//
// Main constructor
//
OOMPI_Environment::OOMPI_Environment()
: alloc(false), buffer(0)
{
  OOMPI_Util a;

  a.Get_sem(SEM_ENVIRONMENT);
  if (!initialized) {
    initialized = true;
    valid = true;
    a.Release_sem(SEM_ENVIRONMENT);
  }
  else {
    valid = false;
    a.Release_sem(SEM_ENVIRONMENT);
    OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_OTHER);
  }
}


//
// Buffer attaches
//
void 
OOMPI_Environment::Buffer_attach(int size)
{
  if (!valid) {
    OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_OTHER);
    return;
  }

  OOMPI_Util a;
  a.Get_sem(SEM_BUFFER);

  if (buffer != 0) {
    a.Release_sem(SEM_BUFFER);
    OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_OTHER);
    return;
  }

  buffer = new char[size];
  alloc = true;

  if (buffer == 0) {
    a.Release_sem(SEM_BUFFER);
    OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_OTHER);
  }
  else {
    buf_size = size;
    a.Release_sem(SEM_BUFFER);
    if (MPI_Buffer_attach(buffer, size) != MPI_SUCCESS) {
      delete[] buffer;
      buffer = 0;
      size = 0;
    }
  }
}


void 
OOMPI_Environment::Buffer_attach(void *buf, int size)
{ 
  if (!valid)
    OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_OTHER);

  OOMPI_Util a;

  a.Get_sem(SEM_BUFFER);
  buffer = (char*) buf;
  alloc = false;
  buf_size = size;
  a.Release_sem(SEM_BUFFER);

  if (MPI_Buffer_attach(buf, size) != MPI_SUCCESS) {
    buffer = 0;
    buf_size = 0;
  }
}


//
// Buffer detach
//
int 
OOMPI_Environment::Buffer_detach()
{
  if (!valid)
    OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_OTHER);

  OOMPI_Util a;
  a.Get_sem(SEM_BUFFER);

  void *temp = buffer;
  int dummy = 0;

  buffer = 0;
  buf_size = 0;
  a.Release_sem(SEM_BUFFER);

  MPI_Buffer_detach(temp, &dummy);

  if( alloc ) {
    delete[] (char*)temp;
  }

  return dummy;
}


//
// Processor dependant stuff
//
char *
OOMPI_Environment::Get_processor_name()
{
  if (!valid)
    OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_OTHER);

  char *name = new char[MPI_MAX_PROCESSOR_NAME];
  int len;

  if (MPI_Get_processor_name(name, &len) == MPI_SUCCESS)
    return name;

  return 0;
}


char *
OOMPI_Environment::Get_processor_name(int &len)
{
  if (!valid)
    OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_OTHER);

  char *name = new char[MPI_MAX_PROCESSOR_NAME];

  if (MPI_Get_processor_name(name, &len) == MPI_SUCCESS)
    return name;

  return 0;
}


void
OOMPI_Environment::Get_processor_name(char name[], int &len)
{
  if (!valid)
    OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_OTHER);

  MPI_Get_processor_name(name, &len);
}


//
// Profiling
//
int
OOMPI_Environment::Pcontrol(int level, ...)
{
  if (!valid)
    OOMPI_ERROR.Handler(&OOMPI_COMM_WORLD, OOMPI_ERR_OTHER);

  va_list ap;
  int ret;

  va_start(ap, level);
  ret = MPI_Pcontrol(level, ap);
  va_end(ap);
  
  return ret;
}


