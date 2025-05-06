////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// Refactored from: Communicate.cpp
////////////////////////////////////////////////////////////////////////////////

#include "AppAbort.h"
#include <iostream>
#include "config.h"

#ifdef ENABLE_GCOV
#ifdef __GNUC__
extern "C" void __gcov_dump();
#endif
#else
#endif

#ifdef HAVE_MPI
#include <mpi.h>

void breakableAppAbort(const std::string& str_msg)
{
  std::cerr << str_msg << '\n';
#ifdef ENABLE_GCOV
#ifdef __GNUC__
  __gcov_dump();
#endif
#endif
  MPI_Abort(MPI_COMM_WORLD, 1);
}
#else
void breakableAppAbort(const std::string& str_msg)
{
  std::cerr << str_msg << '\n';
#ifdef ENABLE_GCOV
#ifdef __GNUC__
  __gcov_dump();
#endif
#endif
  exit(1);
}
#endif
