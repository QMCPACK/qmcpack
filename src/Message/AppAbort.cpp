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

#include "Message/AppAbort.h"
#include <iostream>

#ifdef HAVE_MPI
#include "mpi3/environment.hpp"

void breakableAppAbort(const std::string& str_msg)
{
  std::cerr << str_msg << '\n';
  MPI_Abort(MPI_COMM_WORLD, 1);
}
#else
void breakableAppAbort(const std::string& str_msg)
{
  std::cerr << str_msg << '\n';
  exit(1);
}
#endif
