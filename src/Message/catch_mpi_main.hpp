//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#include "Message/Communicate.h"

// Replacement unit test main function to ensure that MPI is finalized once 
// (and only once) at the end of the unit test.

int main(int argc, char* argv[])
{
#ifdef HAVE_MPI
  mpi3::environment env(argc, argv);
  OHMMS::Controller->initialize(env);
#endif
  int result = Catch::Session().run(argc, argv);
  OHMMS::Controller->finalize();
  return result;
}

