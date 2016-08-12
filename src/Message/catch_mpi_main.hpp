
#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#include "Message/Communicate.h"

// Replacement unit test main function to ensure that MPI is finalized once 
// (and only once) at the end of the unit test.

int main(int argc, const char* argv[])
{
  int result = Catch::Session().run(argc, argv);
  OHMMS::Controller->finalize();
  return result;
}

