// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2023 Alfredo A. Correa

#ifndef BOOST_MPI3_MAIN_ENVIRONMENT_HPP
#define BOOST_MPI3_MAIN_ENVIRONMENT_HPP

#pragma once

#ifdef BOOST_MPI3_MAIN_HPP
#error Include either "mpi3/main.hpp" or "mpi3/main_environment.hpp"
#endif

#include <mpi.h>

#include "../mpi3/communicator.hpp"
#include "../mpi3/environment.hpp"
#include "../mpi3/exception.hpp"

namespace boost {
namespace mpi3 {

static int main(int /*argc*/, char** /*argv*/, boost::mpi3::environment& /*env*/);  // NOLINT(bugprone-exception-escape)

}  // end namespace mpi3
}  // end namespace boost

int main(int argc, char* argv[]) {  // NOLINT(bugprone-exception-escape,misc-definitions-in-headers) main defined in a header file for replacement
	boost::mpi3::environment env{argc, argv};
	return boost::mpi3::main(argc, argv, env);
}
#endif
