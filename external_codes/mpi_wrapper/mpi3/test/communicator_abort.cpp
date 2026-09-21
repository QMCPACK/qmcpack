// Copyright 2018-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

namespace mpi3 = boost::mpi3;

int main(int argc, char** argv) {  // NOLINT(bugprone-exception-escape)
	mpi3::environment env(argc, argv);

	auto world = env.world();
	world.abort(911);

	return 0;
}
