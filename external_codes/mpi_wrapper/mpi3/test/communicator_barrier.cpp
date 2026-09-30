// Copyright 2018-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

#include <iostream>

namespace mpi3 = boost::mpi3;
using std::cout;

int main(int argc, char** argv) {  // NOLINT(bugprone-exception-escape)
	mpi3::environment env(argc, argv);

	auto world = env.world();

	cout << "Before barrier, I am " << world.rank() << " of " << world.size() << '\n';
	world.barrier();
	cout << "After barrier, I am " << world.rank() << " of " << world.size() << '\n';
}
