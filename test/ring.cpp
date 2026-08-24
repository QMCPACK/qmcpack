// Copyright 2021-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/process.hpp>

#include <boost/core/lightweight_test.hpp>

namespace mpi3 = boost::mpi3;

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	auto const size = world.size();
	BOOST_TEST(size != 0);

	mpi3::process&& next  = world[(world.rank() + size + 1) % size];
	mpi3::process&& prior = world[(world.rank() + size - 1) % size];

	int token;  // NOLINT(cppcoreguidelines-init-variables)
	if(!world.is_root()) {
		prior >> token;
	} else {
		token = -1;
	}

	next << token;

	if(world.is_root()) {
		prior >> token;
	}

	BOOST_TEST(token == -1);

	return boost::report_errors();
} catch(...) {
	return 1;
}
