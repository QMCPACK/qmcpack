// Copyright 2018-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

#include <boost/core/lightweight_test.hpp>
#include <iostream>

namespace mpi3 = boost::mpi3;

using std::cout;

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	BOOST_TEST(world.size() == 6);

	mpi3::communicator const fifth = world / 5;

	cout << "I am rank " << world.rank() << " in " << world.name() << ", ";

	if(fifth) {
		cout << "I am also   " << fifth.rank() << " in " << fifth.name() << '\n';
	} else {
		cout << "I am not in " << fifth.name() << '\n';
	}

	return boost::report_errors();
} catch(...) {
	return 1;
}
