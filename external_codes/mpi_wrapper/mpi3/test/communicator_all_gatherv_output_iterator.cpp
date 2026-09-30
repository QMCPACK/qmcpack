// Copyright 2018-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

#include <boost/core/lightweight_test.hpp>
#include <iterator>
#include <vector>

namespace mpi3 = boost::mpi3;

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	using T = double;
	BOOST_TEST(world.size() > 2);

	// initialize local data ///////////////////////////////////////////////////////
	auto v_loc = [&] {
		switch(world.rank()) {
		case 0 : return std::vector<T>{0.0, 0.0, 0.0};
		case 1 : return std::vector<T>{1.0, 1.0, 1.0, 1.0};
		case 2 : return std::vector<T>{2.0, 2.0, 2.0, 2.0, 2.0};
		default: return std::vector<T>{};
		}
	}();

	// gather communication ////////////////////////////////////////////////////////
	// doesn't work on NVHPC 24.5
	// std::vector<T> v;
	// v.reserve(v_local.size()*world.size()); // optional! avoid a few allocations
	// world.all_gatherv(begin(v_loc), end(v_loc), std::back_inserter(v));
	// v.shrink_to_fit();                      // optional! save a few memory

	std::vector<T> v(12);
	world.all_gatherv(begin(v_loc), end(v_loc), v.begin());

	// check communication /////////////////////////////////////////////////////////
	BOOST_TEST((v == std::vector<T>{0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0}));

	return boost::report_errors();
} catch(...) {
	return 1;
}
