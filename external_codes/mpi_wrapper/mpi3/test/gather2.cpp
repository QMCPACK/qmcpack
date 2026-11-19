// Copyright 2022-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

#include <boost/core/lightweight_test.hpp>
#include <cstddef>
#include <iterator>
#include <tuple>
#include <vector>

namespace mpi3 = boost::mpi3;

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	using dd = std::tuple<double, double>;

	std::vector<dd> v_local(10, dd{world.rank(), world.rank() + 1});
	std::vector<dd> v(world.root() ? v_local.size() * static_cast<std::size_t>(world.size()) : 0);

	auto last = world.gather(begin(v_local), end(v_local), begin(v));

	if(world.root()) {
		BOOST_TEST(last == end(v));
		BOOST_TEST((v[0] == dd{0., 1.}));
		BOOST_TEST((v[10] == dd{1., 2.}));
		BOOST_TEST((v[20] == dd{2., 3.}));
	} else {
		BOOST_TEST(last == begin(v));
	}

	return boost::report_errors();
} catch(...) {
	return 1;
}
