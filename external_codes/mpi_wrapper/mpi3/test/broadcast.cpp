// Copyright 2023-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>
// #include <multi/array.hpp>

#include <algorithm>
#include <boost/core/lightweight_test.hpp>
#include <vector>

namespace mpi3 = boost::mpi3;

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();
	{
		auto v = std::vector<double>(10, world.root() ? 3.14 : 99.9);

		auto it = world.broadcast_n(v.begin(), 10);

		BOOST_TEST(it == v.end());
		BOOST_TEST(std::find_if(v.cbegin(), v.cend(), [](auto const& e) { return e != 3.14; }) == v.cend());  // NOLINT(boost-use-ranges) for C++20, use std::ranges::find_if
	}
	// {
	// 	namespace multi = boost::multi;
	// 	auto arr        = multi::array<double, 1>(10, world.root() ? 3.14 : 99.9);

	// 	auto it = world.broadcast_n(arr.begin(), 10);

	// 	BOOST_TEST(it == arr.end());
	// 	BOOST_TEST(std::find_if(arr.cbegin(), arr.cend(), [](auto const& e) { return e != 3.14; }) == arr.cend());  // NOLINT(boost-use-ranges) for C++20, use std::ranges::find_if
	// }

	return boost::report_errors();
} catch(...) {
	return 1;
}
