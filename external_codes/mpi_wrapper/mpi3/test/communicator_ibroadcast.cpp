// Copyright 2020-2024 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

#ifndef EXAMPI
#include <mpi3/ostream.hpp>
#endif

#include <algorithm>
#include <boost/core/lightweight_test.hpp>
#include <iostream>
#include <iterator>
#include <numeric>
#include <thread>
#include <vector>

namespace mpi3 = boost::mpi3;

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	BOOST_TEST(world.size() > 2);

	std::vector<int> large(10);
	if(world.root()) {
		iota(large.begin(), large.end(), 0);  // NOLINT(boost-use-ranges) for C++20, use std::ranges::iota
	}

#ifndef EXAMPI
	mpi3::ostream wout(world);
	wout << "before:\n" << std::flush;
	std::copy(large.begin(), large.end(), std::ostream_iterator<int>(wout, " "));  // NOLINT(boost-use-ranges) for C++20, use std::ranges::copy

	wout << '\n' << std::flush;

	{
		auto req = world.ibroadcast(large.begin(), large.end(), 0);
		using namespace std::chrono_literals;
		std::this_thread::sleep_for(5s);  // NOLINT(misc-include-cleaner) bug in clang-tidy 18.1
	}

	wout << "after:\n" << std::flush;
	std::copy(large.begin(), large.end(), std::ostream_iterator<int>(wout, " "));  // NOLINT(boost-use-ranges) for C++20, use std::ranges::copy
	wout << '\n' << std::flush;
#endif

	return boost::report_errors();
} catch(...) {
	return 1;
}
