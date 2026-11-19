// Copyright 2023-2025 Alfredo A. Correa

// based on http://mpi.deino.net/mpi_functions/MPI_Alltoall.html

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

#ifndef EXAMPI
#include <mpi3/ostream.hpp>
#endif

#include <algorithm>
#include <boost/core/lightweight_test.hpp>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <numeric>
#include <vector>

namespace mpi3 = boost::mpi3;

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	std::size_t const chunk = 5;

	auto sb = std::vector<int>(static_cast<std::size_t>(world.size()) * chunk);
	std::iota(sb.begin(), sb.end(), 40000 + (world.rank() * 100));  // NOLINT(boost-use-ranges) for C++20, use std::ranges::iota

	auto rb = std::vector<int>(static_cast<std::size_t>(world.size()) * chunk);

	auto sz = static_cast<std::size_t>(world.size());
	BOOST_TEST( sz != 0 );
	BOOST_TEST( sb.size() % sz == 0);

	world.all_to_all_n(sb.data(), sb.size() / sz, rb.data());

#ifndef EXAMPI
	mpi3::ostream wout(world);
	std::copy(sb.begin(), sb.end(), std::ostream_iterator<int>(wout << "sb = ", ", "));  // NOLINT(boost-use-ranges) for C++20, use std::ranges::copy
	wout << '\n'
	     << std::flush;
	std::copy(rb.begin(), rb.end(), std::ostream_iterator<int>(wout << "rb = ", ", "));  // NOLINT(boost-use-ranges) for C++20, use std::ranges::copy
	wout << '\n'
	     << std::flush;

	world.all_to_all_inplace_n(sb.data(), sb.size() / sz);
	// do_all_to_all_n(world, sb.data(), sb.size(), sb.data());
	// world.all_to_all_n(sb.data(), sb.size()); //  , sb.data());
	// std::copy(sb.begin(), sb.end(), std::ostream_iterator<int>(wout<<"sb (inplace) = ", ", ")); wout<<std::endl;
	BOOST_TEST(sb == rb);
#endif

	return boost::report_errors();
} catch(...) {
	return 1;
}
