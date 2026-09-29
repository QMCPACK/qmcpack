// Copyright 2018-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

#include <mpi3/process.hpp>
#include <mpi3/vector.hpp>

#include <algorithm>
#include <boost/core/lightweight_test.hpp>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <list>
#include <string>
#include <utility>
#include <vector>

namespace mpi3 = boost::mpi3;

namespace {

void part1(mpi3::communicator& world) {
	std::vector<std::pair<double, int>> small(10, {0.0, world.rank()});
	std::vector<std::pair<double, int>> large(world.root() ? small.size() * static_cast<std::size_t>(world.size()) : 0, std::pair<double, int>(0.0, -10));

	auto it = world.gather_n(small.begin(), small.size(), large.begin(), 0);
	BOOST_TEST(it == large.end());

	if(world.rank() == 0) {
		BOOST_TEST(it != large.begin());
		BOOST_TEST((large[9] == std::pair<double, int>(0.0, 0)));
		BOOST_TEST((large[11] == std::pair<double, int>(0.0, 1)));
	} else {
		BOOST_TEST(it == large.begin());
	}
}

void part2(mpi3::communicator& world) {
	std::vector<std::pair<double, int>> small(10, {0., world.rank()});
	std::vector<std::pair<double, int>> large(world.root() ? small.size() * static_cast<std::size_t>(world.size()) : 0, std::pair<double, int>(0., -1));

	auto it = world.gather_n(small.begin(), small.size(), large.begin());
	BOOST_TEST(it == large.end());

	if(world.root()) {
		BOOST_TEST(it != large.begin());
		BOOST_TEST((large[9] == std::pair<double, int>(0., 0)));
		BOOST_TEST((large[11] == std::pair<double, int>(0., 1)));
	} else {
		BOOST_TEST(it == large.begin());
	}
}

void part3(mpi3::communicator& world) {
	std::list<double> small(10, world.rank());
	std::vector<double> large(world.root()?small.size()*static_cast<std::size_t>(world.size()):0, -1.0);

	world.gather(small.begin(), small.end(), large.begin(), 0);
	if(world.root()) {
		std::cout << "large: ";
		std::copy(large.begin(), large.end(), std::ostream_iterator<double>(std::cout, ", "));  // NOLINT(boost-use-ranges) for C++20, use std::ranges::copy
		std::cout << '\n';

		BOOST_TEST(large[ 1] == 0);
		BOOST_TEST(large[11] == 1);
		BOOST_TEST(large[21] == 2);
	}
}

// void part4(mpi3::communicator& world) {
//  auto val = std::string{"5.1 opa"};

//  using T = decltype(val);

//  std::vector<T> small(10, val);
//  std::vector<T> large(world.root() ? small.size() * static_cast<std::size_t>(world.size()) : 0);

//  world.gather(small.begin(), small.end(), large.begin(), 0);

//  if(world.rank() == 0) {
//      BOOST_TEST(all_of(large.begin(), large.end(), [val](auto& e) { return val == e; }));
//  }
// }

#ifndef EXAMPI
void part5(mpi3::communicator& world) {
	auto Lval  = std::to_string(world.rank() + 1000);
	auto vals0 = (world[0] |= Lval);
	if(world.rank() == 0) {
		BOOST_TEST(vals0.size() - static_cast<std::size_t>(world.size()) == 0);
		BOOST_TEST(vals0[2] == "1002");
	} else {
		BOOST_TEST(vals0.empty());
	}
}
#endif

}  // end namespace 

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	mpi3::vector<double> const v;

	BOOST_TEST( world.size() > 2);

	part1(world);
	part2(world);
	part3(world);
//  part4(world);
#ifndef EXAMPI
	part5(world);
#endif

	return boost::report_errors();
} catch(...) {
	return 1;
}
