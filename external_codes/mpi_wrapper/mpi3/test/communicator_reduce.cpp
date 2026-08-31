// Copyright 2018-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/process.hpp>

#include <boost/core/lightweight_test.hpp>
#include <cstddef>
#include <functional>
#include <numeric>
#include <vector>

namespace mpi3 = boost::mpi3;

namespace {

void part1(mpi3::communicator& world) {
	std::size_t const count = 120;
	std::vector<int> send_buffer(count);
	iota(send_buffer.begin(), send_buffer.end(), 0);  // NOLINT(boost-use-ranges) for C++20, use std::ranges::iota

	std::vector<int> recv_buffer;
	if(world.rank() == 0) {
		recv_buffer.resize(count, -1);
	}

	world.reduce_n(send_buffer.begin(), send_buffer.size(), recv_buffer.begin(), std::plus<>{}, 0);
	if(world.rank() == 0) {
		for(std::size_t i = 0; i != recv_buffer.size(); ++i) {  // NOLINT(altera-unroll-loops) use algorithm
			BOOST_TEST(std::size_t(recv_buffer[i]) == i * static_cast<std::size_t>(world.size()));
		}
	}
}

void part2(mpi3::communicator& world) {
	double const v     = world.rank();
	double       total = 0;

	double const* const it = world.reduce_n(&v, 1, &total, std::plus<>{}); (void)it;
	if(world.rank() == 0) {
		BOOST_TEST(total == static_cast<double>(world.size() * (world.size() - 1)) / 2);
	} else {
		BOOST_TEST(total == 0);
	}
	if(world.rank() == 0) {
		BOOST_TEST(it != &total);
	} else {
		BOOST_TEST(it == &total);
	}
}

void part3(mpi3::communicator& world) {
	mpi3::optional<int> total = (world[0] += world.rank());
	if(world.rank() == 0) {
		BOOST_TEST(total);
	}
	if(world.rank() != 0) {
		BOOST_TEST(not total);
	}
	if(total) {
		BOOST_TEST(*total == static_cast<double>(world.size() * (world.size() - 1)) / 2);
	}
}

}  // end namespace

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	BOOST_TEST(world.size() > 1);

	part1(world);
	part2(world);  // TODO(correaa) fix this
	part3(world);

	return boost::report_errors();
} catch(...) {
	return 1;
}
