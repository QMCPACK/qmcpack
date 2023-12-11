// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#include "../../mpi3/communicator.hpp"
#include "../../mpi3/main.hpp"
#include "../../mpi3/process.hpp"

namespace mpi3 = boost::mpi3;

void part1(mpi3::communicator& world) {
	std::size_t const count = 120;
	std::vector<int> send_buffer(count);
	iota(send_buffer.begin(), send_buffer.end(), 0);

	std::vector<int> recv_buffer;
	if(world.rank() == 0) {
		recv_buffer.resize(count, -1);
	}

	world.reduce_n(send_buffer.begin(), send_buffer.size(), recv_buffer.begin(), std::plus<>{}, 0);
	if(world.rank() == 0) {
		for(std::size_t i = 0; i != recv_buffer.size(); ++i) {  // NOLINT(altera-unroll-loops) use algorithm
			assert(std::size_t(recv_buffer[i]) == i * static_cast<std::size_t>(world.size()));
		}
	}
}

void part2(mpi3::communicator& world) {
	double const v     = world.rank();
	double       total = 0;

	double const* const it = world.reduce_n(&v, 1, &total, std::plus<>{});
	if(world.rank() == 0) {
		assert(total == static_cast<double>(world.size() * (world.size() - 1)) / 2);
	} else {
		assert(total == 0);
	}
	if(world.rank() == 0) {
		assert(it != &total);
	} else {
		assert(it == &total);
	}
}

void part3(mpi3::communicator& world) {
	mpi3::optional<int> total = (world[0] += world.rank());
	if(world.rank() == 0) {
		assert(total);
	}
	if(world.rank() != 0) {
		assert(not total);
	}
	if(total) {
		assert(*total == static_cast<double>(world.size() * (world.size() - 1)) / 2);
	}
}

int mpi3::main(int /*argc*/, char** /*argv*/, mpi3::communicator world) try {
	assert(world.size() > 1);

	part1(world);
	part2(world);  // TODO(correaa) fix this
	part3(world);

	return 0;
} catch(...) { return 1; }
