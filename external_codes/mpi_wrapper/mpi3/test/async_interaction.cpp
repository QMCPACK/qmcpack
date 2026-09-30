// Copyright 2018-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/process.hpp>

#include <boost/core/lightweight_test.hpp>
#include <future>
#include <iostream>

namespace mpi3 = boost::mpi3;

namespace {

auto async_send(mpi3::communicator& comm, int val, int target) {  // NOLINT(bugprone-easily-swappable-parameters)
	return std::async([=, &comm]() {                              // was: [=] () mutable {  in original Joseph's S. posting
		auto value = val + 1;
		comm[target] << value;  // same as comm.send(&value, &value + 1, target);
	});
}

}  // end namespace

auto main(int argc, char** argv) -> int try {
	boost::mpi3::environment env(argc, argv);

	auto world = env.world();

	if(world.rank() == 0) {
		auto fut = async_send(world, 41, 1);
		fut.wait();
	}

	if(world.rank() == 1) {
		int value{};
		world[0] >> value;  // same as world.receive(&value, &value + 1, 0);
		std::cout << "value = " << value << '\n';
		BOOST_TEST(value == 42);
	}

	return boost::report_errors();
} catch(...) {
	return 1;
}
