// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2023 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/main.hpp>
#include <mpi3/process.hpp>

#include <boost/variant.hpp>

#include <complex>
#include <list>
#include <string>

namespace mpi3 = boost::mpi3;

namespace error_ns {
enum code { FAILURE = 0, SUCCESS = 1, RECOVERABLE = 2 };
}  // end namespace error_ns

namespace error_long_ns {
enum code : unsigned long { FAILURE = 0, SUCCESS = 1, RECOVERABLE = 2 };  // NOLINT(google-runtime-int) test unsiged long specifically
}  // end namespace error_long_ns

namespace error_class_ns {
enum class code { FAILURE = 0, SUCCESS = 1, RECOVERABLE = 2 };
}  // end namespace error_class_ns

auto mpi3::main(int /*argc*/, char** /*argv*/, mpi3::communicator world) -> int try {
	assert(world.size() > 1);

	switch(world.rank()) {
	case 0: {
		error_ns::code const ec = error_ns::SUCCESS;
		world.send_n(&reinterpret_cast<int const&>(ec), 1, 1);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) testing not recommended brute force method
	}; break;
	case 1: {
		error_ns::code ec{};
		world.receive_n(&reinterpret_cast<int&>(ec), 1, 0);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) testing not recommended brute force method
		assert(ec == error_ns::SUCCESS);
	};
	}

	switch(world.rank()) {
	case 0: {
		error_ns::code const ec = error_ns::SUCCESS;
		world.send_n(&ec, 1, 1);
	};
	break; case 1: {
		error_ns::code ec = error_ns::FAILURE;
		world.receive_n(&ec, 1, 0);
		assert(ec == error_ns::SUCCESS);
	};
	}

	switch(world.rank()) {
	case 0: {
		error_ns::code const ec = error_ns::SUCCESS;
		world.send_n(&ec, 1, 1);
	};
	break; case 1: {
		error_ns::code ec = error_ns::FAILURE;
		world.receive_n(&ec, 1);
		assert(ec == error_ns::SUCCESS);
	};
	}

	switch(world.rank()) {
	case 0: {
		error_ns::code const ec = error_ns::SUCCESS;
		world[1] << ec;
	};
	break; case 1: {
		error_ns::code ec = error_ns::FAILURE;
		world[0] >> ec;
		assert(ec == error_ns::SUCCESS);
	};
	}

	switch(world.rank()) {
	case 0: {
		error_ns::code const ec = error_ns::SUCCESS;
		world[1] << ec;
	};
	break; case 1: {
		error_ns::code ec = error_ns::FAILURE;
		world >> ec;
		assert(ec == error_ns::SUCCESS);
	};
	}

	switch(world.rank()) {
	case 0: {
		error_long_ns::code const ec = error_long_ns::SUCCESS;
		world[1] << ec;
	};
	break; case 1: {
		error_long_ns::code ec  = error_long_ns::FAILURE;
		world >> ec;
		assert(ec == error_long_ns::SUCCESS);
	};
	}

	switch(world.rank()) {
	case 0: {
		error_class_ns::code const ec = error_class_ns::code::SUCCESS;
		world[1] << ec;
	};
	break; case 1: {
		error_class_ns::code ec = error_class_ns::code::FAILURE;
		world >> ec;
		assert(ec == error_class_ns::code::SUCCESS);
	};
	}

	return 0;
} catch(...) {
	return 1;
}
