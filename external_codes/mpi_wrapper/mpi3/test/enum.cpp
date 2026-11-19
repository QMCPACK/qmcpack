// Copyright 2023-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/process.hpp>

#include <boost/core/lightweight_test.hpp>

namespace mpi3 = boost::mpi3;

namespace error_ns {
	enum code { FAILURE = 0, SUCCESS = 1, RECOVERABLE = 2 };  // NOLINT(cppcoreguidelines-use-enum-class,performance-enum-size)
}  // end namespace error_ns

namespace error_long_ns {
	enum code : unsigned long { FAILURE = 0, SUCCESS = 1, RECOVERABLE = 2 };  // NOLINT(cppcoreguidelines-use-enum-class,google-runtime-int,performance-enum-size) test unsigned long specifically
}  // end namespace error_long_ns

namespace error_class_ns {
	enum class code { FAILURE = 0, SUCCESS = 1, RECOVERABLE = 2 };  // NOLINT(performance-enum-size)
}  // end namespace error_class_ns

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	BOOST_TEST(world.size() > 1);

	switch(world.rank()) {
	case 0: {
		error_ns::code const ec = error_ns::SUCCESS;
		world.send_n(&ec, 1, 1);
	}; break;
	case 1: {
		error_ns::code ec{};
		world.receive_n(&ec, 1, 0);
		BOOST_TEST(ec == error_ns::SUCCESS);
	}; break;
	default: {}
	}

	switch(world.rank()) {
	case 0: {
		error_ns::code const ec = error_ns::SUCCESS;
		world.send_n(&ec, 1, 1);
	};
	break; case 1: {
		error_ns::code ec = error_ns::FAILURE;
		world.receive_n(&ec, 1, 0);
		BOOST_TEST(ec == error_ns::SUCCESS);
	};
	break; default: {}
	}

	switch(world.rank()) {
	case 0: {
		error_ns::code const ec = error_ns::SUCCESS;
		world.send_n(&ec, 1, 1);
	};
	break; case 1: {
		error_ns::code ec = error_ns::FAILURE;
		world.receive_n(&ec, 1);
		BOOST_TEST(ec == error_ns::SUCCESS);
	};
	break; default: {}
	}

	switch(world.rank()) {
	case 0: {
		error_ns::code const ec = error_ns::SUCCESS;
		world[1] << ec;
	};
	break; case 1: {
		error_ns::code ec = error_ns::FAILURE;
		world[0] >> ec;
		BOOST_TEST(ec == error_ns::SUCCESS);
	};
	break; default: {}
	}

	switch(world.rank()) {
	case 0: {
		error_ns::code const ec = error_ns::SUCCESS;
		world[1] << ec;
	};
	break; case 1: {
		error_ns::code ec = error_ns::FAILURE;
		world >> ec;
		BOOST_TEST(ec == error_ns::SUCCESS);
	};
	break; default: {}
	}

	switch(world.rank()) {
	case 0: {
		error_long_ns::code const ec = error_long_ns::SUCCESS;
		world[1] << ec;
	};
	break; case 1: {
		error_long_ns::code ec  = error_long_ns::FAILURE;
		world >> ec;
		BOOST_TEST(ec == error_long_ns::SUCCESS);
	};
	break; default: {}
	}

	switch(world.rank()) {
	case 0: {
		error_class_ns::code const ec = error_class_ns::code::SUCCESS;
		world[1] << ec;
	};
	break; case 1: {
		error_class_ns::code ec = error_class_ns::code::FAILURE;
		world >> ec;
		BOOST_TEST(ec == error_class_ns::code::SUCCESS);
	};
	break; default: {}
	}

	return boost::report_errors();
} catch(...) {
	return 1;
}
