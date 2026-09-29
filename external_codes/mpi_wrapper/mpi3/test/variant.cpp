// Copyright 2023-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/process.hpp>

#include <boost/core/lightweight_test.hpp>
#include <cstddef>
#include <cstdint>
#include <tuple>
#include <variant>

namespace mpi3 = boost::mpi3;

enum class color : std::uint8_t {
	red,
	blue
};

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	BOOST_TEST(world.size() > 1);

	using std::variant;

	switch(world.rank()) {
	case 0: {
		variant<int, double> const v{3.14};
		world[1] << v;
	};
		break;
	case 1: {
		variant<int, double> v;
		world[0] >> v;
		BOOST_TEST((v == variant<int, double>{3.14}));
	};
	default: {
	}
	}

	switch(world.rank()) {
	case 0: {
		variant<color, double> const v{color::blue};
		world[1] << v;
	};
		break;
	case 1: {
		variant<color, double> v;
		world[0] >> v;
		BOOST_TEST((v == variant<color, double>{color::blue}));
	};
	default: {
	}
	}

	switch(world.rank()) {
	case 0: {
		variant<int, double> const v{3.14};
		world[1] << v;
	};
		break;
	case 1: {
		variant<int, double> v;
		world[0] >> v;
		BOOST_TEST((v == variant<int, double>{3.14}));
	};
	default: {
	}
	}

	switch(world.rank()) {
	case 0: {
		variant<int, double> const v{42};
		world.send_n(&v, 1, 1);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) testing not recommended brute force method
	};
		break;
	case 1: {
		variant<int, double> v;
		world.receive_n(&v, 1, 0);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) testing not recommended brute force method
		BOOST_TEST(v.index() == 0);
		BOOST_TEST((v == variant<int, double>{42}));
	};
	default: {
	}
	}

	switch(world.rank()) {
	case 0: {
		variant<int, std::tuple<int, int>> const v{std::make_tuple(5, 42)};
		world.send_n(reinterpret_cast<char const*>(&v), sizeof(v), 1);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) testing not recommended brute force method
	};
		break;
	case 1: {
		variant<int, std::tuple<int, int>> v;
		world.receive_n(reinterpret_cast<char*>(&v), sizeof(v), 0);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) testing not recommended brute force method
		BOOST_TEST((v == variant<int, std::tuple<int, int>>{std::make_tuple(5, 42)}));
	};
	default: {
	}
	}

	switch(world.rank()) {
	case 0: {
		variant<int, std::tuple<int, int>> const v{std::make_tuple(5, 42)};
		world.send_n(reinterpret_cast<char const*>(&v), sizeof(v), 1);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) testing not recommended brute force method
	};
		break;
	case 1: {
		variant<int, std::tuple<int, int>> v;
		world.receive_n(reinterpret_cast<char*>(&v), sizeof(v), 0);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) testing not recommended brute force method
		BOOST_TEST((v == variant<int, std::tuple<int, int>>{std::make_tuple(5, 42)}));
	};
	default: {
	}
	}

	switch(world.rank()) {
	case 0: {
		variant<int, std::tuple<int, int>> v{std::make_tuple(5, 42)};
		world.send_n(reinterpret_cast<std::byte*>(&v), sizeof(v), 1);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) testing not recommended brute force method
	};
		break;
	case 1: {
		variant<int, std::tuple<int, int>> v;
		world.receive_n(reinterpret_cast<std::byte*>(&v), sizeof(v), 0);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast) testing not recommended brute force method
		BOOST_TEST((v == variant<int, std::tuple<int, int>>{std::make_tuple(5, 42)}));
	};
	default: {
	}
	}

	// switch(world.rank()) {
	// break; case 0: {
	//  error_ns::code const ec = error_ns::SUCCESS;
	//  world.send_n(&ec, 1, 1);
	// };
	// break; case 1: {
	//  error_ns::code ec = error_ns::SUCCESS;
	//  world.receive_n(&ec, 1, 0);
	//  BOOST_TEST(ec == error_ns::SUCCESS);
	// };
	// }

	// switch(world.rank()) {
	// break; case 0: {
	//  error_ns::code const ec = error_ns::SUCCESS;
	//  world.send_n(&ec, 1, 1);
	// };
	// break; case 1: {
	//  error_ns::code ec = error_ns::SUCCESS;
	//  world.receive_n(&ec, 1);
	//  BOOST_TEST(ec == error_ns::SUCCESS);
	// };
	// }

	// switch(world.rank()) {
	// break; case 0: {
	//  error_ns::code const ec = error_ns::SUCCESS;
	//  world[1] << ec;
	// };
	// break; case 1: {
	//  error_ns::code ec = error_ns::SUCCESS;
	//  world[0] >> ec;
	//  BOOST_TEST(ec == error_ns::SUCCESS);
	// };
	// }

	// switch(world.rank()) {
	// break; case 0: {
	//  error_ns::code const ec = error_ns::SUCCESS;
	//  world[1] << ec;
	// };
	// break; case 1: {
	//  error_ns::code ec = error_ns::SUCCESS;
	//  world >> ec;
	//  BOOST_TEST(ec == error_ns::SUCCESS);
	// };
	// }

	// switch(world.rank()) {
	// break; case 0: {
	//  error_long_ns::code const ec = error_long_ns::SUCCESS;
	//  world[1] << ec;
	// };
	// break; case 1: {
	//  error_long_ns::code ec  = error_long_ns::SUCCESS;
	//  world >> ec;
	//  BOOST_TEST(ec == error_long_ns::SUCCESS);
	// };
	// }

	// switch(world.rank()) {
	// break; case 0: {
	//  error_class_ns::code const ec = error_class_ns::code::SUCCESS;
	//  world[1] << ec;
	// };
	// break; case 1: {
	//  error_class_ns::code ec = error_class_ns::code::SUCCESS;
	//  world >> ec;
	//  BOOST_TEST(ec == error_class_ns::code::SUCCESS);
	// };
	// }

	return boost::report_errors();
} catch(...) {
	return 1;
}
