// Copyright 2018-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

#include <boost/core/lightweight_test.hpp>
#include <iterator>
#include <list>
#include <sstream>
#include <string>
#include <vector>

namespace mpi3 = boost::mpi3;

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	BOOST_TEST(world.size() > 1);

	switch(world.rank()) {
	case 0: {
		std::list<int> const b = {3, 4, 5};
		world.send(cbegin(b), cend(b), 1);
		break;
	};
	case 1: {
		std::vector<int> b2(3);
		world.receive(begin(b2), end(b2));
		BOOST_TEST(b2[1] == 4.0);
		break;
	};
	default: {
	}
	}
	// switch(world.rank()){
	//  case 0: {
	//      std::vector<std::string> b = {"hola", "blah", "chau"};
	//      world.send(cbegin(b), cend(b), 1);
	//  }; break;
	//  case 1: {
	//      std::list<std::string> b2(3);
	//      world.receive(begin(b2), end(b2));  // TODO(correaa) invesigate why it doesn't work
	//      BOOST_TEST( *begin(b2) == "hola" and *rbegin(b2) == "chau" );
	//  }; break;
	// }
	switch(world.rank()) {
	case 0: {
		std::istringstream iss{"1 2 3"};
		world.send(std::istream_iterator<int>{iss}, {}, 1);
		break;
	};
	case 1: {
		std::vector<int> out(3);
		world.receive(begin(out), end(out));
		BOOST_TEST((out == std::vector<int>{1, 2, 3}));
		break;
	}
	default: {
	}
	}

	return boost::report_errors();
} catch(...) {
	return 1;
}
