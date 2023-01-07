// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/main.hpp>

#include <complex>
#include <list>
#include <string>

namespace mpi3 = boost::mpi3;

auto mpi3::main(int /*argc*/, char** /*argv*/, mpi3::communicator world) -> int try {
	assert(world.size() > 1);

	switch(world.rank()) {
	case 0: {
		std::list<int> b = {3, 4, 5};
		world.send(cbegin(b), cend(b), 1);
		break;
	};
	case 1: {
		std::vector<int> b2(3);
		world.receive(begin(b2), end(b2));
		assert(b2[1] == 4.0);
		break;
	};
	}
	// switch(world.rank()){
	// 	case 0: {
	// 		std::vector<std::string> b = {"hola", "blah", "chau"};
	// 		world.send(cbegin(b), cend(b), 1);
	// 	}; break;
	// 	case 1: {
	// 		std::list<std::string> b2(3);
	// 		world.receive(begin(b2), end(b2));  // TODO(correaa) invesigate why it doesn't work
	// 		assert( *begin(b2) == "hola" and *rbegin(b2) == "chau" );
	// 	}; break;
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
		assert((out == std::vector<int>{1, 2, 3}));
		break;
	}
	}

	return 0;
} catch(...) {
	return 1;
}
