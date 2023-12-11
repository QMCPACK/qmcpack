// Copyright 2020-2023 Alfredo A. Correa

#include "../main.hpp"
#include "../process.hpp"

#include <boost/serialization/vector.hpp>

struct long_long {
	long long  value;  // NOLINT(google-runtime-int,misc-non-private-member-variables-in-classes) testing type
	long_long& operator=(long long v) {  // NOLINT(google-runtime-int) testing type
		value = v;
		return *this;
	}
};

template<class Archive>
void serialize(Archive& ar, long_long& l, unsigned /*version*/= 0) {  // cppcheck-suppress unusedFunction ; false positive in cppcheck 2.11
	ar& l.value;
}

namespace mpi3 = boost::mpi3;

int mpi3::main(int /*argc*/, char** /*argv*/, mpi3::communicator world) try {

	assert( world.size() > 1 );
	long long size = 10000000;  // NOLINT(google-runtime-int) testing type
	switch(world.rank()) {
	case 0: {
		std::vector<long_long> v(static_cast<std::size_t>(size));
		std::iota(v.begin(), v.end(), 0.);
		//  assert(std::accumulate(v.begin(), v.end(), 0.) == size*(size-1)/2 );
		world[1] << v;
	} break;
	case 1: {
		std::vector<long_long> w;
		world[0] >> w;
		assert( w.size() == static_cast<std::size_t>(size) );
		assert( w[45].value == 45 );
		//  assert(std::accumulate(w.begin(), w.end(), 0.) == size*(size-1)/2 );
	} break;
	}
	return 0;
} catch(...) {return 1;}
