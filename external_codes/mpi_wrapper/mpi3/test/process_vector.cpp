// Copyright 2020-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/process.hpp>

#include <boost/serialization/vector.hpp>  // NOLINT(misc-include-cleaner) used indrectly

#include <boost/core/lightweight_test.hpp>
#include <cstddef>
#include <numeric>
#include <vector>

struct long_long {
	long long  value;                    // NOLINT(google-runtime-int,misc-non-private-member-variables-in-classes) testing type
	long_long& operator=(long long v) {  // NOLINT(google-runtime-int) testing type
		value = v;
		return *this;
	}
};

template<class Archive>
static void serialize(Archive& ar, long_long& l, unsigned /*version*/ = 0) {  // NOLINT(misc-use-anonymous-namespace) used by boost.serialization  // cppcheck-suppress [unusedFunction,unmatchedSuppression] ; false positive in cppcheck 2.11
	ar & l.value;
}

namespace mpi3 = boost::mpi3;

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	BOOST_TEST(world.size() > 1);
	constexpr long long size = 10000;  // 000;  // NOLINT(google-runtime-int) testing type
	switch(world.rank()) {
	case 0: {
		std::vector<long_long> v(static_cast<std::size_t>(size));
		std::iota(v.begin(), v.end(), 0LL);  // NOLINT(boost-use-ranges) for C++20, use std::ranges::iota
		//  BOOST_TEST(std::accumulate(v.begin(), v.end(), 0.) == size*(size-1)/2 );
		world[1] << v;
	} break;
	case 1: {
		std::vector<long_long> w;
		world[0] >> w;
		BOOST_TEST(w.size() == static_cast<std::size_t>(size));
		BOOST_TEST(w[45].value == 45);
		//  BOOST_TEST(std::accumulate(w.begin(), w.end(), 0.) == size*(size-1)/2 );
	} break;
	default: {
	}
	}

	return boost::report_errors();
} catch(...) {
	return 1;
}
