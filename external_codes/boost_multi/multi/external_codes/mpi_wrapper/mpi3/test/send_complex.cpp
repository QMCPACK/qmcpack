// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2023 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/main.hpp>

#include <complex>
#include <list>
#include <string>

namespace mpi3 = boost::mpi3;

template<class T>
struct my_complex {
	T _re;  // NOLINT(misc-non-private-member-variables-in-classes) aggregate
	T _im;  // NOLINT(misc-non-private-member-variables-in-classes) aggregate
	T real() const {return _re;}
	T imag() const {return _im;}
	bool operator==(my_complex const& other) const {return _re == other._re and _im == other._im;}
	bool operator!=(my_complex const& other) const {return _re != other._re or  _im != other._im;}

};

auto mpi3::main(int /*argc*/, char** /*argv*/, mpi3::communicator world) -> int try {
	assert(world.size() > 1);

	using complex = my_complex<double>;

	switch(world.rank()) {
	case 0: {
		complex const c{1.0, 2.0};
		world.send_n(&c, 1, 1);
		break;
	};
	case 1: {
		complex c;
		world.receive_n(&c, 1);
		assert((c == complex{1.0, 2.0}));
		break;
	};
	}
	return 0;
} catch(...) {
	return 1;
}
