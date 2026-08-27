// Copyright 2023-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/detail/datatype.hpp>
#include <mpi3/environment.hpp>
#include <mpi3/type.hpp>

#include <boost/core/lightweight_test.hpp>
#include <complex>
#include <iterator>
#include <vector>

namespace mpi3 = boost::mpi3;

struct spinor {
	std::complex<double> up;  // NOLINT(misc-non-private-member-variables-in-classes)
	std::complex<double> dn;  // NOLINT(misc-non-private-member-variables-in-classes)

	bool operator==(spinor const& other) const { return up == other.up && dn == other.dn; }
	bool operator!=(spinor const& other) const { return up != other.up || dn != other.dn; }
};

template<> struct mpi3::datatype<spinor> : mpi3::struct_<std::complex<double>, std::complex<double>> {};

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	using namespace std::complex_literals;  // i

	auto const I = std::complex<double>{0.0, 1.0};

	switch(world.rank()) {
	case 0: {
		std::vector<spinor> v(5);
		v[2] = spinor{3.14 + 6.28 * I, 4.0 + 5.0 * I};
		world.send_n(begin(v), 5, 1);
		break;
	};
	case 1: {
		std::vector<spinor> v(5);
		world.receive_n(begin(v), 5, 0);
		BOOST_TEST((v[2] == spinor{3.14 + 6.28 * I, 4.0 + 5.0 * I}));
		break;
	};
	default: ;
	}

	static_assert(boost::mpi3::has_datatype<spinor>{});

	return boost::report_errors();
} catch(...) {
	return 1;
}
