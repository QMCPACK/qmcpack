// Copyright 2022-2025 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/environment.hpp>

#include <mpi3/detail/datatype.hpp>

#include <boost/core/lightweight_test.hpp>
#include <complex>
#include <string>

#include <mpi3/detail/mpi_impl.h>

namespace mpi3 = boost::mpi3;

auto main(int argc, char** argv) -> int try {
	mpi3::environment env(argc, argv);

	auto world = env.world();

	using mpi3::detail::is_basic_v;

	static_assert(is_basic_v<int>);
	static_assert(is_basic_v<double>);
	static_assert(is_basic_v<mpi3::detail::float_int>);
	static_assert(is_basic_v<std::complex<double>>);

	static_assert(not is_basic_v<std::string>);

	BOOST_TEST( mpi3::detail::basic_datatype<double>{} == MPI_DOUBLE );

	return boost::report_errors();
} catch(...) {
	return 1;
}
