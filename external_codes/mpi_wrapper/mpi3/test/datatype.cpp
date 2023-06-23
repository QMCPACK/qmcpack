// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2022-2023 Alfredo A. Correa

#include <mpi3/communicator.hpp>
#include <mpi3/main.hpp>

#include <complex>
#include <list>
#include <string>

namespace mpi3 = boost::mpi3;

auto mpi3::main(int /*argc*/, char** /*argv*/, mpi3::communicator /*world*/) -> int try {
	using mpi3::detail::is_basic_v;

	static_assert( is_basic_v<int> );
	static_assert( is_basic_v<double> );
	static_assert( is_basic_v<mpi3::detail::float_int> );
	static_assert( is_basic_v<std::complex<double>> );

	static_assert( not is_basic_v<std::string> );

	assert( mpi3::detail::basic_datatype<double>{} == MPI_DOUBLE );

	return 0;
} catch(...) {
	return 1;
}
