// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/adaptors/lapack/hipsolver.hpp>  // first: defines MULTI_USE_HIP for the other adaptors (the CI also passes it globally)

#include <boost/multi/adaptors/lapack/filling.hpp>
#include <boost/multi/adaptors/lapack/potrf.hpp>
#include <boost/multi/adaptors/lapack/syev.hpp>

#include <boost/multi/adaptors/thrust.hpp>
#include <boost/multi/array.hpp>

#include <thrust/complex.h>
#include <thrust/system/hip/memory.h>

#include <boost/core/lightweight_test.hpp>

#include <cmath>  // for abs

namespace multi = boost::multi;

namespace {
template<class T, multi::dimensionality_type D>
using device_array = multi::array<T, D, ::thrust::hip::allocator<T>>;
}  // end unnamed namespace

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	// double potrf on the device, classic Cholesky example: A = Uᵀ.U
	{
		multi::array<double, 2> const A_host = {
			{ 4.0, 12.0, -16.0},
			{ 0.0, 37.0, -43.0},  // the unused (lower) triangle can hold anything
			{ 0.0,  0.0,  98.0},
		};
		device_array<double, 2> A = A_host;

		auto&& res = multi::lapack::potrf(multi::lapack::filling::upper, A);
		BOOST_TEST( res.size() == 3 );  // fully positive definite

		multi::array<double, 2> const A_back = A;  // expected U = {{2, 6, -8}, {0, 1, 5}, {0, 0, 3}}
		BOOST_TEST( std::abs(A_back[0][0] - 2.0) < 1.0e-8 );
		BOOST_TEST( std::abs(A_back[0][1] - 6.0) < 1.0e-8 );
		BOOST_TEST( std::abs(A_back[0][2] + 8.0) < 1.0e-8 );
		BOOST_TEST( std::abs(A_back[1][1] - 1.0) < 1.0e-8 );
		BOOST_TEST( std::abs(A_back[1][2] - 5.0) < 1.0e-8 );
		BOOST_TEST( std::abs(A_back[2][2] - 3.0) < 1.0e-8 );
	}

	// complex (Hermitian) potrf on the device: same example, as complex
	// thrust::complex, not std::complex: the latter currently breaks the host <-> device array copy
	{
		using complex = ::thrust::complex<double>;

		multi::array<complex, 2> const A_host = {
			{ {4.0, 0.0}, {12.0, 0.0}, {-16.0, 0.0}},
			{ {0.0, 0.0}, {37.0, 0.0}, {-43.0, 0.0}},
			{ {0.0, 0.0}, { 0.0, 0.0}, { 98.0, 0.0}},
		};
		device_array<complex, 2> A = A_host;

		auto&& res = multi::lapack::potrf(multi::lapack::filling::upper, A);
		BOOST_TEST( res.size() == 3 );

		multi::array<complex, 2> const A_back = A;
		BOOST_TEST( abs(A_back[0][1] - complex{6.0, 0.0}) < 1.0e-8 );
		BOOST_TEST( abs(A_back[2][2] - complex{3.0, 0.0}) < 1.0e-8 );
	}

	// not-positive-definite potrf reports the order of the offending minor, like the CPU version
	{
		multi::array<double, 2> const A_host = {
			{ 1.0,  0.0, 0.0},
			{ 0.0, -1.0, 0.0},
			{ 0.0,  0.0, 1.0},
		};
		device_array<double, 2> A = A_host;

		auto&& res = multi::lapack::potrf(multi::lapack::filling::upper, A);
		BOOST_TEST( res.size() == 1 );  // only the leading 1x1 minor is positive definite
	}

	// double syev on the device, same data and expected values as the CPU test (test/syev.cpp)
	{
		multi::array<double, 2> const A_host = {
			{167.413, 126.804, 125.114},
			{  0.0  , 167.381, 126.746},
			{  0.0  ,   0.0  , 167.231},
		};
		device_array<double, 2> A = A_host;
		device_array<double, 1> W(A.size());  // NOLINT(readability-identifier-length) conventional lapack name

		auto&& res = multi::lapack::syev(multi::blas::filling::upper, A, W);
		BOOST_TEST( res.size() == A.size() );  // the returned view spans the block that converged

		multi::array<double, 2> const A_back = A;
		multi::array<double, 1> const W_back = W;

		BOOST_TEST( std::abs(W_back[1] - 42.2081) < 1.0e-3 );
		BOOST_TEST( std::abs(std::abs(A_back[2][1]) - 0.579092) < 1.0e-4 );  // eigenvector sign is not fixed
	}

	return boost::report_errors();
}
