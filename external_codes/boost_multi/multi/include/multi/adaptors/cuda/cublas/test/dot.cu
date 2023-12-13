// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2023 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi CUBLAS dot"
#include<boost/test/unit_test.hpp>

#include <multi/adaptors/cuda/cublas.hpp>

#include <multi/adaptors/blas/dot.hpp>
#include <multi/adaptors/blas/axpy.hpp>
#include <multi/adaptors/blas/gemm.hpp>
#include <multi/adaptors/blas/nrm2.hpp>
#include <multi/adaptors/blas/scal.hpp>

#include <multi/adaptors/thrust.hpp>

#include<thrust/complex.h>

#include<numeric>

namespace multi = boost::multi;

// BOOST_AUTO_TEST_CASE(cublas_dot_out_param_complex_C) {
//  namespace blas = multi::blas;
//  using complex = thrust::complex<double>;
//  complex const I{0.0, 1.0};

//  multi::thrust::cuda::array<complex, 1> const x = {1.0 + 0.0*I, 2.0 + 0.0*I, 3.0 + 0.0*I};  // NOLINT(readability-identifier-length) BLAS naming
//  multi::thrust::cuda::array<complex, 1> const y = {1.0 + 0.0*I, 2.0 + 2.0*I, 3.0 + 0.0*I};  // NOLINT(readability-identifier-length) BLAS naming

//  complex res{0.0, 0.0};
//  blas::dot(blas::C(x), y, res);
// //  BOOST_REQUIRE( res == std::inner_product(begin(x), end(x), begin(y), complex{0.0, 0.0}, std::plus<>{}, [](auto const& alpha, auto const& omega) {return conj(alpha)*omega;}) );
// }

BOOST_AUTO_TEST_CASE(cublas_dot_out_array0D_complex_C) {
	namespace blas = multi::blas;
	using complex = thrust::complex<double>;
	complex const I{0.0, 1.0};

	multi::thrust::cuda::array<complex, 1> const x = {1.0 + 0.0*I, 2.0 + 0.0*I, 3.0 + 0.0*I};  // NOLINT(readability-identifier-length) BLAS naming
	multi::thrust::cuda::array<complex, 1> const y = {1.0 + 0.0*I, 2.0 + 2.0*I, 3.0 + 0.0*I};  // NOLINT(readability-identifier-length) BLAS naming

	multi::thrust::cuda::array<complex, 0> res{complex{0.0, 0.0}};
	blas::dot(blas::C(x), y, res);

	{
		multi::array<complex, 0> res_copy{complex{0.0, 0.0}};
		res_copy = res;
		BOOST_REQUIRE(( *res_copy.base() == complex{14.0, 4.0} ));
	}
	{
		multi::array<complex, 0> res_copy{res};
		BOOST_REQUIRE(( *res_copy.base() == complex{14.0, 4.0} ));
	}
}

// BOOST_AUTO_TEST_CASE(blas_dot_functional_complex_C) {
//  namespace blas = multi::blas;
//  using complex = thrust::complex<double>;
//  complex const I{0.0, 1.0};

//  multi::thrust::cuda::array<complex, 1> const x = {1.0 + 0.0*I, 2.0 + 0.0*I, 3.0 + 0.0*I};  // NOLINT(readability-identifier-length) BLAS naming
//  multi::thrust::cuda::array<complex, 1> const y = {1.0 + 0.0*I, 2.0 + 2.0*I, 3.0 + 0.0*I};  // NOLINT(readability-identifier-length) BLAS naming

//  complex res = blas::dot(blas::C(x), y);
//  BOOST_REQUIRE( res == std::inner_product(begin(x), end(x), begin(y), complex{0.0, 0.0}, std::plus<>{}, [](auto const& alpha, auto const& omega) {return conj(alpha)*omega;}) );
// }

// BOOST_AUTO_TEST_CASE(blas_dot_functional_mutate_complex_C) {
//  namespace blas = multi::blas;
//  using complex = thrust::complex<double>;
//  complex const I{0.0, 1.0};

//  multi::thrust::cuda::array<complex, 1> const x = {1.0 + 0.0*I, 2.0 + 0.0*I, 3.0 + 0.0*I};  // NOLINT(readability-identifier-length) BLAS naming
//  multi::thrust::cuda::array<complex, 1> const y = {1.0 + 0.0*I, 2.0 + 2.0*I, 3.0 + 0.0*I};  // NOLINT(readability-identifier-length) BLAS naming

//  complex res;
//  res = blas::dot(blas::C(x), y);
//  BOOST_REQUIRE( res == std::inner_product(begin(x), end(x), begin(y), complex{0.0, 0.0}, std::plus<>{}, [](auto const& alpha, auto const& omega) {return conj(alpha)*omega;}) );
// }
