// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS dot"
#include<boost/test/unit_test.hpp>

#include "../../blas/dot.hpp"

#include "../../../array.hpp"

#include<complex>
#include<numeric>
#include<type_traits>

namespace multi = boost::multi;
namespace blas = multi::blas;

BOOST_AUTO_TEST_CASE(blas_dot_context) {
	multi::array<float, 1> const x = {1., 2., 3.};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<float, 1> const y = {1., 2., 3.};  // NOLINT(readability-identifier-length) BLAS naming
	blas::context ctxt;
	{
		auto  res = +blas::dot(&ctxt, x, y);  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_REQUIRE( res == std::inner_product(begin(x), end(x), begin(y), 0.F) );
	}
	{
		float res = +blas::dot(&ctxt, x, y);  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_TEST_REQUIRE( res == std::inner_product(begin(x), end(x), begin(y), 0.F) );
	}
}

BOOST_AUTO_TEST_CASE(blas_dot_no_context) {
	multi::array<float, 1> const x = {1., 2., 3.};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<float, 1> const y = {1., 2., 3.};  // NOLINT(readability-identifier-length) BLAS naming
	auto res = +blas::dot(x, y);
	BOOST_REQUIRE( res == std::inner_product(begin(x), end(x), begin(y), 0.F) );
}


BOOST_AUTO_TEST_CASE(blas_dot_no_context_out_param) {
	multi::array<float, 1> const x = {1., 2., 3.};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<float, 1> const y = {1., 2., 3.};  // NOLINT(readability-identifier-length) BLAS naming
	float res = NAN;
	blas::dot(x, y, res);
	BOOST_REQUIRE( res == std::inner_product(begin(x), end(x), begin(y), 0.F) );
}

BOOST_AUTO_TEST_CASE(blas_dot_no_context_out_param_complex) {  // if you get a segfaut here, your system may require -DRETURN_BY_STACK
	using complex = std::complex<double>;
	multi::array<complex, 1> const x = {1., 2., 3.};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 1> const y = {1., 2., 3.};  // NOLINT(readability-identifier-length) BLAS naming
	complex res;
	blas::dot(x, y, res);
	BOOST_REQUIRE_EQUAL( real(res) , real(std::inner_product(begin(x), end(x), begin(y), complex{0.}, std::plus<>{}, [](auto const& alpha, auto const& omega) {return alpha*std::conj(omega);})) );
	BOOST_REQUIRE_EQUAL( imag(res) , imag(std::inner_product(begin(x), end(x), begin(y), complex{0.}, std::plus<>{}, [](auto const& alpha, auto const& omega) {return alpha*std::conj(omega);})) );
}

BOOST_AUTO_TEST_CASE(blas_dot_no_context_out_param_complex_C) {
	using complex = std::complex<double>; complex const I{0, 1};  // NOLINT(readability-identifier-length) imag unit
	multi::array<complex, 1> const x = {1., 2.       , 3.};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<complex, 1> const y = {1., 2. + 2.*I, 3.};  // NOLINT(readability-identifier-length) BLAS naming
	complex res;
	blas::dot(blas::C(x), y, res);
	BOOST_REQUIRE( res == std::inner_product(begin(x), end(x), begin(y), complex{0.}, std::plus<>{}, [](auto const& alpha, auto const& omega) {return conj(alpha)*omega;}) );
}

#if defined(CUDA_FOUND) and CUDA_FOUND
#include<thrust/complex.h>
BOOST_AUTO_TEST_CASE(blas_dot_no_context_out_param_complex_C_thrust) {
	using complex = thrust::complex<double>; complex const I{0., 1.};
	multi::array<complex, 1> const A = {1., 2.       , 3.};
	multi::array<complex, 1> const B = {1., 2. + 2.*I, 3.};
	complex C;
	blas::dot(blas::C(A), B, C);
	BOOST_REQUIRE( C == std::inner_product(begin(A), end(A), begin(B), complex{0.}, std::plus<>{}, [](auto& a, auto& b){return conj(a)*b;}) );
}
#endif

BOOST_AUTO_TEST_CASE(multi_blas_dot_strided) {
	multi::array<double, 2> const CA = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	double res = std::numeric_limits<double>::quiet_NaN();
	blas::dot_n(begin(CA[1]), size(CA[1]), begin(CA[2]), &res);
	BOOST_REQUIRE( res == std::inner_product(begin(CA[1]), begin(CA[2]), end(CA[1]), 0.) );

	double res2 = blas::dot(CA[1], CA[2]);
	BOOST_REQUIRE( res == res2 );
}

BOOST_AUTO_TEST_CASE(multi_blas_dot_strided_context) {
	multi::array<double, 2> const CA = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	double res = std::numeric_limits<double>::quiet_NaN();
	blas::context ctxt;
	blas::dot_n(&ctxt, begin(CA[1]), size(CA[1]), begin(CA[2]), &res);
	BOOST_REQUIRE( res == std::inner_product(begin(CA[1]), begin(CA[2]), end(CA[1]), 0.) );

	double res2 = blas::dot(CA[1], CA[2]);
	BOOST_REQUIRE( res == res2 );
}

BOOST_AUTO_TEST_CASE(multi_blas_dot_1d_real) {
	multi::array<float, 1> x = {1., 2., 3.};  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<float, 1> y = {1., 2., 3.};  // NOLINT(readability-identifier-length) BLAS naming

	using blas::dot;
	BOOST_REQUIRE( 14. == dot(x, y) );
	BOOST_REQUIRE( dot(x, y) == 14.F );
}

BOOST_AUTO_TEST_CASE(multi_blas_dot_impl_real) {
	multi::array<double, 2> const cA = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	{
		double res = blas::dot(cA[1], cA[2]);
		BOOST_REQUIRE( res == std::inner_product(begin(cA[1]), begin(cA[2]), end(cA[1]), 0.) );
	}
	{
		double res = NAN;
		blas::dot(cA[1], cA[2], res);
		BOOST_REQUIRE( res == std::inner_product(begin(cA[1]), begin(cA[2]), end(cA[1]), 0.) );
	}
	{
		double res = NAN;
		auto res2 = blas::dot(cA[1], cA[2], res);
		BOOST_REQUIRE( res == res2 );
	}
	 {
		double res = blas::dot(cA[1], cA[2]);
		BOOST_REQUIRE( res == std::inner_product(begin(cA[1]), begin(cA[2]), end(cA[1]), 0.) );
		BOOST_REQUIRE( blas::dot(cA[1], cA[2]) == blas::dot(cA[2], cA[1]) );
	}
}

BOOST_AUTO_TEST_CASE(inq_case) {
	multi::array<double, 1> x(multi::extensions_t<1>{multi::iextension{10}}, +1.0);  // NOLINT(readability-identifier-length) BLAS naming
	multi::array<double, 1> y(multi::extensions_t<1>{multi::iextension{10}}, -1.0);  // NOLINT(readability-identifier-length) BLAS naming

	using blas::dot;
	using blas::hermitized;
	using blas::conj;

	auto res = dot(x, y);
	auto res2 = dot(hermitized(x), y);

	BOOST_REQUIRE(res == res2);

	auto res3 = dot(blas::conj(x), y);  // conjugation doesn't do anything for real array
	BOOST_REQUIRE(res3 == res);

	auto d_arr = dot(blas::C(x), y);
	BOOST_REQUIRE(d_arr == res);

	static_assert( not std::is_same<decltype(d_arr), double>{}, "!" );

	using blas::C;
	double d_doub = dot(C(x), y);

	BOOST_REQUIRE( d_doub == d_arr );
}

BOOST_AUTO_TEST_CASE(multi_blas_dot_impl_complex) {
	namespace blas = multi::blas;

	using complex = std::complex<double>; complex const I{0, 1};  // NOLINT(readability-identifier-length) imag unit
	multi::array<complex, 2> const A = {  // NOLINT(readability-identifier-length) BLAS naming
		{1. +    I,  2. + 3.*I,  3.+2.*I,  4.-9.*I},
		{5. + 2.*I,  6. + 6.*I,  7.+2.*I,  8.-3.*I},
		{9. + 1.*I, 10. + 9.*I, 11.+1.*I, 12.+2.*I}
	};
	{
		complex c; blas::dot(A[1], A[2], c);  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_TEST_REQUIRE( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0}) );
	}
	{
		complex c = blas::dot(A[1], A[2]);  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_TEST_REQUIRE( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0}) );
	}
	{
		complex c = blas::dot(A[1], blas::C(A[2]));  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_TEST_REQUIRE( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0}, std::plus<>{}, [](auto alpha, auto omega) {return alpha*conj(omega);}) );
	}
	{
		complex c = blas::dot(blas::C(A[1]), A[2]);  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_TEST_REQUIRE( c == inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{}, std::plus<>{}, [](auto alpha, auto omega) {return conj(alpha)*omega;}) );
	}
	{
		complex c = blas::dot(blas::conj(A[1]), A[2]);  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_TEST_REQUIRE( c == inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{}, std::plus<>{}, [](auto alpha, auto omega) {return conj(alpha)*omega;}) );
	}
	{
		complex c = blas::dot(blas::C(A[1]), A[2]);  // NOLINT(readability-identifier-length) BLAS naming
		BOOST_TEST_REQUIRE( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0}, std::plus<>{}, [](auto alpha, auto omega) {return conj(alpha)*omega;}) );
	}
}

#include "config.hpp" // cuda found
#if defined(CUDA_FOUND) and CUDA_FOUND
//#include<thrust/complex.h>

BOOST_AUTO_TEST_CASE(multi_blas_dot_impl_complex_thrust) {
	namespace blas = multi::blas;

	using complex = std::complex<double>; complex const I{0, 1};
	multi::array<complex, 2> const A = {
		{1. +    I,  2. + 3.*I,  3.+2.*I,  4.-9.*I},
		{5. + 2.*I,  6. + 6.*I,  7.+2.*I,  8.-3.*I},
		{9. + 1.*I, 10. + 9.*I, 11.+1.*I, 12.+2.*I}
	};
	{
		complex c;
		blas::core::dotu(size(A[1]), A[1].base(), A[1].stride(), A[2].base(), A[2].stride(), &c);
		auto inner = std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.});
		BOOST_REQUIRE( c.real() == inner.real() );
		BOOST_REQUIRE( c.imag() == inner.imag() );
	}
	{
		complex c;
		blas::context::dotu(size(A[1]), A[1].base(), A[1].stride(), A[2].base(), A[2].stride(), &c);
		auto inner = std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.});
		BOOST_REQUIRE( c.real() == inner.real() );
		BOOST_REQUIRE( c.imag() == inner.imag() );
	}
	{
		complex c;
		blas::dot_n(begin(A[1]), size(A[1]), begin(A[2]), &c);
		auto inner = std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.});
		BOOST_REQUIRE( c == inner );
	}
	{
		complex c;
		blas::dot(A[1], A[2], c);
		auto inner = std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.});
		BOOST_REQUIRE( c == inner );
	}
	{
		complex c = blas::dot(A[1], A[2]);
		auto inner = std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.});
		BOOST_REQUIRE( c == inner );
	}
	{
		auto inner = std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0.});
		BOOST_REQUIRE( +blas::dot(A[1], A[2]) == inner );
	}
	{
		complex c; blas::dot(A[1], A[2], c);
		BOOST_REQUIRE( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0}) );
	}
	{
		complex c = blas::dot(A[1], A[2]);
		BOOST_REQUIRE( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0}) );
	}
	 {
		complex c = blas::dot(A[1], blas::C(A[2]));
		BOOST_REQUIRE( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0}, std::plus<>{}, [](auto a, auto b){return a*conj(b);}) );
	}
}
#endif
