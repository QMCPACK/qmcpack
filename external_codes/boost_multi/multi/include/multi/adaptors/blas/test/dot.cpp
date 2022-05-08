// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo A. Correa 2019-2021

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS dot"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../blas/dot.hpp"

#include "../../../array.hpp"

#include<complex>
#include<numeric>
#include<type_traits>

namespace multi = boost::multi;
namespace blas = multi::blas;

BOOST_AUTO_TEST_CASE(blas_dot_context) {
	multi::array<float, 1> const A = {1., 2., 3.};
	multi::array<float, 1> const B = {1., 2., 3.};
	blas::context ctxt;
	auto C = +blas::dot(&ctxt, A, B);
	float c = +blas::dot(&ctxt, A, B);
	BOOST_TEST_REQUIRE( c == std::inner_product(begin(A), end(A), begin(B), 0.F) );
	BOOST_REQUIRE( C == std::inner_product(begin(A), end(A), begin(B), 0.F) );
}

BOOST_AUTO_TEST_CASE(blas_dot_no_context) {
	multi::array<float, 1> const A = {1., 2., 3.};
	multi::array<float, 1> const B = {1., 2., 3.};
	auto C = +blas::dot(A, B);
	BOOST_REQUIRE( C == std::inner_product(begin(A), end(A), begin(B), 0.F) );
}


BOOST_AUTO_TEST_CASE(blas_dot_no_context_out_param) {
	multi::array<float, 1> const A = {1., 2., 3.};
	multi::array<float, 1> const B = {1., 2., 3.};
	float C = NAN;
	blas::dot(A, B, C);
	BOOST_REQUIRE( C == std::inner_product(begin(A), end(A), begin(B), 0.F) );
}

BOOST_AUTO_TEST_CASE(blas_dot_no_context_out_param_complex) { // if you get a segfaut here, your system may require -DRETURN_BY_STACK
	using complex = std::complex<double>;
	multi::array<complex, 1> const A = {1., 2., 3.};
	multi::array<complex, 1> const B = {1., 2., 3.};
	complex C;
	blas::dot(A, B, C);
	BOOST_REQUIRE_EQUAL( real(C) , real(std::inner_product(begin(A), end(A), begin(B), complex{0.}, std::plus<>{}, [](auto const& a, auto const& b) {return a*std::conj(b);})) );
	BOOST_REQUIRE_EQUAL( imag(C) , imag(std::inner_product(begin(A), end(A), begin(B), complex{0.}, std::plus<>{}, [](auto const& a, auto const& b) {return a*std::conj(b);})) );
}

BOOST_AUTO_TEST_CASE(blas_dot_no_context_out_param_complex_C) {
	using complex = std::complex<double>; complex const I{0., 1.};
	multi::array<complex, 1> const A = {1., 2.       , 3.};
	multi::array<complex, 1> const B = {1., 2. + 2.*I, 3.};
	complex C;
	blas::dot(blas::C(A), B, C);
	BOOST_REQUIRE( C == std::inner_product(begin(A), end(A), begin(B), complex{0.}, std::plus<>{}, [](auto const& a, auto const& b){return conj(a)*b;}) );
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
	double d = std::numeric_limits<double>::quiet_NaN();
	blas::dot_n(begin(CA[1]), size(CA[1]), begin(CA[2]), &d);
	BOOST_REQUIRE( d == std::inner_product(begin(CA[1]), begin(CA[2]), end(CA[1]), 0.) );

	double d2 = blas::dot(CA[1], CA[2]);
	BOOST_REQUIRE( d == d2 );
}

BOOST_AUTO_TEST_CASE(multi_blas_dot_strided_context) {
	multi::array<double, 2> const CA = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	double d = std::numeric_limits<double>::quiet_NaN();
	blas::context ctxt;
	blas::dot_n(&ctxt, begin(CA[1]), size(CA[1]), begin(CA[2]), &d);
	BOOST_REQUIRE( d == std::inner_product(begin(CA[1]), begin(CA[2]), end(CA[1]), 0.) );

	double d2 = blas::dot(CA[1], CA[2]);
	BOOST_REQUIRE( d == d2 );
}

BOOST_AUTO_TEST_CASE(multi_blas_dot_1d_real) {
	multi::array<float, 1> V = {1., 2., 3.};
	multi::array<float, 1> W = {1., 2., 3.};

	using blas::dot;
	BOOST_REQUIRE( 14. == dot(V, W) );
	BOOST_REQUIRE( dot(V, W) == 14.F );
}

BOOST_AUTO_TEST_CASE(multi_blas_dot_impl_real) {
	multi::array<double, 2> const cA = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
	{
		double d = blas::dot(cA[1], cA[2]);
		BOOST_REQUIRE( d==std::inner_product(begin(cA[1]), begin(cA[2]), end(cA[1]), 0.) );
	}
	{
		double d = NAN;
		blas::dot(cA[1], cA[2], d);
		BOOST_REQUIRE( d==std::inner_product(begin(cA[1]), begin(cA[2]), end(cA[1]), 0.) );
	}
	{
		double d = NAN;
		auto d2 = blas::dot(cA[1], cA[2], d);
		BOOST_REQUIRE( d==d2 );
	}
	 {
		double d = blas::dot(cA[1], cA[2]);
		BOOST_REQUIRE( d == std::inner_product(begin(cA[1]), begin(cA[2]), end(cA[1]), 0.) );
		BOOST_REQUIRE( blas::dot(cA[1], cA[2]) == blas::dot(cA[2], cA[1]) );
	}
}

BOOST_AUTO_TEST_CASE(inq_case) {
	multi::array<double, 1> v1(multi::extensions_t<1>{multi::iextension{10}}, +1.0);
	multi::array<double, 1> v2(multi::extensions_t<1>{multi::iextension{10}}, -1.0);

	using blas::dot;
	using blas::hermitized;
	using blas::conj;

	auto a = dot(v1, v2);
	auto b = dot(hermitized(v1), v2);

	BOOST_REQUIRE(a == b);

	auto c = dot(blas::conj(v1), v2); // conjugation doesn't do anything for real array
	BOOST_REQUIRE(c == a);

	auto d_arr = dot(blas::C(v1), v2);
	BOOST_REQUIRE(d_arr == a);

	static_assert( not std::is_same<decltype(d_arr), double>{}, "!" );

	using blas::C;
	double d_doub = dot(C(v1), v2);

	BOOST_REQUIRE( d_doub == d_arr );
}

BOOST_AUTO_TEST_CASE(multi_blas_dot_impl_complex) {
	namespace blas = multi::blas;

	using complex = std::complex<double>; complex const I{0, 1};
	multi::array<complex, 2> const A = {
		{1. +    I,  2. + 3.*I,  3.+2.*I,  4.-9.*I},
		{5. + 2.*I,  6. + 6.*I,  7.+2.*I,  8.-3.*I},
		{9. + 1.*I, 10. + 9.*I, 11.+1.*I, 12.+2.*I}
	};
	{
		complex c; blas::dot(A[1], A[2], c);
		BOOST_TEST_REQUIRE( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0}) );
	}
	{
		complex c = blas::dot(A[1], A[2]);
		BOOST_TEST_REQUIRE( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0}) );
	}
	{
		complex c = blas::dot(A[1], blas::C(A[2]));
		BOOST_TEST_REQUIRE( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0}, std::plus<>{}, [](auto a, auto b){return a*conj(b);}) );
	}
	{
		complex c = blas::dot(blas::C(A[1]), A[2]);
		BOOST_TEST_REQUIRE( c == inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{}, std::plus<>{}, [](auto a, auto b){return conj(a)*b;}) );
	}
	{
		complex c = blas::dot(blas::conj(A[1]), A[2]);
		BOOST_TEST_REQUIRE( c == inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{}, std::plus<>{}, [](auto a, auto b){return conj(a)*b;}) );
	}
	 {
		complex c = blas::dot(blas::C(A[1]), A[2]);
		BOOST_TEST_REQUIRE( c == std::inner_product(begin(A[1]), end(A[1]), begin(A[2]), complex{0}, std::plus<>{}, [](auto a, auto b){return conj(a)*b;}) );
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
