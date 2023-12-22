// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

// #define BOOST_TEST_MODULE "C++ Unit Tests for Multi legacy adaptor example"  // NOLINT(cppcoreguidelines-macro-usage) title
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

#include<complex>
#include<iostream>

namespace multi = boost::multi;

namespace fake {

// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) testing a legacy interface
using fftw_complex = double[2];

void fftw_plan_dft(
	int rank, const int* n,
	fftw_complex* in, fftw_complex* out, int sign, unsigned flags
);

void fftw_plan_dft(
	int rank, const int* n,
	fftw_complex* in, fftw_complex* out, int sign, unsigned flags
) {
	(void)rank, (void)n, (void)in, (void)out, (void)sign, (void)flags;
}

}  // end namespace fake

BOOST_AUTO_TEST_CASE(array_legacy_c) {
	using complex = std::complex<double>;
	multi::array<complex, 2> const in = {
		{{150.0, 0.0}, {16.0, 0.0}, {17.0, 0.0}, {18.0, 0.0}, {19.0, 0.0}},
		{{  5.0, 0.0}, { 5.0, 0.0}, { 5.0, 0.0}, { 5.0, 0.0}, { 5.0, 0.0}},
		{{100.0, 0.0}, {11.0, 0.0}, {12.0, 0.0}, {13.0, 0.0}, {14.0, 0.0}},
		{{ 50.0, 0.0}, { 6.0, 0.0}, { 7.0, 0.0}, { 8.0, 0.0}, { 9.0, 0.0}}
	};

	multi::array<std::complex<double>, 2> out(extensions(in));

	BOOST_REQUIRE( dimensionality(out) == dimensionality(in) );
	BOOST_REQUIRE( sizes(out) == sizes(in) );

	static_assert( sizeof(complex) == sizeof(fake::fftw_complex), "!" );
	fake::fftw_plan_dft(
		decltype(in)::dimensionality,
		std::apply([](auto... sizes) {return std::array<int, 2>{{static_cast<int>(sizes)...}};}, in.sizes()).data(),
		// NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-pro-type-const-cast) testing legacy code
		reinterpret_cast<fake::fftw_complex*>(const_cast<complex*>(in .data_elements())),
		// NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast): testing legacy code
		reinterpret_cast<fake::fftw_complex*>(                     out.data_elements() ),
		1, 0
	);

	{
		multi::array<double, 2> d2D = {
			{150.0, 16.0, 17.0, 18.0, 19.0},
			{ 30.0,  1.0,  2.0,  3.0,  4.0},
			{100.0, 11.0, 12.0, 13.0, 14.0},
			{ 50.0,  6.0,  7.0,  8.0,  9.0}
		};

//  #if __has_cpp_attribute(no_unique_address) >=201803L and not defined(__NVCC__) and not defined(__PGI)
//      BOOST_REQUIRE( sizeof(d2D)==sizeof(double*)+7*sizeof(std::size_t) );
//  #endif
		BOOST_REQUIRE( d2D.is_compact() );
		BOOST_REQUIRE( rotated(d2D).is_compact() );
		BOOST_REQUIRE( d2D[3].is_compact() );
		BOOST_REQUIRE( not rotated(d2D)[2].is_compact() );
	}
	{
		multi::array<complex, 2> d2D({5, 3});
		BOOST_REQUIRE( d2D.is_compact() );
		BOOST_REQUIRE( rotated(d2D).is_compact() );
		BOOST_REQUIRE( d2D[3].is_compact() );
		BOOST_REQUIRE( not rotated(d2D)[2].is_compact() );
	}
}

inline constexpr auto f2(multi::array_ref<double, 1>&& array) -> double& {return array[2];}

BOOST_AUTO_TEST_CASE(array_legacy_c_2) {
	double arr[5] = {150.0, 16.0, 17.0, 18.0, 19.0};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
	BOOST_REQUIRE( &f2(arr) == &arr[2] );
}
