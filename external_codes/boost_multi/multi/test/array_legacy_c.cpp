// Copyright 2019-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <array>
#include <complex>
#include <iostream>

// Suppress warnings from boost.test
#if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wold-style-cast"
#  pragma clang diagnostic ignored "-Wundef"
#  pragma clang diagnostic ignored "-Wconversion"
#  pragma clang diagnostic ignored "-Wsign-conversion"
#  pragma clang diagnostic ignored "-Wfloat-equal"
#elif defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wold-style-cast"
#  pragma GCC diagnostic ignored "-Wundef"
#  pragma GCC diagnostic ignored "-Wconversion"
#  pragma GCC diagnostic ignored "-Wsign-conversion"
#  pragma GCC diagnostic ignored "-Wfloat-equal"
#endif

#ifndef BOOST_TEST_MODULE
#  define BOOST_TEST_MAIN
#endif

#include <boost/test/unit_test.hpp>

#if defined(__clang__)
#  pragma clang diagnostic pop
#elif defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif

namespace multi = boost::multi;

namespace fake {

// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) testing a legacy interface
using fftw_complex = double[2];

void fftw_plan_dft(
	int rank, int const* n,
	fftw_complex* in, fftw_complex* out, int sign, unsigned flags
);

void fftw_plan_dft(
	int rank, int const* n,
	fftw_complex* in, fftw_complex* out, int sign, unsigned flags
) {
	(void)rank, (void)n, (void)in, (void)out, (void)sign, (void)flags;
}

}  // end namespace fake

BOOST_AUTO_TEST_CASE(array_legacy_c) {
	using complex                     = std::complex<double>;
	multi::array<complex, 2> const in = {
		{{150.0, 0.0}, {16.0, 0.0}, {17.0, 0.0}, {18.0, 0.0}, {19.0, 0.0}},
		{  {5.0, 0.0},  {5.0, 0.0},  {5.0, 0.0},  {5.0, 0.0},  {5.0, 0.0}},
		{{100.0, 0.0}, {11.0, 0.0}, {12.0, 0.0}, {13.0, 0.0}, {14.0, 0.0}},
		{ {50.0, 0.0},  {6.0, 0.0},  {7.0, 0.0},  {8.0, 0.0},  {9.0, 0.0}},
	};

	multi::array<std::complex<double>, 2> out(extensions(in));

	BOOST_REQUIRE( dimensionality(out) == dimensionality(in) );
	BOOST_REQUIRE( sizes(out) == sizes(in) );

	static_assert(sizeof(complex) == sizeof(fake::fftw_complex), "!");
	fake::fftw_plan_dft(
		decltype(in)::dimensionality,
		std::apply([](auto... sizes) { return std::array<int, 2>{{static_cast<int>(sizes)...}}; }, in.sizes()).data(),
		// NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-pro-type-const-cast) testing legacy code
		reinterpret_cast<fake::fftw_complex*>(const_cast<complex*>(in.data_elements())),  // NOSONAR
		// NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast): testing legacy code
		reinterpret_cast<fake::fftw_complex*>(out.data_elements()),
		1, 0
	);

	{
		multi::array<double, 2> d2D = {
			{150.0, 16.0, 17.0, 18.0, 19.0},
			{ 30.0,  1.0,  2.0,  3.0,  4.0},
			{100.0, 11.0, 12.0, 13.0, 14.0},
			{ 50.0,  6.0,  7.0,  8.0,  9.0},
		};

		//  #if __has_cpp_attribute(no_unique_address) >=201803L and not defined(__NVCC__) and not defined(__PGI)
		//      BOOST_REQUIRE( sizeof(d2D)==sizeof(double*)+7*sizeof(std::size_t) );
		//  #endif
		BOOST_REQUIRE( d2D.is_compact() );
		BOOST_REQUIRE( rotated(d2D).is_compact() );
		BOOST_REQUIRE( d2D[3].is_compact() );
		BOOST_REQUIRE( ! rotated(d2D)[2].is_compact() );
	}
	{
		multi::array<complex, 2> d2D({5, 3});
		BOOST_REQUIRE( d2D.is_compact() );
		BOOST_REQUIRE( rotated(d2D).is_compact() );
		BOOST_REQUIRE( d2D[3].is_compact() );
		BOOST_REQUIRE( ! rotated(d2D)[2].is_compact() );
	}
}

#ifndef _MSC_VER  // TODO(correaa) not supported by MSVC 14.3 in c++17 mode
constexpr auto f2(multi::array_ref<double, 1>&& array) -> double& { return std::move(array)[2]; }

BOOST_AUTO_TEST_CASE(array_legacy_c_2) {
	double arr[5] = {150.0, 16.0, 17.0, 18.0, 19.0};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
	BOOST_REQUIRE( &f2(arr) == &arr[2] );
}
#endif
