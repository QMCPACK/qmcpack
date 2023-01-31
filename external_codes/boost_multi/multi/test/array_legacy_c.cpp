// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi legacy adaptor example"
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

} // end namespace fake

BOOST_AUTO_TEST_CASE(array_legacy_c) {
	using complex = std::complex<double>;
	multi::array<complex, 2> const in = {
		{150., 16., 17., 18., 19.},
		{  5.,  5.,  5.,  5.,  5.},
		{100., 11., 12., 13., 14.},
		{ 50.,  6.,  7.,  8.,  9.}
	};

	multi::array<std::complex<double>, 2> out(extensions(in));

	BOOST_REQUIRE( dimensionality(out) == dimensionality(in) );
	BOOST_REQUIRE( sizes(out) == sizes(in) );

	static_assert( sizeof(complex) == sizeof(fake::fftw_complex), "!" );
	fake::fftw_plan_dft(
		dimensionality(in),
		std::apply([](auto... sizes) {return std::array<int, 2>{{static_cast<int>(sizes)...}};}, in.sizes()).data(),
		reinterpret_cast<fake::fftw_complex*>(const_cast<complex*>(in .data_elements())), // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast, cppcoreguidelines-pro-type-const-cast): testing legacy code
		reinterpret_cast<fake::fftw_complex*>(                     out.data_elements() ), // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast): testing legacy code
		1, 0
	);

//struct basic : multi::layout_t<2>{
//	double* p = {};
//	basic() = default;
//};

//struct ref : basic{
//};


	{
		multi::array<double, 2> d2D = {
			{150, 16, 17, 18, 19},
			{ 30,  1,  2,  3,  4},
			{100, 11, 12, 13, 14},
			{ 50,  6,  7,  8,  9}
		};

//	#if __has_cpp_attribute(no_unique_address) >=201803L and not defined(__NVCC__) and not defined(__PGI)
//		BOOST_REQUIRE( sizeof(d2D)==sizeof(double*)+7*sizeof(std::size_t) );
//	#endif
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
	double arr[5] = {150, 16, 17, 18, 19}; //  NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
	BOOST_REQUIRE( &f2(arr) == &arr[2] );
}

