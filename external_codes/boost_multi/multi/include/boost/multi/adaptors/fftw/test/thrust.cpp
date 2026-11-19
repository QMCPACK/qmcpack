// © Alfredo A. Correa 2020-2024

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi FFTW adaptor (cpu) with thrust complex"
#define BOOST_TEST_DYN_LINK

// #include<boost/test/unit_test.hpp>

#include "../../fftw.hpp"

#include<complex>
#include <thrust/complex.h>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(const fftw_2D_identity){
	using complex = thrust::complex<double>; complex const I{0.0, 1.0};

	multi::array<complex, 2> const in = {
		{  1.0 + 2.0*I,  9.0 - 1.0*I, 2.0 +  4.0*I},
		{  3.0 + 3.0*I,  7.0 - 4.0*I, 1.0 +  9.0*I},
		{  4.0 + 1.0*I,  5.0 + 3.0*I, 2.0 +  4.0*I},
		{  3.0 - 1.0*I,  8.0 + 7.0*I, 2.0 +  1.0*I},
		{ 31.0 - 1.0*I, 18.0 + 7.0*I, 2.0 + 10.0*I}
	};
	auto fwd = multi::fftw::dft({true, true}, in, multi::fftw::forward);
	
	multi::array<thrust::complex<double>, 2> const in_t = in;

	auto fwd_t = multi::fftw::dft({true, true}, in_t, multi::fftw::forward);
	
	BOOST_REQUIRE( fwd == fwd_t );
}
