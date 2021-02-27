#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $0 -o $0x -lfftw3 -lboost_unit_test_framework -ftemplate-backtrace-limit=0&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi FFTW adaptor (cpu) with thrust complex"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../fftw.hpp"

#include<complex>
#include <thrust/complex.h>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(fftw_2D_identity){

	using complex = thrust::complex<double>; complex const I{0, 1};

	multi::array<complex, 2> const in = {
		{  1. + 2.*I,  9. - 1.*I, 2. +  4.*I},
		{  3. + 3.*I,  7. - 4.*I, 1. +  9.*I},
		{  4. + 1.*I,  5. + 3.*I, 2. +  4.*I},
		{  3. - 1.*I,  8. + 7.*I, 2. +  1.*I},
		{ 31. - 1.*I, 18. + 7.*I, 2. + 10.*I}
	};
	auto fwd = multi::fftw::dft({true, true}, in, multi::fftw::forward);
	
	multi::array<thrust::complex<double>, 2> const in_t = in;

	auto fwd_t = multi::fftw::dft({true, true}, in_t, multi::fftw::forward);
	
	BOOST_REQUIRE( fwd == fwd_t );

}

