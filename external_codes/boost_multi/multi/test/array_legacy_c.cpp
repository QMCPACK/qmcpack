#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi legacy adaptor example"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include<iostream>

#include "../array_ref.hpp"
#include "../array.hpp"

#include<iostream>
#include<vector>
#include<complex>

namespace multi = boost::multi;
using std::cout; using std::cerr;

namespace fake{
typedef double fftw_complex[2];
void fftw_plan_dft(
	int rank, const int *n, 
	fftw_complex *in, fftw_complex *out, int sign, unsigned flags){
	(void)rank, (void)n, (void)in, (void)out, (void)sign, (void)flags;
}
}

BOOST_AUTO_TEST_CASE(array_legacy_c){
	using std::complex;

	multi::array<complex<double>, 2> const in = {
		{150., 16., 17., 18., 19.},
		{  5.,  5.,  5.,  5.,  5.}, 
		{100., 11., 12., 13., 14.}, 
		{ 50.,  6.,  7.,  8.,  9.}  
	};
	multi::array<std::complex<double>, 2> out(extensions(in));

	assert( dimensionality(out) == dimensionality(in) );
	assert( sizes(out) == sizes(in) );

	using multi::sizes_as;
	fake::fftw_plan_dft(
		dimensionality(in), sizes_as<int>(in).data(),
		(fake::fftw_complex*)in.data_elements(), (fake::fftw_complex*)out.data_elements(), 1, 0
	);

struct basic : multi::layout_t<2>{
	double* p;
};

struct ref : basic{
};


	{
		multi::array<double, 2> d2D = 
			{
				{150, 16, 17, 18, 19},
				{ 30,  1,  2,  3,  4}, 
				{100, 11, 12, 13, 14}, 
				{ 50,  6,  7,  8,  9} 
			};

	#if __has_cpp_attribute(no_unique_address) >=201803 and not defined(__NVCC__)
		BOOST_REQUIRE( sizeof(d2D)==sizeof(double*)+6*sizeof(std::size_t) );
	#endif
		BOOST_REQUIRE( d2D.is_compact() );
		BOOST_REQUIRE( rotated(d2D).is_compact() );
		BOOST_REQUIRE( d2D[3].is_compact() );
		BOOST_REQUIRE( not rotated(d2D)[2].is_compact() );
	}
	{
		using complex = std::complex<double>;
		multi::array<complex, 2> d2D({5, 3});
		BOOST_REQUIRE( d2D.is_compact() );
		BOOST_REQUIRE( rotated(d2D).is_compact() );
		BOOST_REQUIRE( d2D[3].is_compact() );
		BOOST_REQUIRE( not rotated(d2D)[2].is_compact() );
	}

}

