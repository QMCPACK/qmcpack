#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXXX $CXXFLAGS $0 -o $0.$X `pkg-config --cflags --libs blas cuda-11.0` -lboost_unit_test_framework&&$0.$X&&rm $0.$X;exit
#endif

#include "../../blas/traits.hpp"

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS traits"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "./config.hpp"

#include<complex>

namespace multi = boost::multi;
namespace blas = multi::blas;

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_traits) {
	static_assert( blas::is_d<double>{} );
	static_assert( blas::is_s<float >{} );

	static_assert( blas::is_c<std::complex<float>>{} );
	static_assert( blas::is_z<std::complex<double>>{} );
}

#if CUDA_FOUND
#include<thrust/complex.h>
BOOST_AUTO_TEST_CASE(multi_adaptors_blas_traits_thrust) {
	static_assert( blas::is_c<thrust::complex<float>>{} );
	static_assert( blas::is_z<thrust::complex<double>>{} );
}
#endif

