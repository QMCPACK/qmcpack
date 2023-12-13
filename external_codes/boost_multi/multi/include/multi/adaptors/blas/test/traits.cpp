// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa


#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS traits"
#include<boost/test/unit_test.hpp>

#include "../../blas/traits.hpp"

#include<complex>

namespace multi = boost::multi;
namespace blas = multi::blas;

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_traits) {
	static_assert( blas::is_d<double>{} );
	static_assert( blas::is_s<float >{} );

	static_assert( blas::is_c<std::complex<float>>{} );
	static_assert( blas::is_z<std::complex<double>>{} );
}

#if 0
#if CUDA_FOUND  // TODO(correaa) move test to thrust adaptor
#include<thrust/complex.h>
BOOST_AUTO_TEST_CASE(multi_adaptors_blas_traits_thrust) {
	static_assert( blas::is_c<thrust::complex<float>>{} );
	static_assert( blas::is_z<thrust::complex<double>>{} );
}
#endif
#endif
