// Copyright 2019-2024 Alfredo A. Correa

// #define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS traits"
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
