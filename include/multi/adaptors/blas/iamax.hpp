// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
// Copyright 2019-2021 Alfredo A. Correa
#ifndef MULTI_ADAPTORS_BLAS_IAMAX_HPP
#define MULTI_ADAPTORS_BLAS_IAMAX_HPP

#include "../blas/core.hpp"

namespace boost::multi::blas {

template<class It, class Size>
auto iamax_n(It first, Size n) {
	using core::iamax;
	return iamax(n, base(first), stride(first)); 
	// if you get an error here make sure that you are including (and linking) the appropriate BLAS backend for your memory type
}

template<class It>
auto iamax(It first, It last)
->decltype(iamax_n(first, std::distance(first, last))) {
	return iamax_n(first, std::distance(first, last)); }

template<class X1D>
auto iamax(X1D const& x)
->decltype(iamax(begin(x), end(x))) {assert( not offset(x) );
	return iamax(begin(x), end(x)); }

template<class X1D> auto amax(X1D const& x) {return begin(x) + iamax(x);}

} // end namespace boost::multi::blas

//#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

//#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS iamax"
//#define BOOST_TEST_DYN_LINK
//#include<boost/test/unit_test.hpp>

//#include "../../array.hpp"
//#include "../../utility.hpp"

//#include<complex>
//#include<cassert>

//using std::cout;
//namespace multi = boost::multi;
//namespace blas = multi::blas;

//BOOST_AUTO_TEST_CASE(multi_adaptors_blas_iamax_real){
//	multi::array<double, 1> const A = {1., 2., 3., 4.};

//	auto i = blas::iamax(A);
//	BOOST_REQUIRE( i == 3 );
//	BOOST_REQUIRE( A[blas::iamax(A)] == 4. );

//	BOOST_REQUIRE( *blas::amax(A) == 4. );
//}

//using complex = std::complex<double>;

//BOOST_AUTO_TEST_CASE(multi_adaptors_blas_iamax_complex){
//	multi::array<complex, 1> const A = {1., 2., 3., 4.};
//	auto i = blas::iamax(A);
//	BOOST_REQUIRE( i == 3 );
//	BOOST_REQUIRE( A[blas::iamax(A)] == 4. );
//	BOOST_REQUIRE( *blas::amax(A) == 4. );
//}

//#endif
#endif
