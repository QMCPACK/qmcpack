// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
// Copyright 2019-2021 Alfredo A. Correa

// TODO(correaa) : make it work with thrust complex

#ifndef MULTI_ADAPTORS_BLAS_ASUM_HPP
#define MULTI_ADAPTORS_BLAS_ASUM_HPP

#include "../blas/core.hpp"

namespace boost::multi::blas {

template<class It, typename Size>
auto asum_n(It first, Size n)
->decltype(asum(n, base(first), stride(first))) {
	return asum(n, base(first), stride(first)); }

using std::distance;

template<class It>
auto asum(It f, It last)
->decltype(asum_n(f, distance(f, last))) {assert(stride(f) == stride(last));
	return asum_n(f, distance(f, last)); }

using std::begin; using std::end;

template<class X1D>
auto asum(X1D const& x)
->decltype(asum(begin(x), end(x))) {assert( not offset(x) );
	return asum(begin(x), end(x)); }

} // end namespace boost::multi::blas

//#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

//#define BOOST_TEST_MODULE "C++ Unit Tests for Multi.BLAS asum"
//#define BOOST_TEST_DYN_LINK
//#include<boost/test/unit_test.hpp>
////#include<boost/test/tools/floating_point_comparison.hpp>

//#include "../../array.hpp"
////#include "../../utility.hpp"

//#include<numeric> // accumulate

//namespace multi = boost::multi;
//using multi::blas::asum;

//BOOST_AUTO_TEST_CASE(multi_blas_asum_double){
//	multi::array<double, 2> const A = {
//		{1.,  2.,  3.,  4.},
//		{-5.,  6.,  -7.,  8.},
//		{9., 10., 11., 12.}
//	};
//	BOOST_REQUIRE(asum(A[1]) == std::accumulate(begin(A[1]), end(A[1]), 0., [](auto&& a, auto&& b){return a+std::abs(b);}));
//}

//BOOST_AUTO_TEST_CASE(multi_blas_asum_complex){

//	using complex = std::complex<double>; complex const I{0, 1};
//	multi::array<complex, 2> const A = {
//		{ 1. + 1.*I,  2.,  3.,  4.},
//		{-5. + 3.*I,  6.,  -7.,  8.},
//		{ 9. - 2.*I, 10., 11., 12.}
//	};
//	BOOST_REQUIRE(asum(rotated(A)[0]) == 1.+1. + 5.+3. + 9.+2.);

//}

//BOOST_AUTO_TEST_CASE(multi_blas_asum_double_carray){
////	double A[3][4] = {
////		{1.,  2.,  3.,  4.},
////		{-5.,  6.,  -7.,  8.},
////		{9., 10., 11., 12.}
////	}; (void)A;
////	using std::begin; using std::end;
////	BOOST_REQUIRE(asum(A[1]) == std::accumulate(begin(A[1]), end(A[1]), 0., [](auto&& a, auto&& b){return a+abs(b);}));
//}

//#endif
#endif
