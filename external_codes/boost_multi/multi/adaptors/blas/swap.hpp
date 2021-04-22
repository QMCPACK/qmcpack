#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x `pkg-config --libs blas` -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020
#ifndef MULTI_ADAPTORS_BLAS_SWAP_HPP
#define MULTI_ADAPTORS_BLAS_SWAP_HPP

#include "../blas/core.hpp"

namespace boost{
namespace multi{
namespace blas{

template<class It1, class It2>
It2 swap(It1 first, It2 last, It2 first2){
	assert(stride(first) == stride(last));
	using std::distance;
	auto d = distance(first, last);
	swap(d, base(first), stride(first), base(first2), stride(first2));
	return first2 + d;
}

template<class X1D, class Y1D>
Y1D&& swap(X1D&& x, Y1D&& y){
	assert( size(x) == size(y) );
	assert( offset(x) == 0 and offset(y) == 0 );
	swap( begin(x), end(x), begin(y) );
	return std::forward<Y1D>(y);
}

}}}

#if not __INCLUDE_LEVEL__ // _TEST_MULTI_ADAPTORS_BLAS_SWAP

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS swap"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../array.hpp"
#include "../../utility.hpp"

#include "../blas/dot.hpp"

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_blas_swap, *boost::unit_test::tolerance(0.00001) ){
	multi::array<double, 2> A = {
		{1.,  2.,  3.,  4.},
		{5.,  6.,  7.,  8.},
		{9., 10., 11., 12.}
	};
//	using multi::blas::swap;
	multi::blas::swap(rotated(A)[1], rotated(A)[3]); // can ambiguate with (friend) multi::swap
	BOOST_REQUIRE( A[0][1] == 4. );
	BOOST_REQUIRE( A[0][3] == 2. );
}

#endif
#endif

