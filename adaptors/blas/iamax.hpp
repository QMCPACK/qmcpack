#ifdef COMPILATION_INSTRUCTIONS
(echo '#include"'$0'"'>$0.cpp)&&c++ -Wall -Wextra -Wpedantic -Wfatal-errors -D_TEST_MULTI_ADAPTORS_BLAS_IAMAX $0.cpp -o $0x `pkg-config --libs blas` -lboost_unit_test_framework&&$0x&&rm $0x $0.cpp;exit
#endif
// © Alfredo A. Correa 2019-2020
#ifndef MULTI_ADAPTORS_BLAS_IAMAX_HPP
#define MULTI_ADAPTORS_BLAS_IAMAX_HPP

#include "../blas/core.hpp"

namespace boost{
namespace multi{
namespace blas{

template<class It, class Size>
auto iamax_n(It first, Size n){
	using core::iamax;
	return iamax(n, base(first), stride(first)); 
	// if you get an error here make sure that you are including (and linking) the appropriate BLAS backend for your memory type
}

template<class It>
auto iamax(It first, It last)
->decltype(iamax_n(first, std::distance(first, last))){
	return iamax_n(first, std::distance(first, last));}

template<class X1D>
auto iamax(X1D const& x)
->decltype(iamax(begin(x), end(x))){assert( not offset(x) );
	return iamax(begin(x), end(x));
}

template<class X1D> auto amax(X1D const& x){return begin(x) + iamax(x);}

}}}

#if _TEST_MULTI_ADAPTORS_BLAS_IAMAX

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS iamax"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../array.hpp"
#include "../../utility.hpp"

#include<complex>
#include<cassert>

using std::cout;
namespace multi = boost::multi;
namespace blas = multi::blas;

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_iamax_real){
	multi::array<double, 1> const A = {1., 2., 3., 4.};
	using blas::iamax;
	auto i = iamax(A);
	BOOST_REQUIRE( i == 3 );
	BOOST_REQUIRE( A[iamax(A)] == 4. );
	using blas::amax;
	BOOST_REQUIRE( *amax(A) == 4. );
}

using complex = std::complex<double>;

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_iamax_complex){
	multi::array<complex, 1> const A = {1., 2., 3., 4.};
	using blas::iamax;
	auto i = iamax(A);
	BOOST_REQUIRE( i == 3 );
	BOOST_REQUIRE( A[iamax(A)] == 4. );
	BOOST_REQUIRE( *blas::amax(A) == 4. );
}

#endif
#endif

