#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x `pkg-config --libs blas` -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef MULTI_ADAPTORS_BLAS_SCAL_HPP
#define MULTI_ADAPTORS_BLAS_SCAL_HPP

#include "../blas/core.hpp"
#include "../../config/NODISCARD.hpp"

namespace boost{namespace multi{
namespace blas{

using core::scal;

template<class X1D> // don't do this: ", typename Elem = typename X1D::element_type>"
auto scal(typename std::decay_t<X1D>::element_type a, X1D&& m)
->decltype(scal(size(m), &a, base(m), stride(m)), std::forward<X1D>(m)){
	return scal(size(m), &a, base(m), stride(m)), std::forward<X1D>(m);}

template<class X1D>//, typename Elem = typename X1D::element_type>
NODISCARD("because last argument is const")
typename X1D::decay_type scal(typename X1D::element_type a, X1D const& m){
	return scal(a, typename X1D::decay_type{m});
}

}}}

#if not __INCLUDE_LEVEL__ // _TEST_MULTI_ADAPTORS_BLAS_SCAL

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS scal"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../array.hpp"

namespace multi = boost::multi;
namespace blas = multi::blas;

BOOST_AUTO_TEST_CASE(multi_blas_scal_real){
	{
		multi::array<double, 2> A = {
			{1.,  2.,  3.,  4.},
			{5.,  6.,  7.,  8.},
			{9., 10., 11., 12.}
		};

		auto S = blas::scal(2., rotated(A)[1]);

		BOOST_REQUIRE( A[2][1] == 20 );
		BOOST_REQUIRE( S[0] == 4 );
	}
	{
		multi::array<double, 2> const A = {
			{1.,  2.,  3.,  4.},
			{5.,  6.,  7.,  8.},
			{9., 10., 11., 12.}
		};
		auto rA1_scaled = blas::scal(2., A[1]);
		BOOST_REQUIRE( A[1][1] == 6. );
		BOOST_REQUIRE( size(rA1_scaled) == 4 );
		BOOST_REQUIRE( rA1_scaled[1] == 12 );
	}

}

#endif
#endif

