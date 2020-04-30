#ifdef COMPILATION_INSTRUCTIONS
(echo '#include"'$0'"'>$0.cpp)&&clang++ -D_TEST_MULTI_ADAPTORS_BLAS_SCAL $0.cpp -o $0x `pkg-config --libs blas` -lboost_unit_test_framework&&$0x&&rm $0x $0.cpp;exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef MULTI_ADAPTORS_BLAS_SCAL_HPP
#define MULTI_ADAPTORS_BLAS_SCAL_HPP

#include "../blas/core.hpp"
#include "../../config/NODISCARD.hpp"

namespace boost{namespace multi{
namespace blas{

using core::scal;

template<class X1D>
auto scal(typename std::decay_t<X1D>::element_type a, X1D&& m)
->decltype(scal(size(m), &a, base(m), stride(m)), std::forward<X1D>(m)){
	return scal(size(m), &a, base(m), stride(m)), std::forward<X1D>(m);}

template<class X1D>
NODISCARD("because last argument is const")
auto scal(typename X1D::element_type a, X1D const& m){
	return scal(a, m.decay());
}

}}}

#if _TEST_MULTI_ADAPTORS_BLAS_SCAL

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

		using blas::scal;
		auto S = scal(2., rotated(A)[1]);

		BOOST_REQUIRE( A[2][1] == 20 );
		BOOST_REQUIRE( S[0] == 4 );
	}
	{
		multi::array<double, 2> const A = {
			{1.,  2.,  3.,  4.},
			{5.,  6.,  7.,  8.},
			{9., 10., 11., 12.}
		};
		using multi::blas::scal;
		auto rA1_scaled = scal(2., A[1]);
		BOOST_REQUIRE( size(rA1_scaled) == 4 );
		BOOST_REQUIRE( rA1_scaled[1] == 12 );
	}

}

#endif
#endif

