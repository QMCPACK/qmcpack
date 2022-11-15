#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXXX $CXXFLAGS $0 -o $0x `pkg-config --libs blas` -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2021

#ifndef MULTI_ADAPTORS_BLAS_HPP
#define MULTI_ADAPTORS_BLAS_HPP

#include "../adaptors/blas/asum.hpp"
#include "../adaptors/blas/axpy.hpp"
#include "../adaptors/blas/copy.hpp"
#include "../adaptors/blas/dot.hpp"
#include "../adaptors/blas/gemm.hpp"
#include "../adaptors/blas/gemv.hpp"
//#include "../adaptors/blas/ger.hpp"
#include "../adaptors/blas/herk.hpp"
#include "../adaptors/blas/iamax.hpp"
#include "../adaptors/blas/nrm2.hpp"
#include "../adaptors/blas/scal.hpp"
#include "../adaptors/blas/swap.hpp"
#include "../adaptors/blas/syrk.hpp"
#include "../adaptors/blas/trsm.hpp"

#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"
#include "../utility.hpp"

#include<iostream>
#include<complex>
#include<numeric> // iota
#include<algorithm> // transform

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_blas_herk_complex){
	using complex = std::complex<double>; complex const I{0, 1};
	using multi::blas::herk;
	{
		multi::array<complex, 2> const A = {
			{1. + 3.*I, 9. + 1.*I}, 
			{3. - 2.*I, 7. - 8.*I}, 
			{4. + 1.*I, 1. - 3.*I}
		};
		multi::array<complex, 2> C({3, 3}, 9999.);
		herk(1., A, C); // herk(A, C); // C†=C=AA†=(A†A)†
		BOOST_REQUIRE( C[1][2] == complex(41., 2.) );
		BOOST_REQUIRE( C[2][1] == conj(C[1][2]) );
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_asum_complex){
	using complex = std::complex<double>;
	multi::array<complex, 1> arr(1000, 0.);
//	std::iota(begin(arr), end(arr), -700.);
//	std::transform(cbegin(arr), cend(arr), begin(arr), [](auto&& a){return sqrt(a);});
	{
		using multi::blas::asum;
		BOOST_REQUIRE( asum(arr) == 0 );
	//	std::cout << asum(arr) << std::endl;
	}
}

BOOST_AUTO_TEST_CASE(multi_blas_nrm2_complex){
	multi::array<complex, 1> arr(1000, 0.);
//	std::iota(begin(arr), end(arr), -700.);
//	std::transform(cbegin(arr), cend(arr), begin(arr), [](auto&& a){return sqrt(a);});
	{
		using multi::blas::nrm2;
		BOOST_REQUIRE( nrm2(arr) == 0. );
	}
}

#endif
#endif

