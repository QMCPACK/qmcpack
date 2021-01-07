#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX -std=c++17 $0 -o $0x `pkg-config --libs blas` -lboost_unit_test_framework -lcudart &&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS numeric"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../blas.hpp"
#include "../../blas/cuda.hpp"

#include "../../blas/numeric.hpp"
#include "../../../adaptors/cuda.hpp"

#include<complex>
#include<cassert>

using std::cout;
namespace multi = boost::multi;

using complex = std::complex<double>; constexpr complex I{0, 1};

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_numeric_imag){
	namespace blas = multi::blas;
	multi::array<complex, 1> a = { 1. + 2.*I, 3. + 5.*I, 9. + 2.*I };
	BOOST_REQUIRE( blas::imag(a)[2] == 2. );
}

namespace cuda = multi::cuda;

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_numeric_imag_cuda){
	cuda::array<complex, 1> a = { 1. + 2.*I, 3. + 5.*I, 9. + 2.*I };
	using multi::blas::imag;
	imag(a);
	BOOST_REQUIRE( imag(a)[2] == 2. );
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_numeric_imag_cuda_managed){
	cuda::managed::array<complex, 1> a = { 1. + 2.*I, 3. + 5.*I, 9. + 2.*I };
	using multi::blas::imag;
	BOOST_REQUIRE( imag(a)[2] == 2. );
}

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_test_numeric_hermitized_cuda){
	cuda::array<complex, 2> const a = {
		{ 1. + 2.*I, 3. + 5.*I, 9. + 2.*I },
		{ 1. + 2.*I, 3. + 5.*I, 9. + 2.*I },
		{ 1. + 2.*I, 3. + 5.*I, 9. + 2.*I },
	};
	using multi::blas::hermitized;
	hermitized(a);
}

