#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
echo $X;
mkdir -p build.$X && cd build.$X && cmake .. && make gemv.cppx && ctest;exit
#endif
// © Alfredo A. Correa 2020

// $CXX -D_MULTI_CUBLAS_ALWAYS_SYNC $0 -o $0x `pkg-config --libs blas` -lcudart -lcublas -lboost_unit_test_framework&&$0x&&rm $0x; exit

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS gemv"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../../adaptors/blas/gemv.hpp"
#include "../../../array.hpp"

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_blas_gemv){

	multi::array<double, 2> const M = {
		{ 9., 24., 30., 9.},
		{ 4., 10., 12., 7.},
		{14., 16., 36., 1.}
	};
	multi::array<double, 1> const X = {1.1,2.1,3.1, 4.1};

}

