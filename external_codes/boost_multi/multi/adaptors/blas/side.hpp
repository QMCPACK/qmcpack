#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0x `pkg-config --libs blas` -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2019

#ifndef MULTI_ADAPTORS_BLAS_SIDE_HPP
#define MULTI_ADAPTORS_BLAS_SIDE_HPP

#include    "../blas/core.hpp"
#include    "../blas/operations.hpp"
#include "../../array_ref.hpp"

namespace boost{
namespace multi{
namespace blas{

enum class SIDE : char{L='L', R='R'};

enum side : char{
	left = static_cast<char>(SIDE::R), right = static_cast<char>(SIDE::L),
	pre_multiply = static_cast<char>(SIDE::R), post_multiply = static_cast<char>(SIDE::L)
};

side swap(side s){
	switch(s){
		case side::left: return side::right;
		case side::right: return side::left;
	} __builtin_unreachable();
}

}}}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#if not __INCLUDE_LEVEL__ // _TEST_MULTI_ADAPTORS_BLAS_SIDE

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS adaptors side"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../../array.hpp"
#include "../../utility.hpp"
#include "../blas/nrm2.hpp"

#include<complex>
#include<cassert>
#include<iostream>
#include<numeric>
#include<algorithm>

using std::cout;

template<class M> 
decltype(auto) print(M const& C){
	using boost::multi::size;
	for(int i = 0; i != size(C); ++i){
		for(int j = 0; j != size(C[i]); ++j) cout<< C[i][j] <<' ';
		cout<<std::endl;
	}
	return cout<<"---"<<std::endl;
}

namespace multi = boost::multi;
using complex = std::complex<double>; constexpr complex I{0, 1};

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_side){
	return;
}

#endif
#endif

