#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXXX $CXXFLAGS $0 -o $0.$X `pkg-config --libs blas` -lboost_unit_test_framework&&$0.$X&&rm $0.$X;exit
#endif
// © Alfredo A. Correa 2019-2020

#ifndef MULTI_ADAPTORS_BLAS_SIDE_HPP
#define MULTI_ADAPTORS_BLAS_SIDE_HPP

#include    "../blas/core.hpp"
#include    "../blas/operations.hpp"
#include "../../array_ref.hpp"

namespace boost{
namespace multi{
namespace blas{

//enum class SIDE : char{L='L', R='R'};

enum side : char{
	left  = 'L', 
	right = 'R'//,
//	pre_multiply = 'R', 
//	post_multiply = 'L'
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

#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi BLAS adaptors side"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(multi_adaptors_blas_side){
	return;
}

#endif
#endif

