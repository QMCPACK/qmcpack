#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2019

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi layout"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"
#include "../utility.hpp"

#include<iostream>
#include<tuple>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(linearize){
	multi::array<double, 3> A({10, 20, 30});
	BOOST_REQUIRE(  25 % extensions(A) == std::make_tuple(0, 0, 25) );
	BOOST_REQUIRE(  55 % extensions(A) == std::make_tuple(0, 1, 25) );
	BOOST_REQUIRE( 655 % extensions(A) == std::make_tuple(1, 1, 25) );
	BOOST_REQUIRE(1255 % extensions(A) == std::make_tuple(2, 1, 25) );
	
	std::tuple<multi::index, multi::index, multi::index> p = A.extensions().from_linear(655);
	BOOST_REQUIRE( p == std::make_tuple(1, 1, 25) );
}

BOOST_AUTO_TEST_CASE(layout){
{	
	multi::array<double, 3> A3(
#if defined(__INTEL_COMPILER) or (defined(__GNUC__) and (__GNUC__ < 6))
		multi::iextensions<3>
#endif
		{51, 52, 53}
	);
	BOOST_REQUIRE( size(A3) == 51      ); BOOST_REQUIRE( A3.size() == 51       );
	BOOST_REQUIRE( size(A3[0]) == 52   ); BOOST_REQUIRE( A3[0].size() == 52    );
	BOOST_REQUIRE( size(A3[0][0]) == 53); BOOST_REQUIRE( A3[0][0].size() == 53 );
}
{		
	double DA[50][50][50];
	using multi::size;
	BOOST_REQUIRE( size(DA) == 50 );

	using multi::extension;
	BOOST_REQUIRE(( extension(DA) == multi::index_extension{0, 50} ));
	BOOST_REQUIRE(( extension(DA) == multi::iextension{0, 50}      ));
	BOOST_REQUIRE(( extension(DA) == multi::irange{0, 50} ));
}
{
	multi::array<double, 2> B2(
#if defined(__INTEL_COMPILER) or (defined(__GNUC__) and (__GNUC__ < 6))
		multi::iextensions<2>
#endif
		{50, 50}
	);
	BOOST_REQUIRE( size(B2) == 50 ); BOOST_REQUIRE( B2.size() == 50 );
	BOOST_REQUIRE( B2[0].sliced(10, 20).size() == 10 );
	BOOST_REQUIRE( size(B2[0].sliced(10, 20))  == 10 );

	BOOST_REQUIRE( B2(0, {10, 20}).dimensionality  == 1 );
	BOOST_REQUIRE( dimensionality(B2(0, {10, 20})) == 1 );

	BOOST_REQUIRE( size(B2(0, {10, 20})) == 10 );
}
{
	multi::array<double, 2> A2 = 
//		#if defined(__INTEL_COMPILER)
//		(double[3][3])
//		#endif
		{{1., 2., 3.}, 
		 {4., 5., 6.}, 
		 {7., 8., 9.}}
	;
	multi::array<int, 2> B2(
#if defined(__INTEL_COMPILER) or (defined(__GNUC__) and (__GNUC__ < 6))
		multi::iextensions<2>
#endif
		{4, 4}
	);
	BOOST_REQUIRE( size(B2) == 4 );
	B2[3][3] = 99.;
	
//	decltype(B2({0,2}, {0,2}))::decay_type B2copy = B2({0,2}, {0,2});
	auto B2copy =+ B2({0,2}, {0,2});

	decltype(B2({0,2}, {0,2})) B2blk[2][2] = 
		{{B2({0,2}, {0,2}), B2({0, 2}, {2, 4})},
		 {B2({2,4}, {0,2}), B2({2, 4}, {2, 4})}};
	std::cout << B2blk[1][1][1][1] << std::endl;
	BOOST_REQUIRE( &B2blk[1][1][1][1] == &B2[3][3] );
}
{
	double A[3][4][5] = {};
	using multi::dimensionality;
	static_assert(dimensionality(A)==3, "!");
	using multi::extensions;
	auto xA = extensions(A);
	using std::get;
	BOOST_REQUIRE( size(get<0>(xA)) == 3 );
	BOOST_REQUIRE( size(get<1>(xA)) == 4 );
	BOOST_REQUIRE( size(get<2>(xA)) == 5 );

	multi::array<double, 3> AA({3, 4, 5});
	using multi::layout;
	assert( layout(AA) == layout(A) );

	using multi::stride;
	BOOST_REQUIRE( stride(AA) == 20 );
	static_assert( stride(A) == 20 , "!" );
	static_assert( stride(A[0]) == 5 , "!" );
	static_assert( stride(A[1]) == 5 , "!" );
	static_assert( stride(A[0][0]) == 1 , "!" );
//		assert( stride(A) == 20 );
//		assert( stride(A[0]) == 20 );
}
{
	multi::array<double, 2> B2 = {
		{1.},
		{2.},
		{3.}
	};
	BOOST_REQUIRE( size(B2) == 3 );
	BOOST_REQUIRE( size(rotated(B2)) == 1 ); BOOST_REQUIRE( size(B2[0]) == 1);
	BOOST_REQUIRE( stride(B2) == 1 );
	BOOST_REQUIRE( stride(B2[0]) == 1 );

}

}

