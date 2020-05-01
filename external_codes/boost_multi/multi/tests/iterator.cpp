#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi allocators"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include<iostream>
#include<vector>
#include "../array.hpp"

namespace multi = boost::multi;

template<class MA>
decltype(auto) take(MA&& ma){return ma[0];}

BOOST_AUTO_TEST_CASE(iterator_1d){
	{
		multi::array<double, 1> A(100, 99.); 
		BOOST_REQUIRE( size(A) == 100 );
		BOOST_REQUIRE( begin(A) < end(A) );
		BOOST_REQUIRE( end(A) - begin(A) == size(A) );
		
		multi::array<double, 1>::const_iterator cb = cbegin(A);
		multi::array<double, 1>::iterator b = begin(A);
		BOOST_REQUIRE( cb == b );
		multi::array<double, 1>::const_iterator cb2 = begin(A);
		BOOST_REQUIRE( cb2 == cb );
	}
	{
		multi::array<double, 1> A({100}, 99.); 
		BOOST_REQUIRE( size(A) == 100 );
		BOOST_REQUIRE( begin(A) < end(A) );
	}
}

BOOST_AUTO_TEST_CASE(iterator_2d){
	{
		multi::array<double, 2> A({120, 140}, 99.); 
		BOOST_REQUIRE( size(A) == 120 );
		BOOST_REQUIRE( cbegin(A) < cend(A) );
		BOOST_REQUIRE( cend(A) - cbegin(A) == size(A) );
	}
	{
		std::vector<double> v(10000);
		multi::array_ref<double, 2> A(v.data(), {100, 100}); 
		BOOST_REQUIRE(size(A) == 100);
		begin(std::move(A))[4][3] = 2.; // ok 
		using multi::static_array_cast;
	}
}

BOOST_AUTO_TEST_CASE(iterator_reverse){
	multi::array<double, 3>::reverse_iterator rit;
	BOOST_REQUIRE(( rit.base() == multi::array<double, 3>::reverse_iterator{}.base() ));
	BOOST_REQUIRE(( multi::array<double, 3>::reverse_iterator{}.base() == multi::array<double, 3>::reverse_iterator{}.base() ));
	BOOST_REQUIRE(( multi::array<double, 3>::reverse_iterator{} == multi::array<double, 3>::reverse_iterator{} ));
	BOOST_REQUIRE(( multi::array<double, 3>::reverse_iterator{} == multi::array<double, 3>::reverse_iterator{} ));
}

BOOST_AUTO_TEST_CASE(iterator_interface){
	multi::array<double, 3> A =
		{
			{
				{ 1.2,  1.1}, { 2.4, 1.}
			},
			{
				{11.2,  3.0}, {34.4, 4.}
			},
			{
				{ 1.2,  1.1}, { 2.4, 1.}
			}
		}
	;

	BOOST_REQUIRE( size(A)==3 and size(A[0])==2 and size(A[0][0])==2);
	BOOST_REQUIRE( A[0][0][1] == 1.1 );

	BOOST_REQUIRE( begin(A) < end(A) );
	BOOST_REQUIRE( cbegin(A) < cend(A) );
	BOOST_REQUIRE( begin(A[0]) < end(A[0]) );
	BOOST_REQUIRE( begin(A[0]) < end(A[0]) );

	BOOST_REQUIRE(( multi::array<double, 3>::reverse_iterator{A.begin()} == rend(A) ));

	BOOST_REQUIRE( rbegin(A) < rend(A) );

	BOOST_REQUIRE( end(A) - begin(A) == size(A) );
	BOOST_REQUIRE( rend(A) - rbegin(A) == size(A) );

	BOOST_REQUIRE( size(*begin(A)) == 2 );
	BOOST_REQUIRE( size(begin(A)[1]) == 2 );

	BOOST_REQUIRE( &(A[1][1].begin()[0]) == &A[1][1][0] );
	BOOST_REQUIRE( &A[0][1][0] == &A[0][1][0] );
	BOOST_REQUIRE( &((*A.begin())[1][0]) == &A[0][1][0] );
	BOOST_REQUIRE( &((*A.begin()).operator[](1)[0]) == &A[0][1][0] );
	BOOST_REQUIRE( &(A.begin()->operator[](1)[0]) == &A[0][1][0] );
	BOOST_REQUIRE( &(A.begin()->operator[](1).begin()[0]) == &A[0][1][0] );
	BOOST_REQUIRE( &((A.begin()+1)->operator[](1).begin()[0]) == &A[1][1][0] );
	BOOST_REQUIRE( &((begin(A)+1)->operator[](1).begin()[0]) == &A[1][1][0] );
	BOOST_REQUIRE( &((cbegin(A)+1)->operator[](1).begin()[0]) == &A[1][1][0] );

}

BOOST_AUTO_TEST_CASE(iterator_semantics){

	multi::array<double, 3> A =
		{
			{
				{ 1.2,  1.1}, { 2.4, 1.}
			},
			{
				{11.2,  3.0}, {34.4, 4.}
			},
			{
				{ 1.2,  1.1}, { 2.4, 1.}
			}
		}
	;

	multi::array<double, 3>::iterator it; 
	BOOST_REQUIRE(( it == multi::array<double, 3>::iterator{} ));

	--it;
	it = begin(A);
	BOOST_REQUIRE( it == begin(A) );

	multi::array<double, 3>::iterator it2 = begin(A); 
	BOOST_REQUIRE(it == it2);

	it = end(A);
	BOOST_REQUIRE(it != it2);
	BOOST_REQUIRE(it > it2);

	multi::array<double, 3>::iterator it3{it};
	BOOST_REQUIRE( it3 == it );

	multi::array<double, 3>::const_iterator cit;
	cit = it3;
	BOOST_REQUIRE( cit == it3 );

	BOOST_REQUIRE((begin(A) == multi::array<double, 3>::iterator{rend(A)}));
}

