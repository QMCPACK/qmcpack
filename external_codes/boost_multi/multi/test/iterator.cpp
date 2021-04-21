#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $CXXFLAGS $0 -o $0.$X -lboost_unit_test_framework&&$0.$X $@&&rm $0.$X;exit
#endif
// © Alfredo A. Correa 2018-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi iterators"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

#include<iostream>
#include<numeric>
#include<vector>

namespace multi = boost::multi;

template<class MA> decltype(auto) take(MA&& ma){return ma[0];}

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
		
		auto b = A.begin();
		multi::array<double, 1>::const_iterator cbb = b;
		BOOST_REQUIRE( cbb == b );
		BOOST_REQUIRE( b == cbb );
	}
	
	*begin( multi::array<double, 1>(10, 99.) ) = 44.;
}

BOOST_AUTO_TEST_CASE(iterator_2d){
	{
		multi::array<double, 2> A({120, 140}, 99.); 
		BOOST_REQUIRE( size(A) == 120 );
		BOOST_REQUIRE( cbegin(A) < cend(A) );
		BOOST_REQUIRE( cend(A) - cbegin(A) == size(A) );

		using iter = multi::array<double, 2>::iterator;
		static_assert( std::is_same< iter::element   , double >{}, "!");
		static_assert( std::is_same< iter::value_type, multi::array<double, 1> >{}, "!");
		static_assert( std::is_same< iter::reference, multi::basic_array<double, 1>&&>{}, "!");
		static_assert( std::is_same< iter::element_ptr, double*>{}, "!");

		using citer = multi::array<double, 2>::const_iterator;
		static_assert( std::is_same< citer::element   , double >{}, "!");
		static_assert( std::is_same< citer::value_type, multi::array<double, 1> >{}, "!");
		static_assert( std::is_same< citer::reference, multi::basic_array<double, 1, double const*>&&>{}, "!");
		static_assert( std::is_same< citer::element_ptr, double const* >{}, "!");
	}
	{
		std::vector<double> v(10000);
		multi::array_ref<double, 2> A(v.data(), {100, 100}); 
		BOOST_REQUIRE(size(A) == 100);
		begin(A)[4][3] = 2.;
	}
}

BOOST_AUTO_TEST_CASE(iterator_reverse){
	multi::array<double, 3>::reverse_iterator rit = {};
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
			{{ 1.2,  1.1}, { 2.4, 1.}},
			{{11.2,  3.0}, {34.4, 4.}},
			{{ 1.2,  1.1}, { 2.4, 1.}}
		}
	;

	multi::array<double, 3>::iterator it; 
	BOOST_REQUIRE(( it == multi::array<double, 3>::iterator{} ));

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

	BOOST_REQUIRE(( begin(A) == multi::array<double, 3>::iterator(rend(A)) ));

}

BOOST_AUTO_TEST_CASE(iterator_arrow_operator){

	multi::array<std::string, 2> A = {
		{"00", "01"},
		{"10", "11"},
		{"20", "21"}
	};
	
	BOOST_REQUIRE( A[1][0] == "10" );

	BOOST_REQUIRE( std::is_sorted(begin(A), end(A)) ); // sorted by rows
	BOOST_REQUIRE( std::is_sorted(begin(A.rotated()), end(A.rotated())) ); // sorted by cols
	
	BOOST_REQUIRE( begin( A           )->size() == A[0].size() );
	BOOST_REQUIRE( begin( A.rotated() )->size() == A.size() );
	
	BOOST_REQUIRE( &(begin( A           )->operator[](1)) == &(A[0][1]) );
	BOOST_REQUIRE( &(begin( A.rotated() )->operator[](1)) == &(A[1][0]) );

}

BOOST_AUTO_TEST_CASE(index_range_iteration){
	multi::index_range r{0, 5}; // semiopen interval
	std::ostringstream out;
	std::copy(begin(r), end(r), std::ostream_iterator<int>{out, ","});
	BOOST_REQUIRE( out.str() == "0,1,2,3,4," );

	BOOST_REQUIRE( std::accumulate(begin(r), end(r), 0) == r.size()*(r.size()-1)/2 );

	BOOST_REQUIRE( std::accumulate(begin(r), end(r), 0, [](auto&& acc, auto const& e){return acc + e*e*e;}) > 0 ); // sum of cubes
}

