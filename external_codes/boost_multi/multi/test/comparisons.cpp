#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0.$X -lboost_unit_test_framework&&$0.$X&&rm $0.$X;exit
#endif
//  © Alfredo A. Correa 2018-2019

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi comparisons"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

#include<complex>

BOOST_AUTO_TEST_CASE(comparisons){
	namespace multi = boost::multi;

	using complex = std::complex<double>;
	{
		multi::array<double, 1> A = {1., 2., 3.};
		multi::array<complex, 1> B = {1., 2., 3.};
		BOOST_REQUIRE( A[1] == B[1] );
		BOOST_REQUIRE( A == B );
		BOOST_REQUIRE( B == A );
	}
	{
		multi::array<double, 2> const A = {{1., 2., 3.}, {4.,5.,6.}};
		multi::array<complex, 2> const B = {{1., 2., 3.}, {4.,5.,6.}};
		BOOST_REQUIRE( A[1][1] == B[1][1] );
		BOOST_REQUIRE( A == B );
		BOOST_REQUIRE( B == A );
		BOOST_REQUIRE( std::equal(A[1].begin(), A[1].end(), begin(B[1])) );
	}

	multi::array<double, 3> A =
//		#if __INTEL_COMPILER
//		(double[3][2][2])
//		#endif
		{
			{{ 1.2,  1.1}, { 2.4, 1.}},
			{{11.2,  3.0}, {34.4, 4.}},
			{{ 1.2,  1.1}, { 2.4, 1.}}
		};
	
	BOOST_REQUIRE( A[0] < A[1] );

	swap( A[0], A[1] );
	swap( A[0], A[1] );
	
	multi::array_ref<double, 3>  AR(A.data_elements(), extensions(A));
	multi::array_cref<double, 3> AC(data_elements(A), extensions(A));

	BOOST_REQUIRE( A == A );
	BOOST_REQUIRE( AR == A );
	BOOST_REQUIRE( AR == AC );

	BOOST_REQUIRE( A[0] == A[2] );
	BOOST_REQUIRE( AR[0] == A[2] );
	BOOST_REQUIRE( AR[0] == AC[2] );

	BOOST_REQUIRE( not (A[0] != A[2]) );
	BOOST_REQUIRE( not (AR[0] != AR[2]) );

	BOOST_REQUIRE( not (A[0] != A[2]) );
	BOOST_REQUIRE( not (AR[0] != AR[2]) );

	BOOST_REQUIRE( A[0] <= A[1] );
	BOOST_REQUIRE( AR[0] <= A[1] );
	BOOST_REQUIRE( AC[0] <= AC[1] );
	BOOST_REQUIRE( A[0][0] <= A[0][1] );
	BOOST_REQUIRE( AR[0][0] <= A[0][1] );

	BOOST_REQUIRE_EQUAL( A[1][0][0] , 11.2 );
	BOOST_REQUIRE_EQUAL( AR[1][0][0] , 11.2 );
	BOOST_REQUIRE_EQUAL( AC[1][0][0] , 11.2 );

	BOOST_REQUIRE_EQUAL( A[0][0][0] , 1.2 );
	BOOST_REQUIRE_EQUAL( AR[0][0][0] , 1.2 );
	BOOST_REQUIRE_EQUAL( AC[0][0][0] , 1.2 );

	swap(AR[0], AR[1]);
//	using std::begin;
//	using std::end;

	BOOST_REQUIRE( begin(A) < end(A) );
	BOOST_REQUIRE( cbegin(A) < cend(A) );
//	BOOST_REQUIRE( crbegin(A) < crend(A) );
//	assert( crend(A) > crbegin(A) );
	BOOST_REQUIRE( end(A) - begin(A) == size(A));
	BOOST_REQUIRE( rend(A) - rbegin(A) == size(A));

}

