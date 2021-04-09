#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $CXXFLAGS $0 -o $0.$X -lboost_unit_test_framework&&$0.$X $@&&rm $0.$X;exit
#endif
// © Alfredo A. Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi initializer_list"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

#include<complex>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_tests_initializer_list_1d){
	{
		multi::static_array<double, 1> const A = {1.2, 3.4, 5.6};
		BOOST_REQUIRE( size(A) == 3 );
		BOOST_REQUIRE( A[2] == 5.6 );
	}
	{
		auto il = {1.2, 3.4, 5.6};
		multi::static_array<double, 1> const A(il);
		BOOST_REQUIRE( size(A) == 3 );
		BOOST_REQUIRE( A[2] == il.begin()[2] );
	}
	{
		auto il = {1.2, 3.4, 5.6};
		multi::static_array<double, 1> const A(begin(il), end(il));
		BOOST_REQUIRE( size(A) == 3 );
		BOOST_REQUIRE( A[2] == il.begin()[2] );
	}
	{
		multi::static_array<double, 1> const A = {1.2, 3.4, 5.6};
		BOOST_REQUIRE( size(A) == 3 );
		BOOST_REQUIRE( A[2] == 5.6 );
		BOOST_REQUIRE(( A == multi::static_array<double, 1>{1.2, 3.4, 5.6} ));
		BOOST_REQUIRE(( A == decltype(A){1.2, 3.4, 5.6} ));
	}
	{
	#if defined(__cpp_deduction_guides)
		multi::static_array const A = {1.2, 3.4, 5.6};
		BOOST_REQUIRE( size(A) == 3 );
		BOOST_REQUIRE( A[2] == 5.6 );
		BOOST_REQUIRE(( A == multi::static_array{1.2, 3.4, 5.6} ));
	#endif
	}
	{
		auto il = {1.2, 3.4, 5.6};
		multi::array<double, 1> const A(il.begin(), il.end());
		BOOST_REQUIRE( size(A) == 3 );
		BOOST_REQUIRE( A[2] == 5.6 );
	}
	{
		multi::array<double, 1> const A = {1.2, 3.4, 5.6};
		BOOST_REQUIRE( size(A) == 3 );
		BOOST_REQUIRE( A[2] == 5.6 );
		BOOST_REQUIRE(( A == multi::array<double, 1>{1.2, 3.4, 5.6} ));
		BOOST_REQUIRE(( A == decltype(A){1.2, 3.4, 5.6} ));
		BOOST_REQUIRE(( A == decltype(A)::decay_type({1.2, 3.4, 5.6}) ));
	}
	{
	#if defined(__cpp_deduction_guides)
		multi::array A({1.2, 3.4, 5.6});
		BOOST_REQUIRE( size(A) == 3 );
		BOOST_REQUIRE( A[2] == 5.6 );
		BOOST_REQUIRE(( A == multi::array({1.2, 3.4, 5.6}) ));
	#endif
	}
	{
		std::array<double, 3> const a = {1.1, 2.2, 3.3};
		using multi::num_elements;
		BOOST_REQUIRE( num_elements(a) == 3 );

		using std::begin; using std::end;
		multi::static_array<double, 1> const A(begin(a), end(a));
		BOOST_REQUIRE( size(A) == 3 );
	}
}

BOOST_AUTO_TEST_CASE(multi_tests_initializer_list_array){
//#if not defined (__GNUG__)
#if defined(__INTEL_COMPILER) or (defined(__clang__) and (__clang_major__ >= 10))  // doesn't work on gcc
	{
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc99-designator"
//		double const a[] = { [8] = 8., 9., 10. };
		std::array<double, 11> const a = {{ [8] = 8., 9., 10. }};
#pragma GCC diagnostic pop
		multi::array<double, 1> A = a;
		BOOST_REQUIRE( A.size() == 11 );
		BOOST_REQUIRE( A[9] == 9. );
	}
#endif
}

BOOST_AUTO_TEST_CASE(multi_initialize_from_carray_1d){
	{
		multi::static_array<double, 1> const A = {1.1, 2.2, 3.3};
		BOOST_REQUIRE( size(A) == 3 );
		BOOST_REQUIRE( A[1] == 2.2 );
	}
	{
#if defined(__cpp_deduction_guides)
		multi::array A = {1.1, 2.2, 3.3}; // warning: ISO C++ forbids compound-literals [-Wpedantic]
		BOOST_REQUIRE( size(A)==3 and A[1] == 2.2 );
#endif
	}
	{
		std::array<double, 3> a = {1.1, 2.2, 3.3};
		multi::array<double, 1> const A(begin(a), end(a));
		BOOST_REQUIRE(( A == decltype(A){1.1, 2.2, 3.3} ));
	}
	{
	#if defined(__cpp_deduction_guides)
		std::array a = {1.1, 2.2, 3.3};
		multi::array<double, 1> const A(begin(a), end(a));
		assert(( A == decltype(A){1.1, 2.2, 3.3} ));
	#endif
	}
}

BOOST_AUTO_TEST_CASE(multi_tests_initializer_list_2d){
	{
		multi::static_array<double, 2> const A = {
			{ 1.2,  2.4, 3.6, 8.9},
			{11.2, 34.4, 5.6, 1.1},
			{15.2, 32.4, 5.6, 3.4}
		};
		BOOST_REQUIRE( size(A) == 3 );
		BOOST_REQUIRE( size(A[0]) == 4 );
		BOOST_REQUIRE(( A == decltype(A)
			{
				{ 1.2,  2.4, 3.6, 8.9},
				{11.2, 34.4, 5.6, 1.1},
				{15.2, 32.4, 5.6, 3.4}
			}
		));
	}
	{
		multi::array<double, 2> const A = {
			{ 1.2,  2.4, 3.6},
			{11.2, 34.4, 5.6},
			{15.2, 32.4, 5.6}
		};
		BOOST_REQUIRE( size(A) == 3 );
		BOOST_REQUIRE( size(A) == 3 and size(A[0]) == 3 );
		BOOST_REQUIRE( A[1][1] == 34.4 );
	}
	{
		multi::array<double, 2> A = {
			{ 1.2,  2.4, 3.6},
			{11.2, 34.4, 5.6},
			{15.2, 32.4, 5.6}
		};
		BOOST_REQUIRE( size(A) == 3 );
		BOOST_REQUIRE( size(A) == 3 and size(A[0]) == 3 );
		BOOST_REQUIRE( A[1][1] == 34.4 );
		A = {
			{ 00.,  01., 02.},
			{ 10.,  11., 12.},
			{ 20.,  21., 22.}
		};
		BOOST_REQUIRE( A[1][2] == 12. );
	}
	{
		multi::array<double, 1> vec;
		vec = {4.0, 5.5};
		BOOST_REQUIRE( size(vec) == 2 );
		BOOST_REQUIRE( vec[1] == 5.5 );
	}
	{
		std::array<std::array<double, 2>, 3> const a = {
			{
				{ 1.2,  2.4},
				{11.2, 34.4},
				{15.2, 32.4}
			}
		};
		using std::begin; using std::end;
		multi::static_array<double, 2> A(begin(a), end(a));
		BOOST_REQUIRE( size(A) == 3 );
		BOOST_REQUIRE( size(A[0]) == 2 );
		BOOST_REQUIRE( A[1][0] == 11.2 );
	}
	{
		std::array<std::array<double, 2>, 3> const staticA = {
			{
				{ 1.2,  2.4},
				{11.2, 34.4},
				{15.2, 32.4}
			}
		};
		multi::static_array<double, 2> const A(std::begin(staticA), std::end(staticA));
		BOOST_REQUIRE(( A == multi::array<double, 2>{
				{ 1.2,  2.4},
				{11.2, 34.4},
				{15.2, 32.4}
			}
		));
		BOOST_REQUIRE(( 
			A == decltype(A){
				{ 1.2,  2.4},
				{11.2, 34.4},
				{15.2, 32.4}
			}
		));
	}
	{
//		multi::array<double, 2> A = 
//			(double const[][2]) // may warn with -Wpedantic
//			{
//				{ 1.2,  2.4},
//				{11.2, 34.4},
//				{15.2, 32.4}
//			}
//		;
//		BOOST_REQUIRE( size(A) == 3 );
//		BOOST_REQUIRE( A[1][0] == 11.2 );
	}
	{
		std::array<std::array<double, 2>, 3> a = {{
			{{1.,2.}},
			{{2.,4.}},
			{{3.,6.}}
		}};
		multi::array<double, 2> A(begin(a), end(a));
		BOOST_REQUIRE( num_elements(A) == 6 and A[2][1] == 6. );
	}
	{
		using complex = std::complex<double>; complex const I(0.,1.);
		multi::array<complex, 2> b = {
			{2. + 1.*I, 1. + 3.*I, 1. + 7.*I},
			{3. + 4.*I, 4. + 2.*I, 0. + 0.*I}
		};
		BOOST_REQUIRE( b[1][1] == 4. + 2.*I );
	}
}

BOOST_AUTO_TEST_CASE(multi_tests_static_array_initializer_list){
	multi::static_array<std::complex<double>, 2> SA = {
		{1. , 2.},
		{3. , 4.},
	};
	BOOST_REQUIRE( SA[1][1] == 4. );
}

BOOST_AUTO_TEST_CASE(multi_tests_initializer_list_3d){
	{
		multi::array<double, 3> const A = 
			{
				{
					{ 1.2, 0.}, 
					{ 2.4, 1.}
				},
				{
					{11.2,  3.}, 
					{34.4,  4.}
				},
				{
					{15.2, 99.}, 
					{32.4,  2.}
				}
			}
		;
		BOOST_REQUIRE( A[1][1][0] == 34.4 and A[1][1][1] == 4.   );
	}
}

BOOST_AUTO_TEST_CASE(multi_tests_initializer_list_3d_string){
	{
		using std::string;
		multi::array<string, 3> B3 = {
			{ {"000", "001", "002"}, 
			  {"010", "011", "012"} },
			{ {"100", "101", "102"}, 
			  {"110", "111", "112"} }
		};
		BOOST_REQUIRE( num_elements(B3)==12 and B3[1][0][1] == "101" );
	}

	#if defined(__cpp_deduction_guides)
	{	
		multi::array A({1., 2., 3.});
		static_assert( std::is_same<decltype(A)::element_type, double>{}, "!");
		BOOST_REQUIRE( size(A) == 3 and num_elements(A) == 3 );
		BOOST_REQUIRE( multi::rank<decltype(A)>{}==1 and num_elements(A)==3 and A[1]==2. ); 
		static_assert( typename decltype(A)::rank{}==1 );
	}
	{	
		multi::array A({1., 2.});
		static_assert( std::is_same<decltype(A)::element_type, double>{}, "!");
		BOOST_REQUIRE( size(A) == 2 and num_elements(A) == 2 );
		BOOST_REQUIRE( multi::rank<decltype(A)>{}==1 and num_elements(A)==2 and A[1]==2. ); BOOST_REQUIRE( multi::rank<decltype(A)>{}==1 );
	}
	{
		multi::array A({0, 2}); // 	multi::array A = {0, 2}; not working with CTAD
		static_assert( std::is_same<decltype(A)::element_type, int>{}, "!" );
		BOOST_REQUIRE( size(A) == 2 and num_elements(A) == 2 );
		BOOST_REQUIRE( multi::rank<decltype(A)>{}==1 and num_elements(A)==2 and A[1]==2. ); BOOST_REQUIRE( multi::rank<decltype(A)>{}==1 );
	}
	{
		multi::array A({9.}); // multi::array A = {9.}; not working with CTAD
		static_assert( std::is_same<decltype(A)::element_type, double>{}, "!" );
		BOOST_REQUIRE( multi::rank<decltype(A)>{}==1 and num_elements(A)==1 and A[0]==9. ); BOOST_REQUIRE( multi::rank<decltype(A)>{}==1 );
	}
	{
		multi::array A({9}); // multi::array A = {9}; not working with CTAD
		static_assert( std::is_same<decltype(A)::element_type, int>{}, "!" );
		BOOST_REQUIRE( size(A) == 1 and num_elements(A) == 1 );
		BOOST_REQUIRE( multi::rank<decltype(A)>{}==1 and num_elements(A)==1 and A[0]==9. ); BOOST_REQUIRE( multi::rank<decltype(A)>{}==1 );
	}
	{
		multi::array A({
			{1., 2., 3.}, 
			{4., 5., 6.}
		}); 
		BOOST_REQUIRE( multi::rank<decltype(A)>{}==2 and num_elements(A)==6 );
	}
	#endif
}

