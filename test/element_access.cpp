// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
// © Alfredo A. Correa 2018-2021

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi element access"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>
//#include <boost/test/execution_monitor.hpp>  // for boost::execution_exception

#include "../array.hpp"

namespace multi = boost::multi;

template<class T> void what(T&&) = delete;

BOOST_AUTO_TEST_CASE(empty_intersection){
	multi::array<double, 1> arr(10);
	multi::array<double, 1> arr2;

	auto const is = intersection(arr.extension(), arr2.extension());
	BOOST_REQUIRE( arr(is).is_empty() );
	arr2(is) = arr(is);
	
	BOOST_REQUIRE( arr2(is) == arr(is) );
}

BOOST_AUTO_TEST_CASE(multi_tests_element_access_with_tuple){
	multi::array<double, 2> m({3, 3}, 44.);
	std::array<int, 2> p = {1, 2};
	
	BOOST_REQUIRE( m[p[0]][p[1]] == m(1, 2) );
	BOOST_REQUIRE( &m(p[0], p[1]) == &m[p[0]][p[1]] );

#if 0
	using std::
#if not defined(__cpp_lib_apply) or not(__cpp_lib_apply>=201603)
		experimental::
#endif
		apply;
	BOOST_REQUIRE( &m[p[0]][p[1]] == &apply(m, p) );
#endif

	BOOST_REQUIRE( &m[p[0]][p[1]] == &m(p[0], p[1]) );
	BOOST_REQUIRE( &m(p[0], p[1]) == &m.apply(p) );
}

BOOST_AUTO_TEST_CASE(multi_tests_extension_with_tuple){
	std::tuple<int, int> t = {3, 3};
	multi::array<double, 2> m1(t, 44.);
	BOOST_REQUIRE( size(m1) == 3 );
	
	std::array<int, 3> a = {3, 3};
	multi::array<double, 2> m2(a, 55.);
}

#if 1
BOOST_AUTO_TEST_CASE(multi_test_constness_reference){
	multi::array<double, 2> const m({10, 10}, 99.);

	BOOST_REQUIRE( size( m(1, {0, 3}) ) == 3 );

	BOOST_REQUIRE( m(1, {0, 3})[1] == 99. );
	BOOST_REQUIRE( m({0, 3}, 1).dimensionality == 1 );
	BOOST_REQUIRE( size(m.sliced(0, 3)) == 3 );

	BOOST_REQUIRE( m.range({0, 3}).rotated()[1].unrotated().size() == 3 );

//	BOOST_TEST_REQUIRE( size(m.paren(multi::index_range{0, 3}, 1)) == 3 );
//range(a).rotated().paren(as...).unrotated()
	BOOST_REQUIRE( m({0, 3}, {0, 3})[1][1] == 99. );

	static_assert(not std::is_assignable<decltype(m(1, {0, 3})[1]), double>{}, "!");
// none of these lines should compile because m is read-only
//	m(1, {0, 3})[1] = 88.;
//	m({0, 3}, 1)[1] = 77.;
//	m({0, 3}, {0, 3})[1][1] = 66.;
}

BOOST_AUTO_TEST_CASE(multi_test_non_constness_reference){
	multi::array<double, 2> m({10, 10}, 99.);

	BOOST_REQUIRE( size( m(1, {0, 3}) ) == 3 );
	static_assert(std::is_assignable<decltype(m(1, {0, 3})[1]), double>{}, "!");

	BOOST_REQUIRE( m(1, {0, 3})[1] == 99. );
	BOOST_REQUIRE( m({0, 3}, 1).dimensionality == 1 );
	BOOST_REQUIRE( size(m.sliced(0, 3)) == 3 );

	BOOST_REQUIRE( size(m.range({0, 3}).rotated().paren(1L).unrotated()) == 3 );
	BOOST_REQUIRE( size(m(multi::index_range{0, 3}, 1)) == 3 );

//	BOOST_REQUIRE( size(m({0, 3}, 1)) == 3 );
 
	BOOST_REQUIRE( m({0, 3}, {0, 3})[1][1] == 99. );

	m(1, {0, 3})[1] = 88.; 
	BOOST_REQUIRE( m(1, {0, 3})[1] == 88. );

//	m({0, 3}, 1)[1] = 77.; 
//	BOOST_REQUIRE( m({0, 3}, 1)[1] == 77. );

//	m({0, 3}, {0, 3})[1][1] = 66.; 
//	BOOST_REQUIRE( m({0, 3}, 1)[1] == 66. );
}

BOOST_AUTO_TEST_CASE(multi_test_stencil){

	multi::array<std::string, 2> A = 
		{{"a", "b", "c", "d", "e"},
		 {"f", "g", "h", "f", "g"},
		 {"h", "i", "j", "k", "l"}}
	;
	
	BOOST_REQUIRE(      size(A) == 3                                           );
	BOOST_REQUIRE(           A.num_elements() == 3*5                           );
	BOOST_REQUIRE(           A[1][2] == "h"                                    );

	BOOST_REQUIRE(      size(A          ({1, 3}, {2, 5})) == 2                 );
	BOOST_REQUIRE( extension(A          ({1, 3}, {2, 5})).start() == 0         );
	BOOST_REQUIRE(           A          ({1, 3}, {2, 5}).num_elements() == 2*3 );
	BOOST_REQUIRE(           A          ({1, 3}, {2, 5}).num_elements() == 2*3 );
	BOOST_REQUIRE(           A          ({1, 3}, {2, 5})[0][0] == "h"          );
	BOOST_REQUIRE(          &A          ({1, 3}, {2, 5})[0][0] == &A[1][2]     );

	BOOST_REQUIRE(      size(A.stenciled({1, 3}, {2, 5})) == 2                 );
	BOOST_REQUIRE( extension(A.stenciled({1, 3}, {2, 5})).start() == 1         );
	BOOST_REQUIRE(           A.stenciled({1, 3}, {2, 5}).num_elements() == 2*3 );
	BOOST_REQUIRE(           A.stenciled({1, 3}, {2, 5}) [1][2] == "h"         ); 
	BOOST_REQUIRE(          &A.stenciled({1, 3}, {2, 5}) [1][2] == &A[1][2]    );

}

BOOST_AUTO_TEST_CASE(multi_extension_intersection){
	multi::array<double, 1> A = {{2., 2., 2.}};
	multi::array<double, 1> B = {{3., 3., 3., 3.}};

	BOOST_REQUIRE( intersection( extension(A), extension(B) ).size() == 3 );
	for(auto i : intersection( extension(A), extension(B) ) ){
		B[i] += A[i];
	}

	BOOST_REQUIRE( B[2] == 5. );
	BOOST_REQUIRE( B[3] == 3. );
}

BOOST_AUTO_TEST_CASE(multi_extensions_intersection){
	multi::array<double, 2> A = {{2., 2., 2.}, {2., 2., 2.}};
	multi::array<double, 2> B = {{3., 3., 3., 3.}, {3., 3., 3., 3.}, {3., 3., 3., 3.}} ;

	BOOST_REQUIRE( extensions(A).num_elements() ==  6 );
	BOOST_REQUIRE( extensions(B).num_elements() == 12 );
	
	auto is = intersection( extensions(A), extensions(B) );
	BOOST_REQUIRE( is.num_elements() == 6 );

	multi::array<double, 2> C(is);
	C(std::get<0>(is), std::get<1>(is)) = A;
	BOOST_REQUIRE( C == A );
}

BOOST_AUTO_TEST_CASE(multi_extensions_intersection_2){
	multi::array<double, 2> A({80, 20}, 4.);
	multi::array<double, 2> B({30, 70}, 8.);

	BOOST_REQUIRE( extensions(A).num_elements() == 80*20 );
	BOOST_REQUIRE( extensions(B).num_elements() == 30*70 );
	
	auto is = intersection( extensions(A), extensions(B) );
	BOOST_REQUIRE( is.num_elements() == 30*20 );

	multi::array<double, 2> C(is);
	C(std::get<0>(is), std::get<1>(is)) = A(std::get<0>(is), std::get<1>(is));
	BOOST_REQUIRE( C[16][17] == 4. );
	
	C(std::get<0>(is), std::get<1>(is)) = B(std::get<0>(is), std::get<1>(is));
	BOOST_REQUIRE( C[16][17] == 8. );
}


#endif
