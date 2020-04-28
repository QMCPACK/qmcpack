#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $0 -o $0x -lboost_unit_test_framework&&$0x&&rm $0x;exit
#endif
// © Alfredo A. Correa 2018-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi element access"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

#if __cpp_lib_apply>=201603
#include<tuple> // apply
#else
#include<experimental/tuple>
#endif

namespace multi = boost::multi;

template<class T> void what(T&&) = delete;

BOOST_AUTO_TEST_CASE(multi_tests_element_access_with_tuple){
	multi::array<double, 2> m({3, 3}, 44.);
	std::array<int, 2> p = {1, 2};
	
	BOOST_REQUIRE( m[p[0]][p[1]] == m(1, 2) );
	BOOST_REQUIRE( &m(p[0], p[1]) == &m[p[0]][p[1]] );

	using std::
#if not(__cpp_lib_apply>=201603)
		experimental::
#endif
		apply;
	BOOST_REQUIRE( &m[p[0]][p[1]] == &apply(m, p) );

}

#if 1
BOOST_AUTO_TEST_CASE(multi_test_constness_reference){
	multi::array<double, 2> const m({10, 10}, 99.);

	BOOST_REQUIRE( size( m(1, {0, 3}) ) == 3 );

	BOOST_REQUIRE( m(1, {0, 3})[1] == 99. );
	BOOST_REQUIRE( m({0, 3}, 1).dimensionality == 1 );
	BOOST_TEST_REQUIRE( size(m.sliced(0, 3)) == 3 );

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

	BOOST_REQUIRE( size(m.range({0, 3}).rotated().paren(1l).unrotated()) == 3 );
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
#endif

