// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2023 Alfredo A. Correa

// #define BOOST_TEST_MODULE "C++ Unit Tests for Multi reextent"  // test title NOLINT(cppcoreguidelines-macro-usage)
#include <boost/test/unit_test.hpp>

#include "multi/array.hpp"

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(array_reextent) {
	multi::array<double, 2> arr({2, 3});
	BOOST_REQUIRE( num_elements(arr) == 6 );

	arr[1][2] = 6.0;
	BOOST_REQUIRE( arr[1][2] == 6.0 );

	multi::array<double, 2> arr3({2, 3});
	BOOST_REQUIRE(size(arr3) == 2);
	BOOST_REQUIRE(size(arr3[0]) == 3);

	arr.reextent({5, 4}, 99.);
	BOOST_REQUIRE( num_elements(arr)== 5L*4L );
	BOOST_REQUIRE( arr[1][2] ==  6. );  // reextent preserves values when it can...
	BOOST_REQUIRE( arr[4][3] == 99. );  // ...and gives selected value to the rest
}

BOOST_AUTO_TEST_CASE(array_reextent_noop) {
	multi::array<double, 2> arr({2, 3});
	BOOST_REQUIRE( num_elements(arr) == 6 );

	arr[1][2] = 6.;
	BOOST_REQUIRE( arr[1][2] == 6. );

	multi::array<double, 2> arr3({2, 3});
	BOOST_REQUIRE(size(arr3) == 2);
	BOOST_REQUIRE(size(arr3[0]) == 3);

	auto* const A_base = arr.base();
	arr.reextent({2, 3});
	BOOST_REQUIRE( num_elements(arr)== 2L*3L );
	BOOST_REQUIRE( arr[1][2] ==  6.0 );  // reextent preserves values when it can...

	BOOST_REQUIRE( A_base == arr.base() );
}

BOOST_AUTO_TEST_CASE(array_reextent_noop_with_init) {
	multi::array<double, 2> arr({2, 3});
	BOOST_REQUIRE( num_elements(arr) == 6 );

	arr[1][2] = 6.0;
	BOOST_REQUIRE( arr[1][2] == 6.0 );

	multi::array<double, 2> arr3({2, 3});
	BOOST_REQUIRE(size(arr3) == 2);
	BOOST_REQUIRE(size(arr3[0]) == 3);

	auto* const A_base = arr.base();
	arr.reextent({2, 3}, 99.0);
	BOOST_REQUIRE( num_elements(arr)== 2L*3L );
	BOOST_REQUIRE( arr[1][2] ==  6.0 );  // reextent preserves values when it can...

	BOOST_REQUIRE( A_base == arr.base() );
}

BOOST_AUTO_TEST_CASE(array_reextent_moved) {
	multi::array<double, 2> arr({2, 3});
	BOOST_REQUIRE( num_elements(arr) == 6 );

	arr[1][2] = 6.;
	BOOST_REQUIRE( arr[1][2] == 6.0 );

	auto* const A_base = arr.base();

	arr = std::move(arr).reextent({2, 3});  // "arr = ..." suppresses linter bugprone-use-after-move,hicpp-invalid-access-moved

	BOOST_TEST_REQUIRE( arr.size() == 2 );
	BOOST_REQUIRE( arr.num_elements() == 2L*3L );
	BOOST_REQUIRE( num_elements(arr)== 2L*3L );
	BOOST_TEST(arr[1][2] == 6.0);  // after move the original elments might not be the same

	BOOST_REQUIRE( A_base == arr.base() );
}

BOOST_AUTO_TEST_CASE(array_reextent_moved_trivial) {
	multi::array<double, 2> arr({2, 3});
	BOOST_REQUIRE( num_elements(arr) == 6 );

	arr[1][2] = 6.0;
	BOOST_REQUIRE( arr[1][2] == 6.0 );

	auto* const A_base = arr.base();
	arr                = std::move(arr).reextent({2, 3});  // "arr = ..." suppresses linter bugprone-use-after-move,hicpp-invalid-access-moved
	BOOST_REQUIRE( num_elements(arr)== 2L*3L );
	BOOST_REQUIRE( arr[1][2] ==  6.0 );  // after move the original elments might not be the same

	BOOST_REQUIRE( A_base == arr.base() );
}

BOOST_AUTO_TEST_CASE(array_reextent_moved_trivial_change_extents) {
	multi::array<double, 2> arr({2, 3});
	BOOST_REQUIRE( num_elements(arr) == 6 );

	arr[1][2] = 6.0;
	BOOST_REQUIRE( arr[1][2] == 6.0 );

	auto* const A_base = arr.base();
	arr                = std::move(arr).reextent({4, 5});
	BOOST_REQUIRE( num_elements(arr)== 4L*5L );
	// BOOST_REQUIRE( arr[1][2] !=  6.0 );  // after move the original elements might not be the same, but it is not 100% possible to check

	BOOST_REQUIRE( A_base != arr.base() );
}

BOOST_AUTO_TEST_CASE(array_move_clear) {
	multi::array<int, 1>    const iarr;
	multi::array<double, 2> arr({2, 3});

	arr = multi::array<double, 2>(extensions(arr), 123.0);
	BOOST_REQUIRE( arr[1][2] == 123.0 );

	arr.clear();  // clear(arr);
	BOOST_REQUIRE( num_elements(arr) == 0 );
	BOOST_REQUIRE( size(arr) == 0 );

	arr.reextent({5, 4}, 66.0);
	BOOST_REQUIRE( arr[4][3] == 66.0 );

	BOOST_REQUIRE(iarr.is_empty());
}

BOOST_AUTO_TEST_CASE(array_reextent_1d) {
	multi::array<double, 1> arr(multi::extensions_t<1>{multi::iextension{10}}, 4.0);
	BOOST_REQUIRE( size(arr) == 10 );
	BOOST_REQUIRE( arr[9] == 4.0 );

	arr.reextent(multi::extensions_t<1>{multi::iextension{20}});
	BOOST_REQUIRE( size(arr) == 20 );
	BOOST_REQUIRE( arr[9] == 4.0 );
	// BOOST_REQUIRE( arr[19] == 0.0 );  // impossible to know since it is only sometimes 0.0

	arr.reextent(boost::multi::tuple<int>(22));
	BOOST_REQUIRE( size(arr) == 22 );
	BOOST_REQUIRE( arr[9] == 4.0 );

	arr.reextent({23});
	BOOST_REQUIRE( size(arr) == 23 );
}

BOOST_AUTO_TEST_CASE(tuple_decomposition) {
	boost::multi::tuple<int, int> const tup{1, 2};
	auto [t0, t1] = tup;
	BOOST_REQUIRE( t0 == 1 );
	BOOST_REQUIRE( t1 == 2 );
}

BOOST_AUTO_TEST_CASE(array_reextent_0D) {
	multi::array<double, 0> const arr({}, 4.0);
	//  arr.reextent(arr.extensions());  // TODO(correaa) : fix unused for D = 0
	BOOST_REQUIRE( *arr.data_elements() == 4.0 );
}

BOOST_AUTO_TEST_CASE(array_reextent_1d_with_initialization) {
	multi::array<double, 1> arr(multi::extensions_t<1>{multi::iextension{10}}, 4.0);
	BOOST_REQUIRE( size(arr) == 10 );
	BOOST_REQUIRE( arr[9] == 4.0 );

	arr.reextent(multi::extensions_t<1>{multi::iextension{20}}, 8.0);
	BOOST_REQUIRE( size(arr) == 20 );
	BOOST_REQUIRE( arr[9] == 4.0 );
	BOOST_REQUIRE( arr[19] == 8.0 );
}

BOOST_AUTO_TEST_CASE(array_reextent_2d) {
	multi::array<double, 2> arr({10, 20}, 4.0);
	BOOST_REQUIRE( arr[1][2] == 4.0 );

	arr.clear();
	BOOST_REQUIRE( num_elements(arr) == 0 );
	BOOST_REQUIRE( size(arr) == 0 );

	arr.reextent({20, 30}, 9.0);
	BOOST_REQUIRE( arr[1][2] = 9. );
	BOOST_REQUIRE( arr[11][22] = 9.0 );
}

BOOST_AUTO_TEST_CASE(array_reextent_2d_with_move) {
	multi::array<int, 2> arr = {
		{1, 2, 3},
		{4, 5, 6},
	};
	BOOST_REQUIRE( arr.size() == 2 );

	arr = std::move(arr).reextent({3, 2});

	BOOST_REQUIRE( arr.size() == 3 );
	BOOST_REQUIRE( arr[1][2] = 10 );
}

BOOST_AUTO_TEST_CASE(array_reextent_2d_array) {
	multi::array<double, 2> arr({10, 20}, 4.0);
	BOOST_REQUIRE( arr[1][2] == 4.0 );

	arr.clear();
	BOOST_REQUIRE( num_elements(arr) == 0 );
	BOOST_REQUIRE( size(arr) == 0 );
}

template<class T, class U>
constexpr auto comp_equal(T left, U right) noexcept -> bool {
	using UT = std::make_unsigned_t<T>;
	using UU = std::make_unsigned_t<U>;
	if constexpr(std::is_signed_v<T> == std::is_signed_v<U>) {
		return left == right;
	} else if constexpr(std::is_signed_v<T>) {
		return left < 0 ? false : static_cast<UT>(left) == right;
	} else {
		return right < 0 ? false : left == UU(right);
	}
#if not defined(__INTEL_COMPILER) and not defined(__NVCOMPILER)
	__builtin_unreachable();
#endif
}

BOOST_AUTO_TEST_CASE(array_vector_size) {
	std::vector<double> const vec(100);  // std::vector NOLINT(fuchsia-default-arguments-calls)
	{
		// multi::array<double, 1> a(                             vec.size() );  // warning: sign-conversion
		multi::array<double, 1> const arr(static_cast<multi::size_t>(vec.size()));
		BOOST_REQUIRE( comp_equal(arr.size(), vec.size()) );
	}
	{
		multi::array<double, 1> const arr(multi::iextensions<1>(static_cast<multi::size_t>(vec.size())));  // warning: sign-conversion
		// multi::array<double, 1> a(static_cast<multi::size_t>(v.size()));
		BOOST_REQUIRE( comp_equal(arr.size(), vec.size()) );
	}
}

BOOST_AUTO_TEST_CASE(array_iota) {
	multi::array<double, 1> const Aarr(10);
	multi::array<int, 1>    Barr(Aarr.extension().begin(), Aarr.extension().end());
	BOOST_REQUIRE( Barr[0] == 0 );
	BOOST_REQUIRE( Barr[1] == 1 );
	BOOST_REQUIRE( Barr[9] == 9 );

	multi::array<int, 1> Carr(Aarr.extension());
	BOOST_REQUIRE( Carr[0] == 0 );
	BOOST_REQUIRE( Carr[1] == 1 );
	BOOST_REQUIRE( Carr[9] == 9 );

	multi::array<int, 1> const Darr(Aarr.extensions());
	BOOST_REQUIRE( Darr.extensions() == Aarr.extensions() );
}

#ifndef __INTEL_COMPILER
BOOST_AUTO_TEST_CASE(extension_index_op) {
	multi::array<double, 2> const Aarr({11, 13});
	auto Aext = Aarr.extensions();
	BOOST_REQUIRE( std::get<0>(Aext[3][5]) == 3 );
	BOOST_REQUIRE( std::get<1>(Aext[3][5]) == 5 );

	for(int i = 0; i != 3; ++i) {
		for(int j = 0; j != 5; ++j) {
			auto [ip, jp] = Aext[i][j];
			BOOST_REQUIRE(ip == i);
			BOOST_REQUIRE(jp == j);
		}
	}
}
#endif
