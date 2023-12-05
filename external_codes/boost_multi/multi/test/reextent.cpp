// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi reextent"
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(array_reextent) {
	multi::array<double, 2> arr({2, 3});
	BOOST_REQUIRE( num_elements(arr) == 6 );

	arr[1][2] = 6.;
	BOOST_REQUIRE( arr[1][2] == 6. );

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
	BOOST_REQUIRE( arr[1][2] ==  6. );  // reextent preserves values when it can...

	BOOST_REQUIRE( A_base == arr.base() );
}

BOOST_AUTO_TEST_CASE(array_reextent_noop_with_init) {
	multi::array<double, 2> arr({2, 3});
	BOOST_REQUIRE( num_elements(arr) == 6 );

	arr[1][2] = 6.;
	BOOST_REQUIRE( arr[1][2] == 6. );

	multi::array<double, 2> arr3({2, 3});
	BOOST_REQUIRE(size(arr3) == 2);
	BOOST_REQUIRE(size(arr3[0]) == 3);

	auto* const A_base = arr.base();
	arr.reextent({2, 3}, 99.);
	BOOST_REQUIRE( num_elements(arr)== 2L*3L );
	BOOST_REQUIRE( arr[1][2] ==  6. );  // reextent preserves values when it can...

	BOOST_REQUIRE( A_base == arr.base() );
}

BOOST_AUTO_TEST_CASE(array_reextent_moved) {
	multi::array<double, 2> arr({2, 3});
	BOOST_REQUIRE( num_elements(arr) == 6 );

	arr[1][2] = 6.;
	BOOST_REQUIRE( arr[1][2] == 6. );

	auto* const A_base = arr.base();
	arr = std::move(arr).reextent({2, 3}, 99.);  // "arr = ..." suppresses linter bugprone-use-after-move,hicpp-invalid-access-moved
	BOOST_REQUIRE( num_elements(arr)== 2L*3L );
	BOOST_TEST( arr[1][2] ==  6. );  // after move the original elments might not be the same

	BOOST_REQUIRE( A_base == arr.base() );
}

BOOST_AUTO_TEST_CASE(array_reextent_moved_trivial) {
	multi::array<double, 2> arr({2, 3});
	BOOST_REQUIRE( num_elements(arr) == 6 );

	arr[1][2] = 6.;
	BOOST_REQUIRE( arr[1][2] == 6. );

	auto* const A_base = arr.base();
	arr = std::move(arr).reextent({2, 3});  // "arr = ..." suppresses linter bugprone-use-after-move,hicpp-invalid-access-moved
	BOOST_REQUIRE( num_elements(arr)== 2L*3L );
	BOOST_REQUIRE( arr[1][2] ==  6. );  // after move the original elments might not be the same

	BOOST_REQUIRE( A_base == arr.base() );
}

BOOST_AUTO_TEST_CASE(array_reextent_moved_trivial_change_extents) {
	multi::array<double, 2> arr({2, 3});
	BOOST_REQUIRE( num_elements(arr) == 6 );

	arr[1][2] = 6.;
	BOOST_REQUIRE( arr[1][2] == 6. );

	auto* const A_base = arr.base();
	arr = std::move(arr).reextent({4, 5});
	BOOST_REQUIRE( num_elements(arr)== 4L*5L );
	BOOST_REQUIRE( arr[1][2] !=  6. );  // after move the original elments might not be the same

	BOOST_REQUIRE( A_base != arr.base() );
}

BOOST_AUTO_TEST_CASE(array_move_clear) {
	multi::array<double, 2> arr({2, 3});
	arr = multi::array<double, 2>(extensions(arr), 123.);
	BOOST_REQUIRE( arr[1][2] == 123. );

	arr.clear();  // clear(arr);
	BOOST_REQUIRE( num_elements(arr) == 0 );
	BOOST_REQUIRE( size(arr) == 0 );

	arr.reextent({5, 4}, 66.);
	BOOST_REQUIRE( arr[4][3] == 66. );
}

BOOST_AUTO_TEST_CASE(array_reextent_1d) {
	multi::array<double, 1> arr(multi::extensions_t<1>{multi::iextension{10}}, 4.);
	BOOST_REQUIRE( size(arr) == 10 );
	BOOST_REQUIRE( arr[9] == 4. );

	arr.reextent(multi::extensions_t<1>{multi::iextension{20}});
	BOOST_REQUIRE( size(arr) == 20 );
	BOOST_REQUIRE( arr[9] == 4. );
//  BOOST_REQUIRE( arr[19] == 0. );  // impossible to know since it is sometimes 0.

	arr.reextent( boost::multi::tuple<int>(22) );
	BOOST_REQUIRE( size(arr) == 22 );
	BOOST_REQUIRE( arr[9] == 4. );

	arr.reextent( {23} );
	BOOST_REQUIRE( size(arr) == 23 );

#pragma warning(push)                      // NOLINT(clang-diagnostic-unknown-pragmas)
#pragma warning (disable:1478 1786)        // NOLINT(clang-diagnostic-unknown-pragmas)
#pragma nv_diagnostic push                    // NOLINT(clang-diagnostic-unknown-pragmas)
#pragma nv_diag_suppress 1215,1216,1444,1445  // NOLINT(clang-diagnostic-unknown-pragmas)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
//	arr.reextent( std::make_tuple(24) );
//	BOOST_REQUIRE( size(arr) == 24 );
#pragma GCC diagnostic pop
#pragma nv_diagnostic pop                     // NOLINT(clang-diagnostic-unknown-pragmas)
#pragma warning(pop)                       // NOLINT(clang-diagnostic-unknown-pragmas)
}

//inline void fff(boost::multi::detail::tuple<long> /*t*/) {}  // NOLINT(google-runtime-int) for testing

//#pragma warning(push)                     // NOLINT(clang-diagnostic-unknown-pragmas)
//#pragma warning (disable:1478 1786)       // NOLINT(clang-diagnostic-unknown-pragmas)
//#pragma diagnostic push                   // NOLINT(clang-diagnostic-unknown-pragmas)
//#pragma diag_suppress 1215,1216,1444,1445 // NOLINT(clang-diagnostic-unknown-pragmas)
//#pragma GCC diagnostic push
//#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
//BOOST_AUTO_TEST_CASE(tuple_implicit_test) {
//	fff(1L);
//	fff(1);
//}

//BOOST_AUTO_TEST_CASE(tuple_conversion_deprecated) {
//	boost::multi::tuple<int, int> t{1, 2};
//	BOOST_REQUIRE( std::get<0>(t) == 1 );
//	BOOST_REQUIRE( std::get<1>(t) == 2 );
//}
//#pragma GCC diagnostic pop
//#pragma diagnostic pop                   // NOLINT(clang-diagnostic-unknown-pragmas)
//#pragma warning(pop)                     // NOLINT(clang-diagnostic-unknown-pragmas)

BOOST_AUTO_TEST_CASE(tuple_decomposition) {
	boost::multi::tuple<int, int> tup{1, 2};
	auto [t0, t1] = tup;
	BOOST_REQUIRE( t0 == 1 );
	BOOST_REQUIRE( t1 == 2 );
}

BOOST_AUTO_TEST_CASE(array_reextent_0D) {
	multi::array<double, 0> arr({}, 4.);
//	arr.reextent(arr.extensions()); // TODO(correaa) : fix unused for D = 0
	BOOST_REQUIRE( *arr.data_elements() == 4. );
}

BOOST_AUTO_TEST_CASE(array_reextent_1d_with_initialization) {
	multi::array<double, 1> arr(multi::extensions_t<1>{multi::iextension{10}}, 4.);
	BOOST_REQUIRE( size(arr) == 10 );
	BOOST_REQUIRE( arr[9] == 4. );

	arr.reextent(multi::extensions_t<1>{multi::iextension{20}}, 8.);
	BOOST_REQUIRE( size(arr) == 20 );
	BOOST_REQUIRE( arr[9] == 4. );
	BOOST_REQUIRE( arr[19] == 8. );
}

BOOST_AUTO_TEST_CASE(array_reextent_2d) {
	multi::array<double, 2> arr({10, 20}, 4.);
	BOOST_REQUIRE( arr[1][2] == 4. );

	arr.clear();
	BOOST_REQUIRE( num_elements(arr) == 0 );
	BOOST_REQUIRE( size(arr) == 0 );

	arr.reextent({20, 30}, 9.);
	BOOST_REQUIRE( arr[1][2] = 9. );
	BOOST_REQUIRE( arr[11][22] = 9. );
}

BOOST_AUTO_TEST_CASE(array_reextent_2d_array) {
	multi::array<double, 2> arr({10, 20}, 4.);
	BOOST_REQUIRE( arr[1][2] == 4. );

	arr.clear();
	BOOST_REQUIRE( num_elements(arr) == 0 );
	BOOST_REQUIRE( size(arr) == 0 );
}

#pragma GCC diagnostic push
#pragma GCC diagnostic warning "-Wunknown-pragmas"
#if defined __NVCC__
	#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
		#pragma nv_diagnostic push
		#pragma nv_diag_suppress = implicit_return_from_non_void_function
	#else
		#pragma    diagnostic push
		#pragma    diag_suppress = implicit_return_from_non_void_function
	#endif
#elif defined __NVCOMPILER
	#pragma    diagnostic push
	#pragma    diag_suppress = implicit_return_from_non_void_function
#endif
template< class T, class U >
constexpr auto comp_equal(T left, U right) noexcept -> bool {
    using UT = std::make_unsigned_t<T>;
    using UU = std::make_unsigned_t<U>;
    if constexpr (std::is_signed_v<T> == std::is_signed_v<U>) {
		return left == right;
	} else if constexpr (std::is_signed_v<T>) {
		return left < 0 ? false : static_cast<UT>(left) == right;
	} else {
		return right < 0 ? false : left == UU(right);
	}
	#if not defined(__INTEL_COMPILER) and not defined(__NVCOMPILER)
	__builtin_unreachable();
	#endif
}
#if defined __NVCC__
	#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
		#pragma nv_diagnostic pop
	#else
		#pragma    diagnostic pop
	#endif
#elif defined __NVCOMPILER
	#pragma    diagnostic pop
#endif
#pragma GCC diagnostic pop

BOOST_AUTO_TEST_CASE(array_vector_size) {
	std::vector<double> vec(100);
	{
	//  multi::array<double, 1> a(                             vec.size() );  // warning: sign-conversion
		multi::array<double, 1> arr(static_cast<multi::size_t>(vec.size()));
		BOOST_REQUIRE( comp_equal(arr.size(), vec.size()) );
	}
	{
	 	multi::array<double, 1> arr(multi::iextensions<1>(static_cast<multi::size_t>(vec.size())));  // warning: sign-conversion
	//	multi::array<double, 1> a(static_cast<multi::size_t>(v.size()));
		BOOST_REQUIRE( comp_equal(arr.size(), vec.size()) );
	}
}
