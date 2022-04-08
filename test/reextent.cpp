// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
// Copyright 2018-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi reextent"
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(array_reextent) {
	multi::array<double, 2> A({2, 3});
	BOOST_REQUIRE( num_elements(A) == 6 );

	A[1][2] = 6.;
	BOOST_REQUIRE( A[1][2] == 6. );

	multi::array<double, 2> C({2, 3});
	BOOST_REQUIRE(size(C) == 2);
	BOOST_REQUIRE(size(C[0]) == 3);

	A.reextent({5, 4}, 99.);
	BOOST_REQUIRE( num_elements(A)== 5L*4L );
	BOOST_REQUIRE( A[1][2] ==  6. );  // reextent preserves values when it can...
	BOOST_REQUIRE( A[4][3] == 99. );  // ...and gives selected value to the rest
}

BOOST_AUTO_TEST_CASE(array_move_clear) {
	multi::array<double, 2> A({2, 3});
	A = multi::array<double, 2>(extensions(A), 123.);
	BOOST_REQUIRE( A[1][2] == 123. );

	A.clear();  // clear(A);
	BOOST_REQUIRE( num_elements(A) == 0 );
	BOOST_REQUIRE( size(A) == 0 );

	A.reextent({5, 4}, 66.);
	BOOST_REQUIRE( A[4][3] == 66. );
}

BOOST_AUTO_TEST_CASE(array_reextent_1d) {
	multi::array<double, 1> A(multi::extensions_t<1>{multi::iextension{10}}, 4.);
	BOOST_REQUIRE( size(A) == 10 );
	BOOST_REQUIRE( A[9] == 4. );

	A.reextent(multi::extensions_t<1>{multi::iextension{20}});
	BOOST_REQUIRE( size(A) == 20 );
	BOOST_REQUIRE( A[9] == 4. );
//	BOOST_REQUIRE( A[19] == 0. ); // impossible to know since it is sometimes 0.

//  A.reextent(std::tuple<int>(22) );
	A.reextent( boost::multi::tuple<int>(22) );
	BOOST_REQUIRE( size(A) == 22 );
	BOOST_REQUIRE( A[9] == 4. );

	A.reextent( {23} );
	BOOST_REQUIRE( size(A) == 23 );

#pragma warning(push)                // NOLINT(clang-diagnostic-unknown-pragmas)
#pragma warning (disable:1478 1786)  // NOLINT(clang-diagnostic-unknown-pragmas)
#pragma diagnostic push
#pragma diag_suppress 1215,1216,1444,1445
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
	A.reextent( std::make_tuple(24) );
	BOOST_REQUIRE( size(A) == 24 );
#pragma GCC diagnostic pop
#pragma diagnostic pop
#pragma warning(pop)                 // NOLINT(clang-diagnostic-unknown-pragmas)
}

inline void fff(boost::multi::detail::tuple<long> /*t*/) {}  // NOLINT(google-runtime-int) for testing

#pragma warning(push)                // NOLINT(clang-diagnostic-unknown-pragmas)
#pragma warning (disable:1478 1786)  // NOLINT(clang-diagnostic-unknown-pragmas)
#pragma diagnostic push
#pragma diag_suppress 1215,1216,1444,1445
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
BOOST_AUTO_TEST_CASE(tuple_implicit_test) {
	fff(1L);
	fff(1);
}

BOOST_AUTO_TEST_CASE(tuple_conversion_deprecated) {
	boost::multi::tuple<int, int> t{1, 1};
	BOOST_REQUIRE( t == std::make_tuple(1, 1) );
}
#pragma GCC diagnostic pop
#pragma diagnostic pop
#pragma warning(pop)                 // NOLINT(clang-diagnostic-unknown-pragmas)

BOOST_AUTO_TEST_CASE(array_reextent_0D) {
	multi::array<double, 0> A({}, 4.);
//	A.reextent(A.extensions()); // TODO(correaa) : fix unused for D = 0
	BOOST_REQUIRE( *A.data_elements() == 4. );
}

BOOST_AUTO_TEST_CASE(array_reextent_1d_with_initialization) {
	multi::array<double, 1> A(multi::extensions_t<1>{multi::iextension{10}}, 4.);
	BOOST_REQUIRE( size(A) == 10 );
	BOOST_REQUIRE( A[9] == 4. );

	A.reextent(multi::extensions_t<1>{multi::iextension{20}}, 8.);
	BOOST_REQUIRE( size(A) == 20 );
	BOOST_REQUIRE( A[9] == 4. );
	BOOST_REQUIRE( A[19] == 8. );
}

BOOST_AUTO_TEST_CASE(array_reextent_2d) {
	multi::array<double, 2> A({10, 20}, 4.);
	BOOST_REQUIRE( A[1][2] == 4. );

	A.clear();
	BOOST_REQUIRE( num_elements(A) == 0 );
	BOOST_REQUIRE( size(A) == 0 );

	A.reextent({20, 30}, 9.);
	BOOST_REQUIRE( A[1][2] = 9. );
	BOOST_REQUIRE( A[11][22] = 9. );
}

BOOST_AUTO_TEST_CASE(array_reextent_2d_array) {
	multi::array<double, 2> A({10, 20}, 4.);
	BOOST_REQUIRE( A[1][2] == 4. );

	A.clear();
	BOOST_REQUIRE( num_elements(A) == 0 );
	BOOST_REQUIRE( size(A) == 0 );
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
constexpr auto comp_equal(T t, U u) noexcept -> bool {
    using UT = std::make_unsigned_t<T>;
    using UU = std::make_unsigned_t<U>;
    if constexpr (std::is_signed_v<T> == std::is_signed_v<U>) {
		return t == u;
	} else if constexpr (std::is_signed_v<T>) {
		return t < 0 ? false : UT(t) == u;
	} else {
		return u < 0 ? false : t == UU(u);
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
	std::vector<double> v(100);
	{
	//  multi::array<double, 1> a(                           v.size() );  // warning: sign-conversion
		multi::array<double, 1> a(static_cast<multi::size_t>(v.size()));
		BOOST_REQUIRE( comp_equal(a.size(), v.size()) );
	}
	{
	//  multi::array<double, 1> a({v.size()});                              // warning: sign-conversion
	//	multi::array<double, 1> a({static_cast<multi::size_t>(v.size())});  // semantic problem, thinks it is an initializer_list
	//	BOOST_REQUIRE( comp_equal(a.size(), v.size()) );
	}
	{
	 	multi::array<double, 1> a(multi::iextensions<1>(static_cast<multi::size_t>(v.size())));  // warning: sign-conversion
	//	multi::array<double, 1> a(static_cast<multi::size_t>(v.size()));
		BOOST_REQUIRE( comp_equal(a.size(), v.size()) );
	}
}
