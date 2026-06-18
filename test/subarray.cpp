// Copyright 2019-2025 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for implicit_cast, explicit_cast

#include <boost/core/lightweight_test.hpp>

#include <type_traits>  // for std::is_swappable_v
#include <utility>      // for as_const

namespace multi = boost::multi;

namespace {
auto f_arr(multi::array<int, 1> arr) {
	return arr[2];
}

// auto f_sub(multi::const_subarray<int, 1> const& arr) {
// 	return arr[2];
// }
}  // namespace

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	{
		auto f_arr_ret = f_arr({1, 2, 3});
		BOOST_TEST(f_arr_ret == 3 );

		// auto f_sub_ret = f_sub({1, 2, 3});
		// BOOST_TEST(f_sub_ret == 3 );
	}

	/* subarray_assignment */
	{
		multi::array<int, 3> A1({3, 4, 5}, 99);
		A1[2][1][1] = 88;

		auto constA2 = std::as_const(A1)[2];
		BOOST_TEST( constA2[1][1] == 88 );

		auto A2 = A1[2];
		BOOST_TEST( A2[1][1] == 88 );

		A2[1][1] = 77;
		BOOST_TEST( A2[1][1] == 77 );  // cppcheck-suppress knownConditionTrueFalse;
	}

	/* subarray_assignment */
	{
		multi::array<int, 3> A1 = {
			{{1, 2},
			 {3, 4}},
			{{5, 6},
			 {7, 8}},
		};

		auto const& R0 = std::as_const(A1)[0];
		auto&&      R1 = A1[1];

		R1 = R0;

		BOOST_TEST( A1[0] == A1[1] );
	}

	/* subarray_assignment */
	{
		multi::array<int, 3> A1 = {
			{{1, 2},
			 {3, 4}},
			{{5, 6},
			 {7, 8}},
		};

		auto const& R0 = A1[0];
		auto&&      R1 = A1[1];

		R1 = R0;

		BOOST_TEST( A1[0] == A1[1] );
	}

	/* subarray_base */
	{
		multi::array<int, 3> A1({3, 4, 5}, 99);

		auto&& Asub  = A1();
		*Asub.base() = 88;

		BOOST_TEST( A1[0][0][0] == 88 );

		*A1().base() = 77;

		BOOST_TEST( A1[0][0][0] == 77 );

		// *std::as_const(Asub).base() = 66;  // should not compile, read-only
	}

	/* test ref(begin, end)*/
	{
		multi::array<int, 2> A2D = {
			{1, 2},
			{3, 4},
		};
		BOOST_TEST( A2D[0][0] == 1 );

		// multi::const_subarray<int, 2> R2D(A2D.begin(), A2D.end());
		// BOOST_TEST( R2D.addressof()== A2D.addressof() );
		// BOOST_TEST( &R2D == &A2D );  // TODO(correaa) enable comparison with pointer value (array*)
	}

	// equality 1D
	{
		multi::array<int, 1> const AA = {1, 2, 3};
		multi::array<int, 1> const BB = {2, 3, 4};

		BOOST_TEST(   AA   != BB    );
		BOOST_TEST( !(AA   == BB  ) );

		BOOST_TEST(   AA() != BB()  );
		BOOST_TEST( !(AA() == BB()) );

#if defined(__cpp_multidimensional_subscript) && (__cpp_multidimensional_subscript >= 202110L)
		BOOST_TEST(   AA[] != BB[]  );
#endif
	}

	// equality 2D
	{
		multi::array<int, 2> const AA = {
			{1, 2},
			{3, 4},
		};
		multi::array<int, 2> const BB = {
			{2, 3},
			{4, 5},
		};

		BOOST_TEST(   AA   != BB    );
		BOOST_TEST( !(AA   == BB  ) );

		BOOST_TEST(   AA() != BB()  );
		BOOST_TEST( !(AA() == BB()) );
	}

	// equality 1D
	{
		multi::array<int, 1> const      AA = {1, 2, 3};
		multi::array<unsigned, 1> const BB = {2, 3, 4};

#if __cplusplus >= 202002L
		BOOST_TEST( std::cmp_not_equal(AA[0], BB[0]) );
#else
		BOOST_TEST( AA[0] != static_cast<int>(BB[0]) );
#endif

		auto const to_int = [](auto elem) noexcept {
			return static_cast<int>(elem);
		};

		BOOST_TEST(   AA   != BB.element_transformed(to_int) );

#if !defined(_MSC_VER)  // MSVC would warn deeply in the standard library here, TODO(correaa) make warning evident at top level lib
		BOOST_TEST(   AA   != BB    );
		BOOST_TEST( !(AA   == BB  ) );

		BOOST_TEST(   AA() != BB()  );
		BOOST_TEST( !(AA() == BB()) );
#endif
	}

	// equality 2D
	{
		multi::array<int, 2> const AA = {
			{1, 2},
			{3, 4},
		};
		multi::array<unsigned, 2> const BB = {
			{2, 3},
			{4, 5},
		};

		auto const to_int = [](auto const& elem) noexcept {
			return static_cast<int>(elem);
		};

		BOOST_TEST(   AA   != BB.element_transformed(to_int)  );
		BOOST_TEST( !(AA   == BB.element_transformed(to_int)) );

#if !defined(_MSC_VER)  // MSVC would warn deeply in the standard library here, TODO(correaa) make warning evident at top level lib
		BOOST_TEST(   AA   != BB    );
		BOOST_TEST( !(AA   == BB  ) );

		BOOST_TEST(   AA() != BB()  );
		BOOST_TEST( !(AA() == BB()) );
#endif
	}

	/* test ref(begin, end)*/
	{
		multi::array<int, 2> A2D = {
			{1, 2},
			{3, 4},
		};
		BOOST_TEST( A2D[0][0] == 1 );

		multi::subarray<int, 2> R2D(A2D.begin(), A2D.end());
		// BOOST_TEST( R2D.addressof()== A2D.addressof() );
		// BOOST_TEST( &R2D == &A2D );  // TODO(correaa) allow address comparison with value-pointer (array*)

		R2D[0][0] = 77;
		BOOST_TEST( R2D[0][0] == 77 );  // cppcheck-suppress knownConditionTrueFalse;
	}

	{  // https://godbolt.org/z/a5dr7YvMz
		namespace multi = boost::multi;

		// template<class T, ::boost::multi::dimensionality_t D>
		// struct std::is_trivially_relocatable<::boost::multi::array<T, D>> : std::true_type{};

		// template<class T, ::boost::multi::dimensionality_t D>
		// struct std::is_trivially_relocatable<::boost::multi::array_ref<T, D>> : std::true_type{};

		// template<class T, ::boost::multi::dimensionality_t D>
		// struct std::is_trivially_relocatable<::boost::multi::subarray<T, D>> : std::true_type{};

		static_assert(std::is_nothrow_move_constructible_v<multi::array<int, 4>>);
		static_assert(!std::is_trivially_copyable_v<multi::array<int, 4>>);
		static_assert(std::is_swappable_v<multi::array<int, 4>>);
		static_assert(std::is_nothrow_swappable_v<multi::array<int, 4>>);
		// static_assert( std::is_trivially_relocatable_v<multi::array<int, 4>>);  // <==========

		// static_assert(!std::is_move_constructible_v<multi::array_ref<int, 4>>);
		// static_assert(!std::is_nothrow_move_constructible_v<multi::array_ref<int, 4>>);

#if !defined(__NVCOMPILER) || (__NVCOMPILER_MAJOR__ >= 24) && !defined(__NVCC__)
		static_assert(!std::is_copy_constructible_v<multi::array_ref<int, 4>>);
#endif

#ifndef __circle_build__
		static_assert(!std::is_trivially_copyable_v<multi::array_ref<int, 4>>);
#endif
		static_assert(std::is_copy_assignable_v<multi::array_ref<int, 4>>);
		static_assert(!std::is_trivially_copy_assignable_v<multi::array_ref<int, 4>>);
		static_assert(std::is_swappable_v<multi::array_ref<int, 4>>);
		// static_assert( std::is_trivially_relocatable_v     <multi::array_ref<int, 4>>);  // <==========

		static_assert(std::is_move_constructible_v<multi::subarray<int, 4>>);          // mmm, something strange here
		static_assert(std::is_nothrow_move_constructible_v<multi::subarray<int, 4>>);  // mmm, something strange here

#ifndef __NVCC__
		static_assert(!std::is_copy_constructible_v<multi::subarray<int, 4>>);
#endif

#ifndef __circle_build__
		static_assert(!std::is_trivially_copyable_v<multi::subarray<int, 4>>);
#endif

		static_assert(std::is_copy_assignable_v<multi::subarray<int, 4>>);
		static_assert(!std::is_trivially_copy_assignable_v<multi::subarray<int, 4>>);
		static_assert(std::is_swappable_v<multi::subarray<int, 4>>);  // TODO(correaa) fix?

		// static_assert(    std::is_nothrow_swappable_v      <multi::subarray<int, 4>>);
		// static_assert( std::is_trivially_relocatable_v     <multi::subarray<int, 4>>);  // <==========

		static_assert(std::is_move_constructible_v<multi::array<int, 4>::iterator>);
		static_assert(std::is_nothrow_move_constructible_v<multi::array<int, 4>::iterator>);
		static_assert(std::is_trivially_copyable_v<multi::array<int, 4>::iterator>);
		static_assert(std::is_swappable_v<multi::array<int, 4>::iterator>);
		// static_assert( std::is_trivially_relocatable_v<multi::array<int, 4>::iterator>);  // <==========
	}

	return boost::report_errors();
}
