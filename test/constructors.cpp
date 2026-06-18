// Copyright 2019-2025 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <boost/core/lightweight_test.hpp>

// IWYU pragma: no_include <algorithm>
#include <complex>      // for complex  // IWYU pragma: keep  // bug in iwyu 18.1.7
#include <cstddef>      // for size_t
#include <iterator>     // for size
#include <string>       // for string
#include <type_traits>  // for is_copy_assignable_v, is_copy_co...
#include <utility>      // for move
#include <vector>       // for vector

namespace multi = boost::multi;

struct multiplies_bind1st {  // NOLINT(misc-use-internal-linkage)
	using complex = std::complex<double>;
	explicit multiplies_bind1st(multi::array<complex, 2>&& marr) : m_(std::move(marr)) {}  // this produces a bug in nvcc11.0
 private:
	multi::array<complex, 2> m_;
};

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	// BOOST_AUTO_TEST_CASE(multi_construct_1d)
	{
		multi::dynamic_array<int, 1> arr(multi::extents_t<1>{multi::iextension{10}}, 10);
		//  multi::dynamic_array<int, 1> arr(multi::array<int, 1>::extensions_type{10}, 10);
		BOOST_TEST( size(arr) == 10 );
		BOOST_TEST( arr[1] == 10 );
	}
	{
		multi::array<int, 2> const arr2(multi::extents_t<2>{3, 4});
		BOOST_TEST( arr2.size() == 3 );
	}
	// TODO(correaa) should this work
	// {
	// 	multi::array<int, 2> arr2{multi::extents_t<2>{3, 4}};
	// 	BOOST_TEST( arr2.size() == 3 );
	// 	BOOST_TEST( (~arr2).size() == 4 );
	// }
	// TODO(correaa) should this work
	// {
	// 	std::vector<std::vector<int>> vecs(3, std::vector<int>(4));
	// 	multi::array<int, 2> arr2{vecs.begin(), vecs.end()};
	// 	BOOST_TEST( arr2.size() == 3 );
	// 	BOOST_TEST( (~arr2).size() == 4 );
	// }

	// BOOST_AUTO_TEST_CASE(multi_constructors_inqnvcc_bug)
	{
		using complex = std::complex<double>;

		multi::array<complex, 2> marr({10, 10});
		multiplies_bind1st(std::move(marr));
	}

	// BOOST_AUTO_TEST_CASE(multi_constructors_1d)
	{
		{
			multi::array<double, 1> const arr(multi::extents_t<1>{multi::iextension{10}});
			BOOST_TEST( size(arr) == 10 );
		}
		{
			multi::array<int, 1> arr(multi::extents_t<1>{multi::iextension{10}}, int{});
			BOOST_TEST( size(arr) == 10 );
			BOOST_TEST( arr[5] == int{} );
		}
		{
			multi::array<int, 1> arr(multi::extents_t<1>{multi::iextension{10}}, int{});
			BOOST_TEST( size(arr) == 10 );
			BOOST_TEST( arr[5] == int{} );
		}
#if defined(__cpp_deduction_guides) && !defined(__NVCC__)
		{
			multi::array arr(multi::extents_t<1>({0, 10}), int{});
			BOOST_TEST( size(arr) == 10 );
			BOOST_TEST( arr[5] == int{} );
		}
		{
			// clang-format off
		multi::array arr({{0, 10}}, int{});
			// clang-format on
			BOOST_TEST( size(arr) == 10 );
			BOOST_TEST( arr[5] == int{} );
		}
		{
			multi::array arr({10}, int{});
			BOOST_TEST( size(arr) == 10 );
			BOOST_TEST( arr[5] == int{} );
		}
		{
			multi::array arr(10, int{});
			BOOST_TEST( size(arr) == 10 );
			BOOST_TEST( arr[5] == int{} );
		}
#endif
	}

	// BOOST_AUTO_TEST_CASE(multi_constructors_2d_ctad)
	{
#if defined(__cpp_deduction_guides) && !defined(__NVCC__)
		multi::array arr({10, 20}, int{});
		BOOST_TEST( size(arr) == 10 );
		BOOST_TEST( arr[5][6] == int{} );
#endif
	}

	// BOOST_AUTO_TEST_CASE(multi_constructors)
	{
		{
			// multi::array<double, 1> arr({10}); assert(size(A)==1); // warning in clang
		}
		{
			// multi::array<double, 1> arr({10}, double{}); assert(size(arr)==10); // warning in clang
		}
		{
			// multi::array<double, 1> arr({10}, double{}); assert(size(arr)==10); // warning in clang
		}
		{
			// multi::array<double, 1> arr({10}, 0.); assert(size(arr)==10); // warning in clang
		}
		{
			// multi::array<double, 1> arr({10}, {}); assert(size(arr)==10); // error ambiguous
		}
		{
			multi::array<int, 1> arr = {10};

			BOOST_TEST( size(arr) == 1 );
			BOOST_TEST( arr[0] == 10 );
		}

		{
			multi::array<std::size_t, 1> arr = {10};
			BOOST_TEST( size(arr) == 1 );
			BOOST_TEST( arr[0] == 10 );
		}
		{
			multi::array<int, 1> arr = {10};
			BOOST_TEST( size(arr) == 1 );
			BOOST_TEST( arr[0] == 10 );
		}
		{
			multi::array<int, 1> arr({10});
			BOOST_TEST( size(arr) == 1 );
			BOOST_TEST( arr[0] == 10 );
		}
		{
			multi::array<std::size_t, 1> arr({10});
			BOOST_TEST( size(arr) == 1 );
			BOOST_TEST( arr[0] == 10 );
		}
		{
			multi::array<int, 1> arr({10});
			BOOST_TEST( size(arr) == 1 );
			BOOST_TEST( arr[0] == 10 );
			//}{ multi::array<std::size_t, 1> arr({{10}})   ; assert( size(arr)==1 and arr[0]==10 );  // clang warns about double bracked
			//}{ multi::array<int        , 1> arr({{10}})   ; assert( size(arr)==1 and arr[0]==10 );  // clang warns about double bracked
			//}{ multi::array<double     , 1> arr({{10}})   ; assert( size(arr)==1 and arr[0]==10 );  // clang warns about double bracked
		}
		{
			multi::array<std::size_t, 1> const arr({0, 10});
			BOOST_TEST( size(arr) == 2 );
		}
		{
			multi::array<int, 1> const arr({0, 10});
			BOOST_TEST( size(arr) == 2 );
		}
		{
			multi::array<double, 1> const arr({0, 10});
			BOOST_TEST( size(arr) == 2 );
		}
		{
			using T = multi::array<std::string, 3>;

			static_assert(std::is_nothrow_destructible_v<T>);
			static_assert(std::is_default_constructible_v<T>);
			static_assert(std::is_nothrow_default_constructible_v<T>);

			static_assert(std::is_copy_constructible_v<T>);
			static_assert(std::is_copy_assignable_v<T>);

			// static_assert( std::is_nothrow_copy_constructible_v<T> );
			// static_assert( std::is_nothrow_copy_assignable_v<T> );

			static_assert(std::is_move_constructible_v<T>);
			static_assert(std::is_move_assignable_v<T>);

			static_assert(std::is_nothrow_move_constructible_v<T>);
			static_assert(std::is_nothrow_move_assignable_v<T>);
		}
	}

	// BOOST_AUTO_TEST_CASE(views_are_not_allocable)
	{
		// multi::array<double, 2> const AA = {{1.0, 2.0}, {3.0, 4.0}};
		// [[maybe_unused]] decltype(AA[0])* pp = new decltype(AA[0]){AA[0]};
		// delete pp;
	}

	// BOOST_AUTO_TEST_CASE(views_are_not_placeable)
	{
		// multi::array<double, 2> const AA = {{1.0, 2.0}, {3.0, 4.0}};
		// auto&& A0 = AA[0];
		// new(std::addressof(A0)) decltype(AA[0]){AA[1]};
	}

	// BOOST_AUTO_TEST_CASE(views_cannot_be_elements)
	{
		multi::array<double, 2> const AA = {
			{1.0, 2.0},
			{3.0, 4.0},
		};
		std::vector<decltype(AA[0])> vv;
		vv.emplace_back(AA[0]);
		vv.push_back(AA[0]);
		// auto&& A0 = AA[0];
		// vv.push_back(A0);
	}

	// BOOST_AUTO_TEST_CASE(views_cannot_be_elements2)
	{
		// multi::array<double, 2> const AA = {{1.0, 2.0}, {3.0, 4.0}};
		// std::vector<decltype(AA[0])> vv(3, AA[0]);
	}
	{
		using arri1d = multi::array<int, 1> const;

		// #ifdef __clang__
		// #pragma clang diagnostic push
		// #pragma clang diagnostic ignored "-Wbraced-scalar-init"
		// #endif

		// clang-format off
		// arri1d A1({3}, {11});    BOOST_TEST((A1 == arri1d{11, 11, 11} && A1.size() == 3));  // fair, clang warning: [-Wbraced-scalar-init]
		// arri1d A2{{3}, {11}};    BOOST_TEST((A2 == arri1d{3, 11}      && A2.size() == 2));  // fair, clang warning: [-Wbraced-scalar-init]
		// arri1d A3 = {{3}, {11}}; BOOST_TEST((A3 == arri1d{3, 11}      && A3.size() == 2));  // fair, clang warning: [-Wbraced-scalar-init]
#ifndef CPPCHECK
		arri1d A4({3}, 11);      BOOST_TEST((A4 == arri1d{11, 11, 11} && A4.size() == 3));  // good, no warning  // cppcheck-suppress [internalAstError];
		// arri1d A5{{3}, 11};      BOOST_TEST((A5 == arri1d{3, 11}      && A5.size() == 2));  // ok  , clang warning: [-Wbraced-scalar-init]
		// arri1d A6 = {{3}, 11};   BOOST_TEST((A6 == arri1d{3, 11}      && A6.size() == 2));  // fair, clang warning: [-Wbraced-scalar-init]
		arri1d A7(3, 11);        BOOST_TEST((A7 == arri1d{11, 11, 11} && A7.size() == 3));  // fair, no warning
		arri1d A8{3, 11};        BOOST_TEST((A8 == arri1d{3, 11}      && A8.size() == 2));  // fair, no warning
		arri1d A9 = {3, 11};     BOOST_TEST((A9 == arri1d{3, 11}      && A9.size() == 2));  // good, no warning
#endif

		arri1d B1(multi::extents_t<1>(3), 11); BOOST_TEST((B1.size() == 3));  // good, no warning
#ifndef __circle_build__  // deduced types not allowed in function parameters
		arri1d B2(multi::extents_t(3), 11);    BOOST_TEST((B1.size() == 3));  // good, no warning
#endif
#ifndef CPPCHECK
		// multi::array const C4({3}, 11);  BOOST_TEST((C4 == arri1d{11, 11, 11} && C4.size() == 3));  // good, no warning
		// multi::array const C7(3, 11);    BOOST_TEST((C7 == arri1d{11, 11, 11} && C7.size() == 3));  // fair, no warning
		// multi::array const C8{3, 11};    BOOST_TEST((C8 == arri1d{3, 11}      && C8.size() == 2));  // fair, no warning
		// multi::array const C9 = {3, 11}; BOOST_TEST((C9 == arri1d{3, 11}      && C9.size() == 2));  // good, no warning

		multi::array<int, 1> const C4({3}, 11);  BOOST_TEST((C4 == arri1d{11, 11, 11} && C4.size() == 3));  // good, no warning
		multi::array<int, 1> const C7(3, 11);    BOOST_TEST((C7 == arri1d{11, 11, 11} && C7.size() == 3));  // fair, no warning
		multi::array<int, 1> const C8{3, 11};    BOOST_TEST((C8 == arri1d{3, 11}      && C8.size() == 2));  // fair, no warning
		multi::array<int, 1> const C9 = {3, 11}; BOOST_TEST((C9 == arri1d{3, 11}      && C9.size() == 2));  // good, no warning
#endif
		// clang-format on

		// #ifdef __clang__
		// #pragma clang diagnostic pop
		// #endif
	}
	{
		std::vector<int> const vec = {1, 2, 3, 4, 5, 6};

		multi::array<int, 2> arr({2, 3}, vec.begin());
		BOOST_TEST( arr[1][0] == 4 );
	}
	return boost::report_errors();
}
