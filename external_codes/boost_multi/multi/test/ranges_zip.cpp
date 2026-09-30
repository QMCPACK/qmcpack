// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // IWYU pragma: keep

#include <boost/core/lightweight_test.hpp>

#if defined(TBB_FOUND) && !defined(__NVCC__)
#ifndef __clang__
#include <execution>
#endif
#endif

#if (__cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L)) && __has_include(<ranges>)
#include <ranges>  // IWYU pragma: keep
#endif

namespace multi = boost::multi;  // NOLINT(misc-unused-alias-decls)

auto main() -> int {  // NOLINT(bugprone-exception-escape)
#if defined(__cpp_lib_ranges_zip) and (__cpp_lib_ranges_zip >= 202110L)
	{
		multi::array<int, 1> arr1({10}, 5);
		multi::array<int, 1> arr2({10}, 6);

		auto&& zp = std::views::zip(arr1, arr2);

		auto [e1, e2] = zp[3];

		BOOST_TEST( e1 == 5 );
		BOOST_TEST( e2 == 6 );
	}
	{
		multi::array<int, 1> arr1({10}, 5);
		multi::array<int, 1> arr2({10}, 6);

		auto&& zp = std::views::zip(arr1(), arr2());

		auto&& [e1, e2] = zp[3];

		BOOST_TEST( e1 == 5 );
		BOOST_TEST( e2 == 6 );
	}
	{
		multi::array<int, 1> arr1({10}, 5);
		multi::array<int, 1> arr2({10}, 6);

		auto&& es1 = arr1.elements();
		auto&& es2 = arr2.elements();

		auto&& zp = std::views::zip(es1, es2);

		auto&& [e1, e2] = zp[3];

		BOOST_TEST( e1 == 5 );
		BOOST_TEST( e2 == 6 );
	}
	{
		multi::array<int, 1> arr1({10}, 5);
		multi::array<int, 1> arr2({10}, 6);

		auto&& zp = std::views::zip(arr1.elements(), arr2.elements());

		auto&& [e1, e2] = zp[3];

		BOOST_TEST( e1 == 5 );
		BOOST_TEST( e2 == 6 );
	}
	/////////////////////////////// 2D case
	{
		multi::array<int, 2> arr1({10, 10}, 5);
		multi::array<int, 2> arr2({10, 10}, 6);

		auto&& zp = std::views::zip(arr1.elements(), arr2.elements());

		auto&& [e1, e2] = zp[3];

		BOOST_TEST( e1 == 5 );
		BOOST_TEST( e2 == 6 );
	}
	{
		multi::array<int, 2> arr1({10, 10}, 5);
		multi::array<int, 2> arr2({10, 10}, 6);

		auto&& arr1ro0 = arr1[0];
		BOOST_TEST( arr1ro0[0] == 5 );

		auto&& es1 = arr1().elements();
		auto&& es2 = arr2().elements();

		auto&& zp = std::views::zip(es1, es2);

		auto&& [e1, e2] = zp[3];

		e1 = 55;
		e2 = 66;

		BOOST_TEST( e1 == 55 );
		BOOST_TEST( e2 == 66 );
	}
	{
		multi::array<int, 2> arr1({10, 10}, 5);
		multi::array<int, 2> arr2({10, 10}, 6);

		auto&& zp = std::views::zip(arr1().elements(), arr2().elements());

		auto&& [e1, e2] = zp[3];

		BOOST_TEST( e1 == 5 );
		BOOST_TEST( e2 == 6 );
	}
	{
		multi::array<int, 2> arr1({10, 10}, 5);
		multi::array<int, 2> arr2({10, 10}, 6);

		auto&& zp = std::views::zip(arr1[1], arr2[2]);

		auto&& [e1, e2] = zp[3];

		BOOST_TEST( e1 == 5 );
		BOOST_TEST( e2 == 6 );
	}
#endif

	return boost::report_errors();
}
