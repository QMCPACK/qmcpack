// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>
#include <boost/multi/io.hpp>

#include <boost/core/lightweight_test.hpp>

#include <iostream>
#include <sstream>
// IWYU pragma: no_include <string>                           // for allocator, operator<<

#if __cplusplus >= 202302L || (defined(_MSVC_LANG) && _MSVC_LANG > 202002L)
#if __has_include(<print>)
#include <boost/multi/io/format.hpp>

#include <print>
#endif
#endif

namespace multi = boost::multi;

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	// matlab 1d
	{
		multi::array<double, 1> const arr = {1.0, 2.0, 3.0};

		std::ostringstream oss;
		multi::detail::print(oss, arr, "\t", "\t", "\a", "\a", 0);
		std::cout << "arr = " << oss.str() << "; end\n";
	}
	// matlab 2d
	{
		multi::array<double, 2> const arr = {
			{1,   3},
			{2, -10},
		};

		std::ostringstream oss;
		multi::detail::print(oss, arr, "\t", "\t", "\a", "\a", 0);
		std::cout << "arr = " << oss.str() << "; end\n";
	}
	// matlab 3d
	{
		// clang-format off
		multi::array<double, 3> const arr = {
			{
				{1.0, 2.0, 3.0},
				{4.0, 5.0, 6.0},
			},
			{
				{7.0, 8.0, 9.0},
				{10.0, 11.0, 12.0},
			},
		};
		// clang-format on

		std::ostringstream oss;
		multi::detail::print(oss, arr, "\t", "\t", "\a", "\a", 0);
		std::cout << "arr = " << oss.str() << "; end\n";
	}
	{
		std::ostringstream oss;
		oss << multi::array<int, 3>({3, 4, 5}, int{}).extent();
		BOOST_TEST( oss.str() == "[0, 3)" );
	}
	// fortran 2d
	{
		multi::array<double, 2> const arr = {
			{1,   3},
			{2, -10},
		};

		std::ostringstream oss;
		multi::detail::print(oss, arr, " ", " ", " ", "\a", 0);
		std::cout << "fortran 2d arr = " << oss.str() << "; end\n";
	}
	// julia 1d
	{
		multi::array<double, 1> const arr = {1, 2, 3};

		std::ostringstream oss;
		multi::detail::print(oss, arr, "[", ",", "]", "\t", 0);
		std::cout << "julia 1d arr = " << oss.str() << "; end\n";
	}
	// julia 2d
	{
		multi::array<double, 2> const arr = {
			{1,   3},
			{2, -10},
		};

		std::ostringstream oss;
		multi::detail::print(oss, arr, "[", ",", "]", "\a", 0);
		std::cout << "julia 2d arr = " << oss.str() << "; end\n";
	}
	// numpy 2d
	{
		multi::array<double, 2> const arr = {
			{1,   3},
			{2, -10},
		};
		std::ostringstream oss;
		multi::detail::print(oss, arr, "[", " ", "]", " ", 0);
		std::cout << "numpy 2d arr = " << oss.str() << "; end\n";
	}
	{
		multi::array<double, 1> const arr = {1.0, 2.0, 3.0};

		std::cout << "A1D = " << arr << "; no more, no less.\n";
	}
	{
		multi::array<double, 2> const arr = {
			{1.0, 2.0, 3.0},
			{4.0, 5.0, 6.0}
		};
		std::ostringstream oss;
		oss << "A2D = " << arr;
		BOOST_TEST( oss.str() == "A2D = {\n"
			"\t{1, 2, 3},\n"
			"\t{4, 5, 6}, \n"
			"}"
		);

		std::cout << "A2D = " << arr << "; no more, no less\n";
	}
	{
		// clang-format off
		multi::array<double, 3> const arr = {
			{
				{1.0, 2.0, 3.0},
				{4.0, 5.0, 6.0},
			},
			{
				{7.0, 8.0, 9.0},
				{10.0, 11.0, 12.0},
			},
		};
		// clang-format on

		std::cout << "A3D = " << arr << "; no more, no less\n";
	}
	{
		multi::array<double, 0> const arr{5.0};
		std::cout << "A0D = " << arr << "; no more, no less\n";
	}
	{
		multi::array<double, 4> const arr = {
			{
             {
             {1.0, 2.0, 3.0},
             {4.0, 5.0, 6.0},
             },
             {
             {7.0, 8.0, 9.0},
             {10.0, 11.0, 12.0},
             },
			 },
			{
             {
             {1.0, 2.0, 3.0},
             {4.0, 5.0, 6.0},
             },
             {
             {7.0, 8.0, 9.0},
             {10.0, 11.0, 12.0},
             },
			 }
		};
	}
	{
#if __cplusplus >= 202302L || (defined(_MSVC_LANG) && _MSVC_LANG > 202002L)
#if __has_include(<print>)
		multi::array<double, 2> const arr = {
			{1,   3},
			{2, -10},
		};

		std::print("{}", arr);
#endif
#endif
	}

	return boost::report_errors();
}
