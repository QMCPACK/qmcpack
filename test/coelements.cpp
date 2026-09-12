// Copyright 2025 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <boost/core/lightweight_test.hpp>

#if defined(__cpp_lib_generator) && (__cpp_lib_generator >= 202207L)
#include <iostream>
#endif

namespace multi = boost::multi;

// https://godbolt.org/z/7MqxhWvz3

#if defined(__cpp_lib_generator) && (__cpp_lib_generator >= 202207L)
#if defined(__GNUC__)
#pragma GCC diagnostic ignored "-Wnull-dereference"
#endif

#include <generator>

template<class Arr2D>
std::generator<typename Arr2D::indexes>
co_extensions_elements(Arr2D const& arr2d) {
	auto const [is, js] = arr2d.extents();
	for(auto const i : is) {
		for(auto const j : js) {
			co_yield typename Arr2D::indexes{i, j};
		}
	}
}

template<class Arr2D>
std::generator<typename Arr2D::element_cref>
co_celements(Arr2D const& arr2d) {
	auto const [is, js] = arr2d.extents();
	for(auto const i : is) {
		for(auto const j : js) {
			co_yield arr2d[i][j];
		}
	}
}

#endif

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	multi::array<int, 2> const arr = {
		{0, 1, 2},
		{3, 4, 5}
	};

	BOOST_TEST( arr.extent()[1] == 1 );
	{
		auto const [i, j] = arr.extents()[1][2];
		BOOST_TEST( i == 1 );
		BOOST_TEST( j == 2 );
	}

#if defined(__cpp_lib_generator) && (__cpp_lib_generator >= 202207L)
	for(auto const& [i, j] : co_extensions_elements(arr)) {
		std::cout << i << ' ' << j << '\n';
	}
	{
		auto const [i, j] = *co_extensions_elements(arr).begin();
		BOOST_TEST( i == 0 );
		BOOST_TEST( j == 0 );
	}
	{
		auto const [i, j] = *(++co_extensions_elements(arr).begin());
		BOOST_TEST( i == 0 );
		BOOST_TEST( j == 1 );
	}
	{
		BOOST_TEST( *(++co_celements(arr).begin()) == 1 );
	}
#endif

	return boost::report_errors();
}
