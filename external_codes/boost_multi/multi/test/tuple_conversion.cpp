// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for array, dynamic_array, num_elements

#include <boost/core/lightweight_test.hpp>

#include <tuple>

namespace multi = boost::multi;

auto main() -> int {
	{
		multi::array<int, 1> const arr = {1, 2, 3, 4};
		BOOST_TEST( arr.size() == 4 );
	}
	{
		auto const tup = std::make_tuple(1, 2, 3);

		multi::array<int, 1> const arr(tup);

		BOOST_TEST( arr.size() == 3 );
	}
	{
		multi::array<int, 2> const arr2d({3, 2}, 0);
		multi::array<int, 1> const arr2d_sizes(arr2d.sizes());
		BOOST_TEST( arr2d_sizes.size() == 2 );
	}
	{
		auto const arr2d = multi::array<int, 2>({3, 2}, 0);

		auto const arr2d_sizes = std::apply([](auto... szs) { return multi::array<multi::array<int, 2>::size_type, 1>({szs...}); }, arr2d.sizes());

		BOOST_TEST( arr2d_sizes.size() == 2 );
	}
	{
		auto const arr2d = multi::array<int, 2>({3, 2}, 0);

		auto const arr2d_sizes = std::apply([](auto... szs) { return multi::array<multi::ssize_t, 1>({szs...}); }, arr2d.sizes());

		BOOST_TEST( arr2d_sizes.size() == 2 );
	}
	{
		auto const arr2d = multi::array<int, 2>({3, 2}, 0);

		auto const arr2d_sizes = apply([](auto... ss) { return multi::array<multi::ssize_t, 1>({ss...}); }, arr2d.sizes());

		BOOST_TEST( arr2d_sizes.size() == 2 );
	}

	return boost::report_errors();
}
