// Copyright 2019-2025 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <boost/core/lightweight_test.hpp>

#include <tuple>    // for get  // NOLINT(misc-include-cleaner)
#include <utility>  // for forward

namespace multi = boost::multi;

namespace {
template<class Array>
auto flatted_last(Array&& arr) {
	return std::forward<Array>(arr).reversed().transposed().flatted().reversed();
}

template<class Array>
auto partitioned_last(Array&& arr, multi::ssize_t n) {
	return std::forward<Array>(arr).reversed().partitioned(n).transposed().transposed().reversed();
}
}  // end unnamed namespace

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	// BOOST_AUTO_TEST_CASE(multi_reversed_3d)
	{
		multi::array<double, 3> arr({30, 40, 50});

		BOOST_TEST( arr.reversed().size() == 50 );

		BOOST_TEST( & arr.reversed()[3][5][7] == &arr[7][5][3] );
	}

	// BOOST_AUTO_TEST_CASE(multi_reversed_4d)
	{
		using std::get;  // workaround no prior declaration in function call with explicit template arguments is a C++20 extension [-Wc++20-extensions]

		multi::array<double, 4> arr({13, 5, 7, 11});

		BOOST_TEST( arr.reversed().size() == 11 );

		BOOST_TEST( &arr.reversed()[1][2][3][4] == &arr[4][3][2][1] );

		BOOST_TEST( get<0>( arr.reversed().transposed().flatted().reversed().sizes() ) == 13 );
		BOOST_TEST( get<1>( arr.reversed().transposed().flatted().reversed().sizes() ) ==  5 );
		BOOST_TEST( get<2>( arr.reversed().transposed().flatted().reversed().sizes() ) == 77 );

		BOOST_TEST(( sizes(arr.reversed().transposed().flatted().reversed()) == decltype(sizes(arr.reversed().transposed().flatted().reversed())){13, 5, 77} ));
		BOOST_TEST( &arr.reversed().transposed().flatted().reversed()[1][2][ 5] == & arr[1][2][0][ 5] );
		BOOST_TEST( &arr.reversed().transposed().flatted().reversed()[1][2][10] == & arr[1][2][0][10] );
		BOOST_TEST( &arr.reversed().transposed().flatted().reversed()[1][2][11] == & arr[1][2][1][ 0] );
		BOOST_TEST( &arr.reversed().transposed().flatted().reversed()[1][2][12] == & arr[1][2][1][ 1] );

		BOOST_TEST( & flatted_last(arr)[1][2][12] == & arr[1][2][1][1] );
	}

	// BOOST_AUTO_TEST_CASE(multi_reversed_4d_partition_last)
	{
		multi::array<double, 4> arr({11, 5, 7, 12});

		BOOST_TEST( arr.reversed().size() == 12 );

		BOOST_TEST( & arr.reversed()[1][2][3][4] == &arr[4][3][2][1] );

		BOOST_TEST( & arr.reversed().partitioned(3).transposed().reversed()[1][2][3][0][1] == & arr[1][2][3][1] );
		BOOST_TEST( & arr.reversed().partitioned(3).transposed().reversed()[1][2][3][1][0] == & arr[1][2][3][4] );
		BOOST_TEST( & arr.reversed().partitioned(3).transposed().reversed()[1][2][3][1][1] == & arr[1][2][3][5] );

		BOOST_TEST( & partitioned_last(arr, 3)[1][2][3][1][1] == & arr[1][2][3][5] );
	}

	return boost::report_errors();
}
