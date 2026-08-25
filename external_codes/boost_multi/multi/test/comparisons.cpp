// Copyright 2018-2025 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for array, apply, operator!=, operat...

#include <boost/core/lightweight_test.hpp>

#include <algorithm>  // for equal
#include <complex>    // for complex, operator==
#include <iterator>   // for begin, end, cbegin, cend, size
#include <utility>    // for swap // IWYU pragma: keep  // NOLINT(misc-include-cleaner)

namespace multi = boost::multi;

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	// BOOST_AUTO_TEST_CASE(comparison_complex)
	{
		using complex = std::complex<double>;
		{
			multi::array<double, 1>  arr  = {1.0, 2.0, 3.0};
			multi::array<complex, 1> arr2 = {
				{1.0, 0.0},
				{2.0, 0.0},
				{3.0, 0.0},
			};

			BOOST_TEST( arr[1] == arr2[1] );
			BOOST_TEST( arr == arr2 );
			BOOST_TEST( !(arr != arr2) );
			BOOST_TEST( arr2 == arr );
			BOOST_TEST( !(arr2 != arr) );
		}
		{
			multi::array<double, 2> const arr = {
				{1.0, 2.0, 3.0},
				{4.0, 5.0, 6.0},
			};

			multi::array<complex, 2> const arr2 = {
				{{1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}},
				{{4.0, 0.0}, {5.0, 0.0}, {6.0, 0.0}},
			};

			BOOST_TEST( arr[1][1] == arr2[1][1] );
			BOOST_TEST( arr == arr2 );
			BOOST_TEST( !(arr != arr2) );
			BOOST_TEST( arr2 == arr );
			BOOST_TEST( !(arr2 != arr) );
			BOOST_TEST( std::equal(arr[1].begin(), arr[1].end(), begin(arr2[1]), end(arr2[1])) );
		}
	}

	// BOOST_AUTO_TEST_CASE(multi_comparisons_swap)
	{
		multi::array<double, 3> arr = {
			{ {1.2, 1.1},  {2.4, 1.0}},
			{{11.2, 3.0}, {34.4, 4.0}},
			{ {1.2, 1.1},  {2.4, 1.0}},
		};
		BOOST_TEST( arr[0] < arr[1] );

		swap(arr[0], arr[1]);
		BOOST_TEST( arr[1] < arr[0] );

		swap(arr[0], arr[1]);
		BOOST_TEST( arr[0] < arr[1] );
	}

	// BOOST_AUTO_TEST_CASE(comparisons_equality)
	{
		multi::array<double, 3> arr = {
			{ {1.2, 1.1},  {2.4, 1.0}},
			{{11.2, 3.0}, {34.4, 4.0}},
			{ {1.2, 1.1},  {2.4, 1.0}},
		};

		multi::array_ref<double, 3>  ref(arr.data_elements(), extents(arr));
		multi::array_cref<double, 3> cref(data_elements(arr), extents(arr));

		BOOST_TEST( arr ==  arr );
		BOOST_TEST( !(arr !=  arr) );
		BOOST_TEST( ref ==  arr );
		BOOST_TEST( !(ref !=  arr) );
		BOOST_TEST( ref == cref );
		BOOST_TEST( !(ref != cref) );

		BOOST_TEST( arr[0] ==  arr[2] );
		BOOST_TEST( ref[0] ==  arr[2] );
		BOOST_TEST( ref[0] == cref[2] );

		BOOST_TEST( !( arr[0] != arr[2]) );
		BOOST_TEST( !( ref[0] != ref[2]) );

		BOOST_TEST( !( arr[0] != arr[2]) );
		BOOST_TEST( !( ref[0] != ref[2]) );
	}

	// BOOST_AUTO_TEST_CASE(comparisons_ordering)
	{
		multi::array<int, 3> arr = {
			{ {12, 11},  {24, 10}},
			{{112, 30}, {344, 40}},
			{ {12, 11},  {24, 10}},
		};

		multi::array_ref<int, 3> ref(arr.data_elements(), extents(arr));

		multi::array_cref<int, 3> cref(data_elements(arr), extents(arr));

		BOOST_TEST(  arr[0]    <=  arr[1] );
		BOOST_TEST(  ref[0]    <=  arr[1] );

#if !defined(_MSC_VER)  // not working on msvc + cuda
		BOOST_TEST( cref[0]    <= cref[1] );
#endif

		BOOST_TEST(  arr[0][0] <= arr[0][1] );
		BOOST_TEST(  ref[0][0] <= arr[0][1] );

		BOOST_TEST(  arr[1][0][0] == 112 );
		BOOST_TEST(  ref[1][0][0] == 112 );
		BOOST_TEST( cref[1][0][0] == 112 );

		BOOST_TEST(  arr[0][0][0] == 12 );
		BOOST_TEST(  ref[0][0][0] == 12 );
		BOOST_TEST( cref[0][0][0] == 12 );

		swap(ref[0], ref[1]);

		BOOST_TEST(  begin(arr) <  end(arr) );
		BOOST_TEST( !(begin(arr) <  begin(arr)) );

		BOOST_TEST( arr.cbegin() < arr.cend() );

		BOOST_TEST( arr.end() -  arr.begin() == arr.size() );
	}
	{
		class range_1d {  // a type that has .extensions() and .elements() but is NOT a const_subarray
			multi::array<int, 1> impl_;

		 public:  // NOLINT(whitespace/indent)
			explicit range_1d(multi::array<int, 1> impl) : impl_(std::move(impl)) {}

			auto extents() const { return impl_.extents(); }
			auto elements() const { return impl_.elements(); }
		};

		multi::array<int, 1> const arr = {1, 2, 3};

		range_1d const brr(multi::array<int, 1>{1, 2, 3});
		range_1d const crr(multi::array<int, 1>{1, 2, 4});

		BOOST_TEST( !(arr != brr) );
		BOOST_TEST(   arr != crr  );
	}

	return boost::report_errors();
}
