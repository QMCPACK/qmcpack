// Copyright 2018-2025 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for array, layout_t, subarray, range

#include <boost/core/lightweight_test.hpp>

#include <algorithm>    // for copy
#include <array>        // for array, get
#include <deque>        // for deque, operator==
#include <iterator>     // for size, back_inserter
#include <numeric>      // for iota
#include <string>       // for operator""s, allocator, basic_st...
#include <tuple>        // for apply, make_tuple, tuple
#include <type_traits>  // for is_assignable_v
#include <utility>      // for move, forward
#include <vector>       // for vector, operator==

namespace multi = boost::multi;

namespace {
template<class Array1D>
void assign_elements_from_to(Array1D&& arr, std::deque<std::vector<double>>& dest) {  // NOLINT(google-runtime-references) dest is mutated
	// NOLINTNEXTLINE(bugprone-use-after-move,hicpp-invalid-access-moved)
	std::copy(std::forward<Array1D>(arr).begin(), std::forward<Array1D>(arr).end(), std::back_inserter(dest));
}
}  // end namespace

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	// BOOST_AUTO_TEST_CASE(empty_intersection)
	{
		multi::array<double, 1> arr({10});
		multi::array<double, 1> arr2;

		auto const is = intersection(arr.extent(), arr2.extent());
		BOOST_TEST( arr(is).is_empty() );
		arr2(is) = arr(is);

		BOOST_TEST( arr2(is) == arr(is) );
	}

	// BOOST_AUTO_TEST_CASE(multi_tests_element_access_with_tuple)
	{
		multi::array<char, 2> arr({3, 3}, 'k');

		std::array<int, 2> point = {
			{1, 2}
		};

		BOOST_TEST(  arr[point[0]][point[1]] ==  arr(1, 2) );
		BOOST_TEST( &arr(point[0], point[1]) == &arr[point[0]][point[1]] );

		BOOST_TEST( &arr[point[0]][point[1]] == &arr(point[0], point[1]) );
		BOOST_TEST( &arr(point[0], point[1]) == &arr.apply(point) );

		BOOST_TEST( &arr[point[0]][point[1]] == &std::apply(arr, point) );
		BOOST_TEST( &arr[point[0]][point[1]] == &     apply(arr, point) );
	}

	// BOOST_AUTO_TEST_CASE(multi_tests_extension_with_tuple)
	{
		{
			multi::array<double, 2>::extents_type const ext = {3, 4};

			multi::array<double, 2> const arr(ext, 44.0);

			BOOST_TEST( size(arr) == 3 );
		}
		{
			auto const [en, em] = std::make_tuple(3, 4);
			multi::array<double, 2> const arr({en, em}, 44.0);
			BOOST_TEST( size(arr) == 3 );
		}
		{
			auto arr = std::apply([](auto const&... szs) { return multi::array<double, 2>({szs...}, 55.0); }, std::make_tuple(3, 4));
			BOOST_TEST( size(arr) == 3 );

			using std::get;

			BOOST_TEST( get<0>(arr.sizes()) == 3 );
			BOOST_TEST( get<1>(arr.sizes()) == 4 );
		}
	}

	// BOOST_AUTO_TEST_CASE(multi_test_constness_reference)
	{
		multi::array<char, 2> const carr({10, 10}, '9');

		using std::size;

		BOOST_TEST( size( carr(1, {0, 3}) ) == 3 );

		BOOST_TEST( carr(1, {0, 3})[1] == '9' );
		static_assert(decltype(carr({0, 3}, 1))::rank_v == 1);
		BOOST_TEST( carr.sliced(0, 3).size() == 3 );

		BOOST_TEST( carr({0, 3}).rotated()[1].unrotated().size() == 3 );

		BOOST_TEST( carr({0, 3}, {0, 3})[1][1] == '9' );

		static_assert(!std::is_assignable_v<decltype(carr(1, {0, 3})[1]), double>);
	}

	// comparison for elements iterator
	{
		multi::array<int, 2> const arr({3, 3}, 99);
		auto const&                subarr = arr({0, 3}, {0, 3});

		BOOST_TEST(   subarr.elements().begin() == subarr.elements().begin()  );
		BOOST_TEST( !(subarr.elements().begin() != subarr.elements().begin()) );
		BOOST_TEST( !(subarr.elements().begin() < subarr.elements().begin())  );  // cppcheck-suppress duplicateExpression ; for testing purposes
	}

	// BOOST_AUTO_TEST_CASE(multi_test_stencil)
	{
		using namespace std::string_literals;  // NOLINT(build/namespaces) ""s

		// NOLINTBEGIN(misc-include-cleaner) bug in clang-tidy 18
		multi::array<std::string, 2> arr = {
			{"a"s, "b"s, "c"s, "d"s, "e"s},
			{"f"s, "g"s, "h"s, "f"s, "g"s},
			{"h"s, "i"s, "j"s, "k"s, "l"s},
		};
		// NOLINTEND(misc-include-cleaner) bug in clang-tidy 18

		BOOST_TEST(      size(arr) == 3                                            );
		BOOST_TEST(           arr.num_elements() == 3*5L                           );
		BOOST_TEST(           arr[1][2] == "h"                                     );

		BOOST_TEST(           arr          ({1, 3}, {2, 5}).size() == 2                  );
		BOOST_TEST(           arr          ({1, 3}, {2, 5}).extent().first() == 0          );
		BOOST_TEST(           arr          ({1, 3}, {2, 5}).num_elements() == 2*3L );
		BOOST_TEST(           arr          ({1, 3}, {2, 5}).num_elements() == 2*3L );
		BOOST_TEST(           arr          ({1, 3}, {2, 5})[0][0] == "h"           );
		BOOST_TEST(          &arr          ({1, 3}, {2, 5})[0][0] == &arr[1][2]    );

		BOOST_TEST(  arr().elements().size() == arr.num_elements() );

		BOOST_TEST( &arr({1, 3}, {2, 5}).elements()[0] == &arr(1, 2) );
		BOOST_TEST( &arr({1, 3}, {2, 5}).elements()[arr({1, 3}, {2, 5}).elements().size() - 1] == &arr(2, 4) );

		BOOST_TEST( &arr({1, 3}, {2, 5}).elements().front() == &arr(1, 2) );
		BOOST_TEST( &arr({1, 3}, {2, 5}).elements().back()  == &arr(2, 4) );

		auto beg = arr({1, 3}, {2, 5}).elements().begin();
		beg += (arr({1, 3}, {2, 5}).elements().size() - 1);
		BOOST_TEST( &*beg  == &arr(2, 4) );  // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay) bug in clang-tidy 14?

		{
			auto beg1 = arr({1, 3}, {2, 5}).elements().begin();
			auto end1 = arr({1, 3}, {2, 5}).elements().end();
			auto end2 = arr({1, 3}, {2, 5}).elements().end();

			for(; end1 != beg1; --end1) {  // NOLINT(altera-id-dependent-backward-branch,altera-unroll-loops)
			}
			BOOST_TEST( end1 == beg1 );  // cppcheck-suppress knownConditionTrueFalse ; for testing purposes

			for(; end1 != end2; ++end1) {  // NOLINT(altera-id-dependent-backward-branch,altera-unroll-loops)
			}
			BOOST_TEST( end1 == end2 );  // cppcheck-suppress knownConditionTrueFalse ;
		}

		BOOST_TEST( arr.elements().size() == arr.elements().end() - arr.elements().begin() );
		BOOST_TEST( arr.elements().size() - 1 == arr.elements().end() - (arr.elements().begin() + 1) );
		BOOST_TEST( arr.elements().size() - 1 == (arr.elements().end() - 1) - arr.elements().begin() );
		BOOST_TEST( arr.elements().size() - 2 == (arr.elements().end() - 1) - (arr.elements().begin() + 1) );
	}

	// BOOST_AUTO_TEST_CASE(empty_elements)
	{
		multi::array<int, 2> arr1;
		multi::array<int, 2> arr2;

		BOOST_TEST( arr1.elements().size() == 0 );
		BOOST_TEST( arr2.elements().size() == 0 );
		BOOST_TEST(   arr1.elements() == arr2.elements()  );
		BOOST_TEST( !(arr1.elements() != arr2.elements()) );
	}

	// BOOST_AUTO_TEST_CASE(multi_test_elements_1D)
	{
		multi::array<int, 1> arr = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
		BOOST_TEST( arr.size() == 10 );

		BOOST_TEST(  arr.elements().size() == 10 );
		BOOST_TEST( &arr.elements()[0] == &arr[0] );
		BOOST_TEST( &arr.elements()[9] == &arr[9] );

		BOOST_TEST(    arr.elements().begin() <  arr.elements().begin() + 1 );
		BOOST_TEST( !(arr.elements().begin() <  arr.elements().begin()) );
		BOOST_TEST( !(arr.elements().begin() + 1 <  arr.elements().begin() + 1) );

		BOOST_TEST(    arr.elements().begin() <  arr.elements().end()     );
		BOOST_TEST(    arr.elements().end()   >  arr.elements().begin()   );
		BOOST_TEST(    arr.elements().begin() != arr.elements().end()     );
		BOOST_TEST( !( arr.elements().begin() == arr.elements().end()   ) );

		BOOST_TEST(  arr().elements().begin() <  arr().elements().end() );
		BOOST_TEST(  arr().elements().begin() == arr().elements().begin() );

		BOOST_TEST( arr().elements().begin() <  arr().elements().end() || arr().elements().begin() == arr().elements().end() );
		BOOST_TEST( arr().elements().begin() <= arr().elements().end() );

		BOOST_TEST(  arr().elements().end()  >  arr().elements().begin() );
		BOOST_TEST(  arr().elements().end()  >= arr().elements().begin() );

		arr.elements() = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
		BOOST_TEST( arr[2] == 7 );
		BOOST_TEST( arr.elements()[2] == 7 );
		BOOST_TEST( arr.elements().end() - arr.elements().begin() == arr.size() );
		BOOST_TEST( arr.elements().begin() + arr.size() == arr.elements().end() );
		BOOST_TEST( &(*(arr.elements().begin() + arr.size())) == &(*arr.elements().end()) );

		auto beg = arr.elements().begin();
		beg += arr.size();
		BOOST_TEST( &(*beg) == &(*arr.elements().end()) );
	}

	// BOOST_AUTO_TEST_CASE(multi_test_elements_1D_as_range)
	{
		multi::array<int, 1> arr = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
		BOOST_TEST( arr.size() == 10 );

		arr().elements() = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
		BOOST_TEST( arr[2] == 7 );
		BOOST_TEST( arr.elements()[2] == 7 );

		arr(2) = 9;
		BOOST_TEST( arr[2] == 9 );
	}

	// BOOST_AUTO_TEST_CASE(elements_from_init_list_2D)
	{
		multi::array<int, 2> arr({3, 2});
		arr().elements() = {1, 2, 3, 4, 5, 6};
		BOOST_TEST(arr[1][0] == 3);

		arr.elements() = {10, 20, 30, 40, 50, 60};
		BOOST_TEST(arr[1][0] == 30);
	}

	// BOOST_AUTO_TEST_CASE(front_back_2D)
	{
		multi::array<int, 2> arr({3, 4});
		std::iota(arr.elements().begin(), arr.elements().end(), int{});

		BOOST_TEST(  arr.front()[2] ==  arr[0][2] );
		BOOST_TEST( &arr.front()[2] == &arr[0][2] );

		BOOST_TEST(  (*(arr.begin() + 2)).base() ==  arr[2].base() );
		BOOST_TEST(    (arr.begin() + 2)->base() ==  arr[2].base() );

		BOOST_TEST(  (*(arr.end() - 1)).base() ==  arr[2].base() );
		BOOST_TEST(    (arr.end() - 1)->base() ==  arr[2].base() );

		// auto const prv = std::prev(arr.end());
		// BOOST_TEST(  (*(prv)).base() ==  arr[2].base() );  // TODO(correaa) investigate why this fails in NVCC

		// BOOST_TEST(  (*(std::prev(arr.end()))).base() ==  arr[2].base() );  // TODO(correaa) investigate why this fails in NVCC
		// BOOST_TEST(  (*(std::prev(arr.end(), 1))).base() ==  arr[2].base() );  // TODO(correaa) investigate why this fails in NVCC

		BOOST_TEST(  arr.back ().base() ==  arr[2].base() );
		BOOST_TEST(  arr.back () ==  arr[2] );

		BOOST_TEST(  arr.back ()[2] ==  arr[2][2] );
		BOOST_TEST( &arr.back ()[2] == &arr[2][2] );
	}

	// BOOST_AUTO_TEST_CASE(front_back_1D)
	{
		multi::array<int, 1> arr({30}, int{});
		std::iota(arr.elements().begin(), arr.elements().end(), 0);

		BOOST_TEST(  arr.front() ==  arr[ 0] );
		BOOST_TEST( &arr.front() == &arr[ 0] );

		BOOST_TEST(  arr.back () ==  arr[29] );
		BOOST_TEST( &arr.back () == &arr[29] );
	}

	// BOOST_AUTO_TEST_CASE(elements_rvalues)
	{
		using movable_type = std::vector<int>;
		movable_type const movable_value(5, 99);  // NOLINT(fuchsia-default-arguments-calls)

		multi::array<movable_type, 1> arr = {movable_value, movable_value, movable_value};
		BOOST_TEST( arr.size() == 3 );

		movable_type const front = std::move(arr)[0];  // cppcheck-suppress accessMoved ; for testing purposes

		BOOST_TEST( front == movable_value );

		// cppcheck-suppress accessMoved ; for testing purposes
		BOOST_TEST( arr[0].empty()           );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing purposes

		BOOST_TEST( arr[1] == movable_value  );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing purposes

		std::move(arr)[1] = movable_value;
	}

	// BOOST_AUTO_TEST_CASE(elements_rvalues_nomove)
	{
		using movable_type = std::vector<double>;
		movable_type const movable_value(5, 99.0);  // NOLINT(fuchsia-default-arguments-calls)

		multi::array<movable_type, 1> arr = {movable_value, movable_value, movable_value};
		BOOST_TEST( arr.size() == 3 );

		std::deque<std::vector<double>> q1;

		assign_elements_from_to(arr, q1);

		BOOST_TEST( arr[0] == movable_value );

		std::deque<std::vector<double>> q2;

		assign_elements_from_to(std::move(arr), q2);

		//  BOOST_TEST( arr[0].empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing purposes

		BOOST_TEST( q1 == q2 );
	}

	// BOOST_AUTO_TEST_CASE(elements_rvalues_assignment)
	{
		std::vector<int> vec = {10, 20, 30};  // NOLINT(fuchsia-default-arguments-calls)

		std::move(vec) = std::vector<int>{30, 40, 50};  // NOLINT(fuchsia-default-arguments-calls)

		// cppcheck-suppress accessMoved ; for testing purposes
		std::move(vec)[1] = 990;  // it compiles  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing purposes

		// cppcheck-suppress accessMoved ; for testing purposes
		BOOST_TEST( vec[1] == 990 );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing purposes

		multi::array<int, 1>       arr1 = {10, 20, 30};
		multi::array<int, 1> const arr2 = {10, 20, 30};

		std::move(arr1) = arr2;  // this compiles TODO(correaa) should it?
	}

	// BOOST_AUTO_TEST_CASE(range_2)
	{
		multi::array<int, 3>       arr3({3, 4, 5}, 99);
		multi::array<int, 3> const brr3({2, 2, 5}, 88);

		// what(arr3, arr3({0, 2}, {0, 2}));
		// what(arr3, arr3.range({0, 2}), arr3.paren_aux_({0, 2}), arr3({0, 2}), arr3({0, 2}, {0, 2}));
		arr3({0, 2}, {0, 2}) = brr3;

		BOOST_TEST( arr3[0][0][0] == 88 );  // should not compile
	}

	{
		multi::array<double, 2> const A2D({3, 3}, 11);
		multi::array<double, 2>       B2D({2, 2}, 22);
		multi::array<double, 1>       v1D(3, 33);

		using boost::multi::_;
		v1D(_)    = A2D(_, 0);            // v1D() = A2D( _ , 0);
		v1D(_)    = A2D(0, _);            // v1D() = A2D( 0 )   ;
		B2D(_, _) = A2D({0, 2}, {0, 2});  // B2D() = A2D({0, 2}, {0, 2});
	}
	{
		auto A2D = multi::array<int, 2>{
			{1, 2},
			{3, 4}
		};
		BOOST_TEST( A2D[1][1] == 4 );

		A2D[1][1] = 44;

		BOOST_TEST( A2D[1][1] == 44 );  // cppcheck-suppress knownConditionTrueFalse ;  // test syntax

#if defined(__cpp_multidimensional_subscript) && (__cpp_multidimensional_subscript >= 202110L)
		BOOST_TEST(( A2D[1, 1] == 44 ));

		A2D[1, 1] = 444;
		BOOST_TEST(( A2D[1, 1] == 444 ));

		using boost::multi::_;
		BOOST_TEST(( &A2D[_, 1][1] == &A2D[1, 1] ));
#endif
	}
	{
		multi::array<int, 2> arr({(5 * 5) + 3, (7 * 7) + 11});

		auto&& barr = arr({0, 5L * 5L}, {0, 7L * 7L}).strided(5);

		BOOST_TEST( barr.size() == 5 );
		BOOST_TEST( barr.stride() == 5L*((7L*7L) + 11L) );

		{
			auto i0 = 3;
			auto j0 = 13;

			auto* ptr = &(barr[i0][j0]);

			auto dist = ptr - barr.base();

			auto i = dist / barr.layout().stride();  // get<1>(barr.layout().nelemss());

			BOOST_TEST( i == i0 );

			dist = dist % barr.layout().stride();

			using std::get;
			auto j = dist / get<1>(barr.layout().strides());

			BOOST_TEST( j == j0 );
		}
		{
			using std::get;
			auto [is, js] = barr.extents();
			for(auto i : is) {
				for(auto j : js) {  // NOLINT(altera-unroll-loops)
					BOOST_TEST(
						&barr[i][j] ==
						&barr[
							(&barr[i][j] - barr.base()) / get<0>(barr.strides())
						]
						[
							(&barr[i][j] - barr.base()) % get<0>(barr.strides()) / get<1>(barr.strides())
						]
					);
				}
			}
		}
	}
	{
		multi::array<int, 2> arr = {
			{1, 2, 3},
			{4, 5, 6}
		};

		BOOST_TEST( arr[1][1] == 5 );
		BOOST_TEST( multi::detail::invoke_square(arr, 1, 1) == 5 );
		BOOST_TEST( multi::detail::invoke_square(arr[1], 1) == 5 );
		BOOST_TEST( multi::detail::invoke_square(arr[1][1]) == 5 );
	}
	return boost::report_errors();
}
