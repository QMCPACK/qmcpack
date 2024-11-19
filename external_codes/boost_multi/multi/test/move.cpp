// Copyright 2020-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 10.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/core/lightweight_test.hpp>

#include <boost/multi/array.hpp>  // for array, apply, array_types<>::ele...

#include <algorithm>  // for copy, equal, fill_n, move
#include <iterator>   // for size, back_insert_iterator, back...
#include <memory>     // for unique_ptr, make_unique, allocat...
// IWYU pragma: no_include <type_traits>  // for remove_reference<>::type
// IWYU pragma: no_include <map>
// IWYU pragma: no_include <set>
// IWYU pragma: no_include <stack>
#include <utility>  // for move, swap
#include <vector>   // for vector, operator==, vector<>::va...

// IWYU pragma: no_include <pstl/glue_algorithm_defs.h>       // for move

namespace multi = boost::multi;

#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

namespace {

void move_element_1d_array() {
	{
		std::vector<multi::array<int, 1> > varr(3, multi::array<int, 1>({5}, 99), {});
		multi::array<int, 2> marr({3, 5}, 99);
		marr[0] = std::move(varr[0]);
		marr[1] = std::move(varr[1]);
		marr[2] = std::move(varr[2]);

		BOOST_TEST( marr[0][0] == 99 );
		BOOST_TEST( !varr[0].empty() );
	}
	{
		multi::array<std::vector<double>, 1> arr(10, std::vector<double>(5, {}, {}));
		multi::array<std::vector<double>, 1> brr(10, {}, {});

		std::copy_n( std::move(arr).begin(), brr.size(), brr.begin() );
		BOOST_TEST( arr[0].empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved)
		BOOST_TEST( brr[0].size() == 5 );
	}
	{
		multi::array<std::vector<double>, 1> arr(10, std::vector<double>(5, {}, {}));
		multi::array<std::vector<double>, 1> brr(10, {}, {});

		std::copy_n( arr.mbegin(), arr.size(), brr.begin() );
		BOOST_TEST( arr[0].empty() );
		BOOST_TEST( brr[0].size() == 5 );
	}
	{
		multi::array<std::vector<double>, 1> arr(10, std::vector<double>(5, {}, {}));

		using std::move;  // not necessary, just testing if it works
		auto vec = move(arr({2, 6}))[0];
		BOOST_TEST( vec.size() == 5 );
		BOOST_TEST( arr[2].empty() );
	}
	{
		multi::array<std::vector<double>, 1> arr(10, std::vector<double>(5, {}, {}));

		auto mbeg = arr({2, 6}).mbegin();
		auto vec = *mbeg;
		BOOST_TEST( vec.size() == 5 );
		BOOST_TEST( arr[2].empty() );
	}
	{
		multi::array<std::vector<double>, 1> arr(10, std::vector<double>(5, {}, {}));

		auto mbeg = arr({2, 6}).mbegin();
		auto vec = *mbeg;
		BOOST_TEST( vec.size() == 5 );
		BOOST_TEST( arr[2].empty() );
	}
	{
		multi::array<std::vector<double>, 1> arr(10, std::vector<double>(5, {}, {}));

		std::vector<std::vector<double>> out_vec(4, {}, {});
		std::copy(arr({2, 6}).mbegin(), arr({2, 6}).mend(), out_vec.begin());
		BOOST_TEST( out_vec[0].size() == 5 );
		BOOST_TEST( arr[2].empty() );
	}
	{
		multi::array<std::vector<double>, 1> arr(10, std::vector<double>(5, {}, {}));

		std::vector<std::vector<double>> out_vec(4, {}, {});
		std::copy(multi::move(arr({2, 6})).begin(), multi::move(arr({2, 6})).end(), out_vec.begin());
		BOOST_TEST( out_vec[0].size() == 5 );
		BOOST_TEST( arr[2].empty() );
	}
	{
		multi::array<std::vector<double>, 1> arr(10, std::vector<double>(5, {}, {}));

		std::vector<std::vector<double>> out_vec(4, {}, {});
		auto&& marr62 = multi::move(arr({2, 6}));
		std::copy(std::move(marr62).begin(), std::move(marr62).end(), out_vec.begin());  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved)
		BOOST_TEST( out_vec[0].size() == 5 );
		BOOST_TEST( arr[2].empty() );
	}
	{
		multi::array<std::vector<double>, 1> arr(10, std::vector<double>(5, {}, {}));

		std::vector<std::vector<double>> out_vec(4, {}, {});
		auto&& marr62 = arr({2, 6});
		std::copy(multi::move(marr62).begin(), multi::move(marr62).end(), out_vec.begin());  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved)
		BOOST_TEST( out_vec[0].size() == 5 );
		BOOST_TEST( out_vec[1].size() == 5 );
		BOOST_TEST( arr[2].empty() );
	}
	{
		multi::array<std::vector<double>, 1> arr(10, std::vector<double>(5, {}, {}));

		std::vector<std::vector<double>> out_vec(4, {}, {});
		auto&& marr62 = multi::move(arr({2, 6}));
		std::copy(marr62.begin(), marr62.end(), out_vec.begin());
		BOOST_TEST( out_vec[0].size() == 5 );
		BOOST_TEST( !arr[2].empty() );
	}
}

void move_element_2d_array() {
	multi::array<std::vector<double>, 2> arr({10, 10}, std::vector<double>(5, {}, {}));

	using std::move;
	auto vec = move(arr({2, 6}, {2, 6}))[0][0];
	BOOST_TEST( vec.size() == 5 );
	BOOST_TEST( arr[2][2].empty() );
}

void move_element_1d_total_array() {
	{
		multi::array<std::vector<double>, 1> arr(10, std::vector<double>(5, {}, {}));

		auto vec = std::move(arr)[2];
		BOOST_TEST( vec.size() == 5 );
		BOOST_TEST( arr[2].empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved)
	}
	{
		multi::array<std::vector<double>, 1> arr(10, std::vector<double>(5, {}, {}));

		using std::move;
		auto vec = move(arr)[2];
		BOOST_TEST( vec.size() == 5 );
		BOOST_TEST( arr[2].empty() );
	}
}
}  // namespace

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	move_element_1d_array();
	move_element_2d_array();

	move_element_1d_total_array();

	BOOST_AUTO_TEST_CASE(move_unique_ptr_1D) {
		{
			multi::array<std::unique_ptr<int>, 1> arr(multi::extensions_t<1>{10});
			arr[1] = std::make_unique<int>(42);

			multi::array<std::unique_ptr<int>, 1> arr2(multi::extensions_t<1>{10});
			std::move(arr.begin(), arr.end(), arr2.begin());

			BOOST_TEST( !arr[1] );
			BOOST_TEST(  arr2[1] );
			BOOST_TEST( *arr2[1] == 42 );
		}
		{
			multi::array<std::unique_ptr<int>, 1> arr(multi::extensions_t<1>{10});
			arr[1] = std::make_unique<int>(42);

			multi::array<std::unique_ptr<int>, 1> arr2 = std::move(arr);
			BOOST_TEST(  arr.is_empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved,clang-analyzer-cplusplus.Move)
			BOOST_TEST(  arr2[1] );
			BOOST_TEST( *arr2[1] == 42 );
		}
		{
			multi::array<std::unique_ptr<int>, 1> arr(multi::extensions_t<1>{10});
			arr[1] = std::make_unique<int>(42);

			multi::array<std::unique_ptr<int>, 1> arr2;  // (multi::extensions_t<1>{10});
			arr2 = std::move(arr);
			BOOST_TEST(  arr.is_empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved,clang-analyzer-cplusplus.Move)
			BOOST_TEST(  arr2[1] );
			BOOST_TEST( *arr2[1] == 42 );
		}
		// {
		//  multi::array<std::unique_ptr<int>, 1> arr(multi::extensions_t<1>{10});
		//  arr[1] = std::make_unique<int>(42);

		//  multi::array<std::unique_ptr<int>, 1> arr2(multi::extensions_t<1>{10});
		//  //  arr2() = arr();  // fails to compile, elements are not copy assignable
		//  arr2() = arr().element_moved();
		//  BOOST_TEST( !arr[1] );
		//  BOOST_TEST(  arr2[1] );
		//  BOOST_TEST( *arr2[1] == 42 );
		// }
	}

	BOOST_AUTO_TEST_CASE(multi_swap) {
#ifndef _MSC_VER  // problems with 14.3 c++17
		multi::array<int, 2> arr({3, 5}, 990);
		multi::array<int, 2> arr2({7, 11}, 880);
#else
		multi::array<int, 2> arr(multi::extensions_t<2>{3, 5}, 990);
		multi::array<int, 2> arr2(multi::extensions_t<2>{7, 11}, 880);
#endif

		swap(arr, arr2);

		BOOST_TEST( size(arr) == 7 );
		BOOST_TEST( arr[1][2] == 880 );
		BOOST_TEST( arr2[1][2] == 990 );
	}

	BOOST_AUTO_TEST_CASE(multi_std_swap) {
#ifndef _MSC_VER  // problems with 14.3 c++17
		multi::array<int, 2> arr({3, 5}, 990);
		multi::array<int, 2> arr2({7, 11}, 880);
#else
		multi::array<int, 2> arr(multi::extensions_t<2>{3, 5}, 990);
		multi::array<int, 2> arr2(multi::extensions_t<2>{7, 11}, 880);
#endif

		using std::swap;
		swap(arr, arr2);

		BOOST_TEST( size(arr) == 7 );
		BOOST_TEST( arr[1][2] == 880 );
		BOOST_TEST( arr2[1][2] == 990 );
	}

	BOOST_AUTO_TEST_CASE(multi_array_clear) {
		multi::array<int, 2> arr({10, 10}, 990);

		arr.clear();

		BOOST_TEST(arr.is_empty());

		arr.reextent({20, 20}, 990);
		// BOOST_TEST(! arr.is_empty());

		// clear(arr).reextent({30, 30}, 88.0);
		// BOOST_TEST(arr[15][15] == 88.0);
	}

	BOOST_AUTO_TEST_CASE(multi_array_move) {
		std::vector<multi::array<int, 2>> Av(10, multi::array<int, 2>({4, 5}, 990));  // std::vector NOLINT(fuchsia-default-arguments-calls)
		multi::array<int, 2>              arr2(std::move(Av[0]), std::allocator<int>{});

		BOOST_TEST( is_empty(Av[0]) );
		BOOST_TEST( size(arr2) == 4 );
		BOOST_TEST( arr2[1][2] == 990 );
	}

	BOOST_AUTO_TEST_CASE(multi_array_move_into_vector) {
		std::vector<multi::array<int, 2>> Av(10, multi::array<int, 2>({4, 5}, 990));  // NOLINT(fuchsia-default-arguments-calls)
		std::vector<multi::array<int, 2>> Bv;
		Bv.reserve(Av.size());  // NOLINT(fuchsia-default-arguments-calls)

		std::move(begin(Av), end(Av), std::back_inserter(Bv));

		BOOST_TEST( size(Bv) == size(Av) );
		BOOST_TEST( is_empty(Av[4]) );
		BOOST_TEST( size(Bv[5]) == 4 );
		BOOST_TEST( Bv[5][1][2] == 990 );
	}

	BOOST_AUTO_TEST_CASE(multi_array_move_into_vector_reserve) {
		std::vector<multi::array<int, 2>> Av(10, multi::array<int, 2>({4, 5}, 990));  // NOLINT(fuchsia-default-arguments-calls)
		std::vector<multi::array<int, 2>> Bv;
		Bv.reserve(Av.size());

		//  for(auto& v: Av) Bv.emplace_back(std::move(v), std::allocator<int>{});  // segfaults nvcc 11.0 but not nvcc 11.1
		std::move(begin(Av), end(Av), std::back_inserter(Bv));

		BOOST_TEST( size(Bv) == size(Av) );
		BOOST_TEST( is_empty(Av[4]) );
		BOOST_TEST( size(Bv[5]) == 4 );
		BOOST_TEST( Bv[5][1][2] == 990 );
	}

	BOOST_AUTO_TEST_CASE(multi_array_move_into_vector_move) {
		std::vector<multi::array<int, 2>> Av(10, multi::array<int, 2>({4, 5}, 990));  // std::vector NOLINT(fuchsia-default-arguments-calls)
		std::vector<multi::array<int, 2>> Bv = std::move(Av);

		Av.clear();

		BOOST_TEST( size(Av) == 0 );
		BOOST_TEST( size(Bv) == 10 );
		BOOST_TEST( size(Bv[5]) == 4 );
		BOOST_TEST( Bv[5][1][2] == 990 );
	}

	BOOST_AUTO_TEST_CASE(multi_array_move_array) {
		multi::array<std::vector<int>, 2> arr({10, 10}, std::vector<int>(5));  // std::vector NOLINT(fuchsia-default-arguments-calls)
		auto                              arr2 = std::move(arr);

		// NOLINTNEXTLINE(bugprone-use-after-move,hicpp-invalid-access-moved,clang-analyzer-cplusplus.Move) test deterministic moved from state
		BOOST_TEST( arr .   empty() );

		// NOLINTNEXTLINE(bugprone-use-after-move,hicpp-invalid-access-moved,clang-analyzer-cplusplus.Move) test deterministic moved from state
		BOOST_TEST( arr .is_empty() );
		BOOST_TEST( arr2.size() == 10 );
	}

	BOOST_AUTO_TEST_CASE(multi_array_move_elements) {
		multi::array<std::vector<int>, 1> arr({10}, std::vector<int>(5));  // std::vector NOLINT(fuchsia-default-arguments-calls)

		std::vector<std::vector<int>> sink(5);  // std::vector NOLINT(fuchsia-default-arguments-calls)

		auto* ptr1 = arr[1].data();

		std::copy(arr({0, 5}).element_moved().begin(), arr({0, 5}).element_moved().end(), sink.begin());
		BOOST_TEST(  arr[1].empty() );
		BOOST_TEST( !arr[5].empty() );

		BOOST_TEST( sink[1].data() == ptr1 );
	}

	BOOST_AUTO_TEST_CASE(multi_array_move_elements_range) {
		multi::array<std::vector<int>, 1> arr({10}, std::vector<int>(5));  // std::vector NOLINT(fuchsia-default-arguments-calls)

		std::vector<std::vector<int>> sink(5);  // NOLINT(fuchsia-default-arguments-calls)

		auto* ptr1 = arr[1].data();

		std::copy(arr({0, 5}).element_moved().elements().begin(), arr({0, 5}).element_moved().elements().end(), sink.begin());
		BOOST_TEST(     arr[1].empty() );
		BOOST_TEST( !arr[5].empty() );

		BOOST_TEST( sink[1].data() == ptr1 );
	}

	BOOST_AUTO_TEST_CASE(multi_array_move_elements_to_array) {
		multi::array<std::vector<int>, 1> arr({10}, std::vector<int>(5, 990));  // std::vector NOLINT(fuchsia-default-arguments-calls)
		BOOST_TEST( arr.size() == 10 );
		multi::array<std::vector<int>, 1> arr2({5}, {}, {});  // std::vector NOLINT(fuchsia-default-arguments-calls)

		auto* ptr1 = arr[1].data();

		arr2().elements() = arr({0, 5}).element_moved().elements();

		BOOST_TEST( arr2[1].size() == 5 );
		BOOST_TEST( arr2[1][4] == 990 );

		BOOST_TEST(     arr[1].empty() );
		BOOST_TEST( !arr[5].empty() );

		BOOST_TEST( arr2[1].data() == ptr1 );
	}

	BOOST_AUTO_TEST_CASE(move_range_vector_1D) {
		std::vector<std::vector<int>> arr(10, std::vector<int>{10, 20, 30});  // NOLINT(fuchsia-default-arguments-calls)
		std::vector<std::vector<int>> arr2(10);                               // NOLINT(fuchsia-default-arguments-calls)
		std::move(arr.begin(), arr.end(), arr2.begin());

		BOOST_TEST( arr2[0] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)
		BOOST_TEST( arr2[1] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)

		BOOST_TEST( arr[0].empty() );
		BOOST_TEST( arr[1].empty() );
	}

	BOOST_AUTO_TEST_CASE(copy_range_1D) {
		multi::array<std::vector<int>, 1> arr({3}, std::vector<int>{10, 20, 30});  // std::vector NOLINT(fuchsia-default-arguments-calls)
		BOOST_TEST( arr.size() == 3 );
		multi::array<std::vector<int>, 1> arr2({3}, std::vector<int>{});
		std::copy(arr.begin(), arr.end(), arr2.begin());

		BOOST_TEST( arr2[0] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)
		BOOST_TEST( arr2[1] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)

		BOOST_TEST( arr [0] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)
		BOOST_TEST( arr [1] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)
	}

	BOOST_AUTO_TEST_CASE(move_range_1D) {
		multi::array<std::vector<int>, 1> arr({3}, std::vector<int>{10, 20, 30});  // std::vector NOLINT(fuchsia-default-arguments-calls)
		BOOST_TEST( arr.size() == 3 );
		multi::array<std::vector<int>, 1> arr2({3}, std::vector<int>{});  // std::vector NOLINT(fuchsia-default-arguments-calls)
		std::move(arr.begin(), arr.end(), arr2.begin());

		BOOST_TEST( arr2[0] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)
		BOOST_TEST( arr2[1] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)

		BOOST_TEST( arr[0].empty() );
		BOOST_TEST( arr[1].empty() );
	}

	// BOOST_AUTO_TEST_CASE(move_range_1D_moved_begin) {
	//  multi::array<std::vector<int>, 1> arr({3}, std::vector<int>{10, 20, 30});  // std::vector NOLINT(fuchsia-default-arguments-calls)
	//  BOOST_TEST( arr.size() == 3 );
	//  multi::array<std::vector<int>, 1> arr2({3}, std::vector<int>{});  // std::vector NOLINT(fuchsia-default-arguments-calls)
	//  std::copy(arr.mbegin(), arr.mend(), arr2.begin());

	//  BOOST_TEST( arr2[0] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)
	//  BOOST_TEST( arr2[1] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)

	//  BOOST_TEST( arr[0].empty() );
	//  BOOST_TEST( arr[1].empty() );
	// }

	// BOOST_AUTO_TEST_CASE(copy_move_range) {
	//  multi::array<std::vector<int>, 2> arr({10, 20}, std::vector<int>{10, 20, 30});  // std::vector NOLINT(fuchsia-default-arguments-calls)
	//  multi::array<std::vector<int>, 2> arr2({10, 20}, std::vector<int>{});           // std::vector NOLINT(fuchsia-default-arguments-calls)

	//  std::copy(arr.mbegin(), arr.mend(), arr2.begin());

	//  BOOST_TEST( arr2[0][0] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)
	//  BOOST_TEST( arr2[0][1] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)

	//  BOOST_TEST( arr2[1][0] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)
	//  BOOST_TEST( arr2[1][1] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)

	//  BOOST_TEST( arr[0][0].empty() );
	//  BOOST_TEST( arr[0][1].empty() );

	//  BOOST_TEST( arr[1][0].empty() );
	//  BOOST_TEST( arr[1][1].empty() );
	// }

	BOOST_AUTO_TEST_CASE(copy_move_range_moved_begin) {
		multi::array<std::vector<int>, 2> arr({10, 20}, std::vector<int>{10, 20, 30});  // NOLINT(fuchsia-default-arguments-calls)
		multi::array<std::vector<int>, 2> arr2({10, 20}, std::vector<int>{});           // NOLINT(fuchsia-default-arguments-calls)

		std::copy(arr.element_moved().begin(), arr.element_moved().end(), arr2.begin());

		BOOST_TEST( arr2[0][0] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)
		BOOST_TEST( arr2[0][1] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)

		BOOST_TEST( arr2[1][0] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)
		BOOST_TEST( arr2[1][1] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)

		// BOOST_TEST( arr[0][0].size() == 0);
		// BOOST_TEST( arr[0][0].empty() );
		// BOOST_TEST( arr[0][1].empty() );

		// BOOST_TEST( arr[1][0].empty() );
		// BOOST_TEST( arr[1][1].empty() );
	}

	BOOST_AUTO_TEST_CASE(copy_move_range_moved_begin_block) {
		multi::array<std::vector<int>, 2> arr({10, 20}, std::vector<int>{10, 20, 30});  // NOLINT(fuchsia-default-arguments-calls)
		multi::array<std::vector<int>, 2> arr2({3, 5}, std::vector<int>{});

		std::copy(arr({5, 8}, {10, 15}).element_moved().begin(), arr({5, 8}, {10, 15}).element_moved().end(), arr2.begin());

		BOOST_TEST( arr2[0][0] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)
		BOOST_TEST( arr2[0][1] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)

		BOOST_TEST( arr2[1][0] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)
		BOOST_TEST( arr2[1][1] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)

		// BOOST_TEST( arr[5][10].empty() );
		// BOOST_TEST( arr[5][11].empty() );

		// BOOST_TEST( arr[6][10].empty() );
		// BOOST_TEST( arr[6][11].empty() );
	}

	BOOST_AUTO_TEST_CASE(move_reference_range) {
		multi::array<std::vector<int>, 2> arr({10, 20}, std::vector<int>{10, 20, 30});  // std::vector NOLINT(fuchsia-default-arguments-calls)
		multi::array<std::vector<int>, 2> arr2({10, 20}, std::vector<int>{});           // std::vector NOLINT(fuchsia-default-arguments-calls)

		//  arr2() = arr().element_moved();
		std::copy(arr().element_moved().begin(), arr().element_moved().end(), arr2().begin());

		BOOST_TEST( arr2[0][0] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)
		BOOST_TEST( arr2[0][1] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)

		BOOST_TEST( arr2[1][0] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)
		BOOST_TEST( arr2[1][1] == std::vector<int>({10, 20, 30}) );  // NOLINT(fuchsia-default-arguments-calls)

		// BOOST_TEST( arr[0][0].empty() );
		// BOOST_TEST( arr[0][1].empty() );

		// BOOST_TEST( arr[1][0].empty() );
		// BOOST_TEST( arr[1][1].empty() );
	}

	BOOST_AUTO_TEST_CASE(move_array_elements) {  // NOLINT(readability-function-cognitive-complexity)
		{
			auto arr  = multi::array<std::vector<int>, 1>({5}, std::vector<int>(7));  // std::vector NOLINT(fuchsia-default-arguments-calls)
			auto arr2 = std::move(arr);
			BOOST_TEST( arr2.size() == 5 );
			BOOST_TEST( arr2[0].size() == 7 );
			BOOST_TEST( arr.is_empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved,clang-analyzer-cplusplus.Move)
		}
		{
			auto arr = multi::array<std::vector<int>, 1>({5}, std::vector<int>(7));  // std::vector NOLINT(fuchsia-default-arguments-calls)

			std::vector<int> const v0 = std::move(arr[0]);
			BOOST_TEST( v0.size() == 7 );
			BOOST_TEST( arr[0].empty() );

			std::vector<int> const v1 = std::move(arr)[1];
			BOOST_TEST( v1.size() == 7 );
			BOOST_TEST( arr[1].empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved,clang-analyzer-cplusplus.Move) for test

			auto arr2 = multi::array<std::vector<int>, 1>({1}, std::vector<int>{});

			arr2({0, 1}) = arr({2, 3});
			BOOST_TEST( arr2[0].size() == 7 );
			BOOST_TEST( arr [2].size() == 7 );
		}
		{
			auto arr  = multi::array<std::vector<int>, 1>({5}, std::vector<int>(7));  // std::vector NOLINT(fuchsia-default-arguments-calls)
			auto arr2 = multi::array<std::vector<int>, 1>({5}, std::vector<int>{});

			arr2() = arr();
			BOOST_TEST( arr2[0].size() == 7 );
			BOOST_TEST( arr [2].size() == 7 );

			arr2() = std::move(arr)();
			BOOST_TEST( arr2[0].size() == 7 );
			// BOOST_TEST( arr [2].empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
		}
		{
			auto arr  = multi::array<std::vector<int>, 1>({5}, std::vector<int>(7));  // std::vector NOLINT(fuchsia-default-arguments-calls)
			auto arr2 = multi::array<std::vector<int>, 1>({5}, std::vector<int>{});

			arr2() = arr();
			BOOST_TEST( arr2[0].size() == 7 );
			BOOST_TEST( arr [0].size() == 7 );
		}
		{
			auto arr  = multi::array<std::vector<int>, 1>({5}, std::vector<int>(7));  // std::vector NOLINT(fuchsia-default-arguments-calls)
			auto arr2 = multi::array<std::vector<int>, 1>({5}, std::vector<int>{});

			arr2() = std::move(arr)();
			BOOST_TEST( arr2[0].size() == 7 );
			// BOOST_TEST( arr [0].empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
		}
		{
			auto arr  = multi::array<std::vector<int>, 1>({5}, std::vector<int>(7));  // std::vector NOLINT(fuchsia-default-arguments-calls)
			auto arr2 = multi::array<std::vector<int>, 1>({5}, std::vector<int>{});

			auto&& mAp = std::move(arr)();
			arr2()     = mAp;
			BOOST_TEST( arr2[0].size() == 7 );
			BOOST_TEST( arr [0].size() == 7 );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
		}
		{
			auto arr  = multi::array<std::vector<int>, 1>({5}, std::vector<int>(7));  // std::vector NOLINT(fuchsia-default-arguments-calls)
			auto arr2 = multi::array<std::vector<int>, 1>({5}, std::vector<int>{});

			arr2({0, 5}) = std::move(arr)();
			BOOST_TEST( arr2[0].size() == 7 );
			// BOOST_TEST( arr [0].empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
		}
		{
			auto arr  = multi::array<std::vector<int>, 1>({5}, std::vector<int>(7));  // std::vector NOLINT(fuchsia-default-arguments-calls)
			auto arr2 = multi::array<std::vector<int>, 1>({5}, std::vector<int>{});

			arr2() = arr.taked(5);
			BOOST_TEST( arr2[0].size() == 7 );
			BOOST_TEST( arr [0].size() == 7);  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
		}
		{
			auto arr  = multi::array<std::vector<int>, 1>({5}, std::vector<int>(7));  // std::vector NOLINT(fuchsia-default-arguments-calls)
			auto arr2 = multi::array<std::vector<int>, 1>({5}, std::vector<int>{});

			arr2() = std::move(arr).taked(5);
			BOOST_TEST( arr2[0].size() == 7 );
			// BOOST_TEST( arr [0].empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
		}
		{
			auto   arr  = multi::array<std::vector<int>, 1>({5}, std::vector<int>(7));  // std::vector NOLINT(fuchsia-default-arguments-calls)
			auto   arr2 = multi::array<std::vector<int>, 1>({5}, std::vector<int>{});
			auto&& mAt5 = std::move(arr).taked(5);
			arr2()      = mAt5;
			BOOST_TEST( arr2[0].size() == 7 );
			BOOST_TEST( arr [0].size() == 7 );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
		}
		{
			auto   arr  = multi::array<std::vector<int>, 1>({5}, std::vector<int>(7));  // std::vector NOLINT(fuchsia-default-arguments-calls)
			auto   arr2 = multi::array<std::vector<int>, 1>({5}, std::vector<int>{});
			auto&& mAt5 = std::move(arr).taked(5);
			arr2()      = mAt5;
			BOOST_TEST( arr2[0].size() == 7 );
			BOOST_TEST( arr [0].size() == 7 );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
		}
		{
			auto   arr  = multi::array<std::vector<int>, 1>({5}, std::vector<int>(7));  // std::vector NOLINT(fuchsia-default-arguments-calls)
			auto   arr2 = multi::array<std::vector<int>, 1>({5}, std::vector<int>{});
			auto&& mAt5 = std::move(arr).taked(5);
			arr2()      = std::move(mAt5);  // NOLINT(hicpp-move-const-arg,performance-move-const-arg) just testing
			BOOST_TEST( arr2[0].size() == 7 );
			// BOOST_TEST( arr[0].empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
		}
		{
			auto   arr  = multi::array<std::vector<int>, 1>({5}, std::vector<int>(7));  // std::vector NOLINT(fuchsia-default-arguments-calls)
			auto   arr2 = multi::array<std::vector<int>, 1>({5}, std::vector<int>{});
			auto&& mAt5 = std::move(arr).taked(5);
			arr2()      = std::move(mAt5).taked(5);
			BOOST_TEST( arr2[0].size() == 7 );
			// BOOST_TEST( arr[0].empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
		}
		{
			auto   arr    = multi::array<std::vector<int>, 1>({5}, std::vector<int>(7));  // std::vector NOLINT(fuchsia-default-arguments-calls)
			auto   arr2   = multi::array<std::vector<int>, 1>({5}, std::vector<int>{});
			auto&& mAt5   = std::move(arr).taked(5);
			auto&& mAt5t5 = std::move(mAt5).taked(5);
			arr2()        = mAt5t5;
			BOOST_TEST( arr2[0].size() == 7 );
			BOOST_TEST( arr[0].size() == 7 );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
		}
		{
			auto   arr  = multi::array<std::vector<int>, 1>({5}, std::vector<int>(7));  // std::vector NOLINT(fuchsia-default-arguments-calls)
			auto   arr2 = multi::array<std::vector<int>, 1>({5}, std::vector<int>{});
			auto&& mAt5 = std::move(arr).taked(5);
			arr2()      = std::move(mAt5).dropped(0);
			BOOST_TEST( arr2[0].size() == 7 );
			// BOOST_TEST( arr[0].empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
		}
		{
			auto arr  = multi::array<std::vector<int>, 1>({5}, std::vector<int>(7));  // std::vector NOLINT(fuchsia-default-arguments-calls)
			auto arr2 = multi::array<std::vector<int>, 1>({4}, std::vector<int>{});   // std::vector NOLINT(fuchsia-default-arguments-calls)
			arr2()    = std::move(arr).dropped(1);
			BOOST_TEST( arr2[0].size() == 7 );
			BOOST_TEST( arr [0].size() == 7 );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
										   // BOOST_TEST( arr [1].empty()     );    // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved) for testing
		}
	}

	BOOST_AUTO_TEST_CASE(multi_array_view_swap) {
		multi::array<int, 2> arrA({4, 5}, 99);
		multi::array<int, 2> arrB({4, 5}, 88);

		arrA().swap(arrB());

		BOOST_TEST( arrA[0][0] == 88 );
		BOOST_TEST( arrB[0][0] == 99 );
	}

	BOOST_AUTO_TEST_CASE(multi_array_view_swap_dimension_1) {
		multi::array<int, 2> arrA({4, 5}, 99);
		multi::array<int, 2> arrB({4, 5}, 88);

		arrA[0].swap(arrB[0]);

		BOOST_TEST( arrA[0][0] == 88 );
		BOOST_TEST( arrB[0][0] == 99 );

		BOOST_TEST( arrA[1][0] == 99 );
		BOOST_TEST( arrB[1][0] == 88 );
	}

	BOOST_AUTO_TEST_CASE(multi_array_view_swap_dimension_1_free) {
		multi::array<int, 2> arrA({4, 5}, 99);
		multi::array<int, 2> arrB({4, 5}, 88);

		swap(arrA[0], arrB[0]);

		BOOST_TEST( arrA[0][0] == 88 );
		BOOST_TEST( arrB[0][0] == 99 );

		BOOST_TEST( arrA[1][0] == 99 );
		BOOST_TEST( arrB[1][0] == 88 );
	}

	BOOST_AUTO_TEST_CASE(move_array_vector_1d) {
		multi::array<std::vector<double>, 1> arrA(10, std::vector<double>(5, 0.0, {}));

		BOOST_TEST( arrA[2].size() == 5 );
		{
			multi::array<std::vector<double>, 1> arrB = std::move(arrA);

			BOOST_TEST( arrA.empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved,clang-analyzer-cplusplus.Move)
			BOOST_TEST( arrB.size() == 10 );
			BOOST_TEST( arrB[2].size() == 5 );
		}
	}

	BOOST_AUTO_TEST_CASE(move_subarray_vector_1d) {
		multi::array<std::vector<double>, 1> arrA(10, std::vector<double>(5, 0.0));  // NOLINT(fuchsia-default-arguments-calls)

		BOOST_TEST( arrA[2].size() == 5 );
		{
			// this is not wrong but it is misleading since std::move is not doing anything on a reference-like type
			// using std::move;
			// multi::array<std::vector<double>, 1> arrB = move(arrA());  // NOLINT(clang-diagnostic-unqualified-std-cast-call)

			multi::array<std::vector<double>, 1> arrB = arrA();  // NOLINT(clang-diagnostic-unqualified-std-cast-call)

			BOOST_TEST( arrA.size() == 10 );
			BOOST_TEST( arrB.size() == 10 );
			BOOST_TEST( arrA[2].size() == 5 );
			BOOST_TEST( arrB[2].size() == 5 );
		}
	}

	// BOOST_AUTO_TEST_CASE(smart_move_subarray_vector_1d) {
	//  multi::array<std::vector<double>, 1> arrA(10, std::vector<double>(5));

	//  BOOST_TEST( arrA[2].size() == 5 );
	//  {
	//      using std::move;
	//      multi::array<std::vector<double>, 1> arrB = move(arrA());

	//      BOOST_TEST( arrA.size() == 10 );
	//      BOOST_TEST( arrB.size() == 10 );
	//      BOOST_TEST( arrA[2].size() == 0 );
	//      BOOST_TEST( arrB[2].size() == 5 );
	//  }
	// }

	return boost::report_errors();
}
