// Copyright 2019-2025 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for implicit_cast, explicit_cast

#include <boost/core/lightweight_test.hpp>

#include <algorithm>  // for for_each, equal
#include <array>      // for array
#include <cstdint>    // for int64_t
#include <iostream>   // for char_traits, operator<<, basic_o...
#include <iterator>   // for size
#include <memory>     // for allocator, unique_ptr
#include <numeric>    // for accumulate, iota
#ifdef BOOST_MULTI_HAS_SPAN
#include <span>  // for span
#endif
#include <string>       // for basic_string, operator""s, string
#include <tuple>        // for std::tie
#include <type_traits>  // for remove_reference, remove_const
// #include <typeinfo>     // for bad_cast
#include <utility>  // for as_const, move
#include <vector>   // for vector

namespace multi = boost::multi;

namespace boost::multi {

// NOLINTBEGIN(whitespace/indent_namespace) bug in cpplint
template<class T, boost::multi::dimensionality_type D, class Alloc = std::allocator<std::decay_t<T>>>
using Array = std::conditional_t<
	std::is_reference_v<T>,
	std::conditional_t<
		std::is_const_v<std::remove_reference_t<T>>,
		boost::multi::array_ref<std::remove_const_t<std::remove_reference_t<T>>, D> const&,
		boost::multi::array_ref<std::remove_reference_t<T>, D>&>,
	multi::array<T, D, Alloc>>;
// NOLINTEND(whitespace/indent_namespace)

}  // end namespace boost::multi

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
#endif

namespace {
auto f1d5(int const (&carr)[5]) -> int;   // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
auto f1d5(int const (&carr)[5]) -> int {  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
	return carr[1];
}

void f2d54(int const (&carr)[5][4]);   // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
void f2d54(int const (&carr)[5][4]) {  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
	BOOST_TEST(carr[0][1] == 1);
}

#ifdef __clang__
#pragma clang diagnostic pop
#endif

template<class T>
auto trace_array_deduce(multi::array<T, 2> const& arr) -> T {
	auto const& diag = arr.diagonal();
	return std::accumulate(diag.begin(), diag.end(), T{0});
}

template int trace_array_deduce(multi::array<int, 2> const&);

inline auto trace_separate_ref(multi::array_ref<int, 2> const& arr) -> int {
	auto const& diag = arr.diagonal();
	return std::accumulate(diag.begin(), diag.end(), 0);
}

inline auto trace_separate_sub(multi::subarray<int, 2> const& arr) -> int {
	auto const& diag = arr.diagonal();
	return std::accumulate(diag.begin(), diag.end(), 0);
}

}  // end unnamed namespace

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	{
		// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) for test
		int icarr[5] = {};  // cppcheck-suppress constVariable;

		multi::array_ref<int, 1> const iarrr(std::data(icarr), 5);
		static_assert(std::is_base_of_v<
					  std::random_access_iterator_tag,
					  std::iterator_traits<multi::array_ref<int, 1>::iterator>::iterator_category>);

		BOOST_TEST( &iarrr[0] == std::data(icarr) );
		// static_assert( std::is_base_of_v<
		//  std::contiguous_iterator_tag,
		//  std::iterator_traits<multi::array_ref<int, 1>::iterator>::iterator_category
		// >);
	}

	// BOOST_AUTO_TEST_CASE(array_ref_from_carray)
	{
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif

		// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) for test
		int arr[4][5] = {
			{  0,  10,  20,  30,  40},
			{ 50,  60,  70,  80,  90},
			{100, 110, 120, 130, 140},
			{150, 160, 170, 180, 190},
		};

		auto map = &multi::array_ref<int, 2>(arr);
		// multi::array_ptr<int, 2> const map{&arr};

		BOOST_TEST( &(*map).operator[](1)[1] == &arr[1][1] );  // cppcheck-suppress danglingTemporaryLifetime;
		BOOST_TEST( &map->operator[](1)[1] == &arr[1][1] );

		BOOST_TEST( (*&arr)[1][1] == 60 );  // cppcheck-suppress redundantPointerOp;

		multi::array_ref<int, 2>&& mar = *map;  // cppcheck-suppress [danglingTemporaryLifetime,danglingTempReference];

		BOOST_TEST( &mar[1][1] == &arr[1][1] );  // cppcheck-suppress [danglingTemporaryLifetime,danglingTempReference];

		mar[1][1] = 90;  // cppcheck-suppress [danglingTemporaryLifetime,danglingTempReference];

		BOOST_TEST( mar[1][1] == 90 );           // cppcheck-suppress [danglingTemporaryLifetime,danglingTempReference];
		BOOST_TEST( &mar[1][1] == &arr[1][1] );  // cppcheck-suppress danglingTempReference;

		auto const& a_const = arr;
		//  int const(&a_const)[4][5] = a;
		BOOST_TEST(&a_const[1][1] == &arr[1][1]);  // cppcheck-suppress knownConditionTrueFalse;

		static_assert(decltype(mar(2, {1, 3}))::rank_v == 1);  // cppcheck-suppress danglingTempReference;

		BOOST_TEST( size(mar(2, {1, 3})) == 2 );         // cppcheck-suppress danglingTempReference;
		BOOST_TEST( &mar(2, {1, 3})[1] == &arr[2][2] );  // cppcheck-suppress danglingTempReference;

		[[maybe_unused]] multi::array_ref<int, 2> const& cmar = *map;  // cppcheck-suppress danglingTempReference;

		// *(cmar.base()) = 99.0;
		// *(cmar[0].base()) = 88.0;
		// *(cmar.data_elements()) = 99.0;

#ifdef __clang__
#pragma clang diagnostic pop
#endif
	}

#ifndef __NVCC__
	// BOOST_AUTO_TEST_CASE(array_ref_test_ub)
	{
#if defined(__GNUC__) || defined(__NVCC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#endif

		// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays): test
		int arr[][4] = {
			{  0,  10,  20,  30},
			{ 50,  60,  70,  80},
			{100, 110, 120, 130},
			{150, 160, 170, 180},
		};

		multi::array_ref<int, 2> const map{arr};  // multi::array_ref<double, 2> const map(&arr[0][0], {4, 4});

		auto const& diag = map.diagonal();

		BOOST_TEST( diag.begin() != diag.end() );
// BOOST_TEST( std::accumulate(diag.begin(), diag.end(), 0) == 0 + 6 + 12 + 18 );
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
	}
#endif

	// BOOST_AUTO_TEST_CASE(array_ref_test_no_ub)
	{
		// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays): test
		int arr[5][4] = {
			{ 00,  10,  20,  30},
			{ 50,  60,  70,  80},
			{100, 110, 120, 130},
			{150, 160, 170, 180},
		};

		multi::array_ref<int, 2> const map(&arr[0][0], {4, 4});

		auto const& diag = map.diagonal();
		BOOST_TEST( diag.begin() != diag.end() );
		BOOST_TEST( std::accumulate(diag.begin(), diag.end(), 0) == 0 + 60 + 120 + 180 );
	}

	/*use iterator as pointer*/
	{
		multi::array<int, 1> const                                           arr = {1, 2, 3, 4, 5, 6};
		multi::array_ref<int, 2, multi::array<int, 1>::const_iterator> const arr2({2, 3}, arr.begin());
		BOOST_TEST(( arr2 == multi::array<int, 2>{
		{1, 2, 3},
		{4, 5, 6}
	}));
	}

	/*use iterator as pointer*/
	{
		// multi::array<int, 2> const arr = { {1, 2}, {3, 4}, {5, 6}, {7, 8}};
		// multi::array_ref<int, 2, multi::array<int, 2>::const_iterator> const arr2({2, 2}, arr.begin());
		// BOOST_TEST(( arr2 == multi::array{
		//  {1, 2, 3},
		//  {4, 5, 6}
		// }));
	}

	// BOOST_AUTO_TEST_CASE(array_ref_test_no_ub2)
	{
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif

		// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays): test
		int arr[][4] = {
			{},
			{00, 10, 20, 30},
			{50, 60, 70, 80},
			{100, 110, 120, 130},
			{150, 160, 170, 180},
			{},
		};

		multi::array_ref<int, 2> const map(&arr[1][0], {4, 4});

#ifdef __clang__
#pragma clang diagnostic pop
#endif

		auto const& diag = map.diagonal();
		BOOST_TEST( diag.begin() != diag.end() );
		BOOST_TEST( std::accumulate(diag.begin(), diag.end(), 0) == 0 + 60 + 120 + 180 );
	}

	// BOOST_AUTO_TEST_CASE(array_ref_test_allocated_ub_unique_ptr)
	{
		// NOLINTNEXTLINE(clang-analyzer-cplusplus.NewDeleteLeaks,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) for illustration
		std::unique_ptr<int const[]> const arrp(new int const[4UL * 4UL]{0, 10, 20, 30, 50, 60, 70, 80, 100, 110, 120, 130, 150, 160, 170, 180});  // cppcheck-suppress leakReturnValNotUsed;  // NOLINT(whitespace/line_length)

		BOOST_TEST( arrp[3] == 30 );
		{
			multi::array_ref<int, 2, int const*> const map(arrp.get(), {4, 4});

			auto const& diag = map.diagonal();

			BOOST_TEST( diag.begin() != diag.end() );
			BOOST_TEST( std::accumulate(diag.begin(), diag.end(), 00) == 00 + 60 + 120 + 180 );  // is this UB?
		}
	}

	// BOOST_AUTO_TEST_CASE(array_ref_1D_reindexed)
	{
		using namespace std::string_literals;  // NOLINT(build/namespaces) for literal "string"s

		// clang-format off
	std::array<std::string, 5> stdarr = {{"a"s, "b"s, "c"s, "d"s, "e"s}};  // NOLINT(misc-include-cleaner) bug in clang-tidy 18.1.3
		// clang-format on

		multi::array_ref<std::string, 1> mar = *&multi::array_ref<std::string, 1>(stdarr);  // *multi::array_ptr<std::string, 1>(&stdarr);

		BOOST_TEST( &mar[1] == &stdarr[1] );
		BOOST_TEST( sizes(mar.reindexed(1)) == sizes(mar) );

		auto diff = &(mar.reindexed(1)[1]) - &mar[0];
		BOOST_TEST( diff == 0 );

		BOOST_TEST( &mar.blocked(2, 4)[2] == &mar[2] );
		for(auto idx : mar.stenciled({2, 4}).extent()) {  // NOLINT(altera-unroll-loops)
			BOOST_TEST( &mar.stenciled({2, 4})[idx] == &mar[idx] );
		}

		// clang-format off
		multi::array<std::string, 1> arr({{2, 7}}, std::string{"xx"});  // NOLINT(fuchsia-default-arguments-calls) std::string
		// clang-format on

		BOOST_TEST( arr.size() == 5 );
		BOOST_TEST( arr.extent() == multi::iextension(2, 7) );
		arr[2] = "a";
		arr[3] = "b";
		arr[4] = "c";
		arr[5] = "d";
		arr[6] = "e";
		BOOST_TEST( std::equal(arr.begin(), arr.end(), mar.begin(), mar.end()) );  // NOLINT(llvm-use-ranges,modernize-use-ranges) for C++20
	}

	// BOOST_AUTO_TEST_CASE(array_ref_of_nested_std_array_reindexed)
	{
		// clang-format off
	std::array<std::array<double, 5>, 4> arr = {{
		{ { 0.0, 1.0, 2.0, 3.0, 4.0 } },
		{ { 5.0, 6.0, 7.0, 8.0, 9.0 } },
		{ { 10.0, 11.0, 12.0, 13.0, 14.0 } },
		{ { 15.0, 16.0, 17.0, 18.0, 19.0 } }
	}};
		// clang-format on

		// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays): test type
		multi::array_ref<double, 2> mar = *&multi::array_ref<double, 2>(arr);
		BOOST_TEST( &mar[1][1] == &arr[1][1] );
	}

	// BOOST_AUTO_TEST_CASE(array_ref_reindexed)
	{
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif

		// NOLINTNEXTLINE(hicpp-avoid-c-arrays, modernize-avoid-c-arrays, cppcoreguidelines-avoid-c-arrays): test
		double arr[4][5] = {
			{ 0.0,  1.0,  2.0,  3.0,  4.0},
			{ 5.0,  6.0,  7.0,  8.0,  9.0},
			{10.0, 11.0, 12.0, 13.0, 14.0},
			{15.0, 16.0, 17.0, 18.0, 19.0},
		};

		// NOLINTNEXTLINE(hicpp-avoid-c-arrays, modernize-avoid-c-arrays, cppcoreguidelines-avoid-c-arrays): special type
		multi::array_ref<double, 2> mar = *&multi::array_ref<double, 2>(arr);

		BOOST_TEST( &mar[1][1] == &arr[1][1] );

#ifdef __clang__
#pragma clang diagnostic pop
#endif

		BOOST_TEST( mar   .reindexed(1).size() == mar.size() );
		BOOST_TEST( mar[0].reindexed(1).size() == mar[0].size() );

		BOOST_TEST( mar.reindexed(1).sizes() == mar.sizes() );

		BOOST_TEST( &mar.reindexed(1)[1][0] == &mar[0][0] );

		BOOST_TEST( mar[0].reindexed(1).sizes() == mar[0].sizes() );
		BOOST_TEST( mar[0].reindexed(1).extent().first() == mar[0].extent().first () + 1 );
		BOOST_TEST( mar[0].reindexed(1).extent().last() == mar[0].extent().last() + 1 );

		auto diff = &mar[0].reindexed(1)[1] - &mar[0][0];
		BOOST_TEST( diff == 0 );

		//  BOOST_TEST( &(((mar<<1).reindexed(2)>>1).reindexed(1))[1][2] == &mar[0][0] );
		BOOST_TEST( &mar.reindexed(1, 2)[1][2] == &mar[0][0] );

		BOOST_TEST( &mar.reindexed(1)({1, 5})[1][0] == &mar[0][0] );

		BOOST_TEST(( sizes(mar.stenciled({2, 4})) == decltype(sizes(mar.stenciled({2, 4}))){2, 5} ));
		BOOST_TEST( &mar.stenciled({2, 4})[2][0] == &mar[2][0] );
		BOOST_TEST( &mar.stenciled({2, 4}, {1, 3})[2][1] == &mar[2][1] );

		//  BOOST_TEST( &mar[0][0] == mar.origin() ); // origin changed meaning in on 2020/Dec/16
		//  BOOST_TEST( mar.base() == mar.origin() );

		//  BOOST_TEST( mar.stenciled({2, 4}).origin() == mar.origin() );  // origin changed meaning in on 2020/Dec/16
		BOOST_TEST( mar.stenciled({2, 4}).base()   != mar.base()   );

		BOOST_TEST( &mar.stenciled({2, 4})[2][0] == mar.stenciled({2, 4}).base() );

		{
			// NOLINTBEGIN(fuchsia-default-arguments-calls) std::string ctor
			// multi::array<std::string, 2> arrB = {
			//  {"a", "b", "c", "d", "e"},
			//  {"f", "g", "h", "f", "g"},
			//  {"h", "i", "j", "k", "l"},
			// };
			// NOLINTEND(fuchsia-default-arguments-calls) std::string ctor
			// arrB.reindex(2);
			// BOOST_TEST( size(arrB) == 3 );
			// BOOST_TEST( arrB[2][0] == "a" );
		}
		{
			// NOLINTBEGIN(fuchsia-default-arguments-calls) std::string ctor
			// multi::array<std::string, 2> arrB = {
			//  {"a", "b", "c", "d", "e"},
			//  {"f", "g", "h", "f", "g"},
			//  {"h", "i", "j", "k", "l"},
			// };
			// NOLINTEND(fuchsia-default-arguments-calls) std::string ctor
			// arrB.reindex(2, 1);
			// BOOST_TEST( size(arrB) == 3 );
			// BOOST_TEST( arrB[2][1] == "a" );
		}
		{
			// using namespace std::string_literals;  // NOLINT(build/namespaces) for literal "string"s
			// multi::array<std::string, 2> arrB = (multi::array<std::string, 2>{
			//  {"a"s, "b"s, "c"s, "d"s, "e"s},
			//  {"f"s, "g"s, "h"s, "f"s, "g"s},
			//  {"h"s, "i"s, "j"s, "k"s, "l"s},
			// });  // .reindex(2, 1);  // std::string NOLINT(fuchsia-default-arguments-calls)

			// BOOST_TEST( arrB.reindex(2).extension() == multi::iextension(2, 5) );
			// auto exts = arrB.reindexed(2).extensions();

			// multi::array<std::string, 2> const arrC(exts);
			// BOOST_TEST( size(arrC) == 3 );
			// BOOST_TEST( size(arrC) == size(arrB) );

			// BOOST_TEST( arrC.extension().first()  == 2 );
			// BOOST_TEST( arrC.extension().last() == 5 );
		}
	}

	// BOOST_AUTO_TEST_CASE(array_ref_with_stencil)
	{
		std::array<std::array<double, 5>, 4> arr = {
			{{{0.0, 1.0, 2.0, 3.0, 4.0}},
			 {{5.0, 6.0, 7.0, 8.0, 9.0}},
			 {{10.0, 11.0, 12.0, 13.0, 14.0}},
			 {{15.0, 16.0, 17.0, 18.0, 19.0}}},
		};

		// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test type
		auto const& mar = *&multi::array_ref<double, 2>(arr);
		BOOST_TEST( mar.size() == 4 );

		// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test type
		multi::array<double, 2> ss = {
			{ 0.0, +1.0,  0.0},
			{+1.0, -4.0, +1.0},
			{ 0.0, +1.0,  0.0},
		};
		auto const& stencil = ss.reindexed(-1, -1);

		BOOST_TEST( stencil.size() == 3 );
		BOOST_TEST( &stencil[-1][-1] == stencil.base() );
	}

	// BOOST_AUTO_TEST_CASE(array_ref_1D_from_vector)
	{
		std::vector<double> vec = {1.0, 2.0, 3.0};  // std::vector NOLINT(fuchsia-default-arguments-calls)
		// clang-format off
	multi::array_ref<double, 1> aref({{1, 3}}, vec.data());
		// clang-format on
		BOOST_TEST( aref.extent() == multi::iextension(1, 3) );
		BOOST_TEST( &aref[1] == vec.data() );
	}

	// BOOST_AUTO_TEST_CASE(array_ref_2D_from_vector)
	{
		std::vector<double> vec = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};  // std::string NOLINT(fuchsia-default-arguments-calls)

		multi::array_ref<double, 2> aref({2, 3}, vec.data());

		BOOST_TEST( &aref[1][0] == &vec[3] );
	}

	// BOOST_AUTO_TEST_CASE(array_ref_from_subarray)
	{
		std::vector<std::int64_t> vec = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};

		multi::array_ref<std::int64_t, 2> aref({4, 4}, vec.data());

		multi::array<int, 2> barr = {
			{ 1,  2,  3,  4},
			{ 5,  6,  7,  8},
			{ 9, 10, 11, 12},
			{13, 14, 15, 16}
		};

		aref = barr();
	}

	// BOOST_AUTO_TEST_CASE(array_ref_2D_from_vector_with_offset)
	{
		// NOLINTNEXTLINE(fuchsia-default-arguments-calls)
		std::vector<double> vec = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

		multi::array_ref<double, 2> aref({multi::iextension(1, 3), multi::iextension(1, 4)}, vec.data());

		{
			auto exts                 = aref.extents();
			auto const [exts0, exts1] = exts;
			BOOST_TEST( exts0 == multi::iextension(1, 3) );

			BOOST_TEST( exts1.first()  == 1 );
			BOOST_TEST( exts1.last () == 4 );

			BOOST_TEST( exts1 == multi::iextension(1, 4) );

			BOOST_TEST( exts == decltype(exts)(multi::iextension(1, 3), multi::iextension(1, 4)) );
		}
		{
			auto exts = aref.extents();

			using std::get;

			BOOST_TEST( get<0>(exts) == multi::iextension(1, 3) );
			BOOST_TEST( get<1>(exts).first()  == 1 );
			BOOST_TEST( get<1>(exts).last () == 4 );
			BOOST_TEST( get<1>(exts) == multi::iextension(1, 4) );
			BOOST_TEST( exts == decltype(exts)(multi::iextension(1, 3), multi::iextension(1, 4)) );
		}
		{
			auto const exts = aref.extents();

			using std::get;

			BOOST_TEST( get<0>(exts) == multi::iextension(1, 3) );
			BOOST_TEST( get<1>(exts).first()  == 1 );
			BOOST_TEST( get<1>(exts).last () == 4 );
			BOOST_TEST( get<1>(exts) == multi::iextension(1, 4) );

			BOOST_TEST( exts == decltype(exts)(multi::iextension(1, 3), multi::iextension(1, 4)) );
		}
		{
			auto const exts = aref.extents();
			BOOST_TEST( exts.get<0>() == multi::iextension(1, 3) );
			BOOST_TEST( exts.get<1>().first()  == 1 );
			BOOST_TEST( exts.get<1>().last () == 4 );
			BOOST_TEST( exts.get<1>() == multi::iextension(1, 4) );
			BOOST_TEST(( exts == decltype(exts){multi::iextension(1, 3), multi::iextension(1, 4)} ));
		}
		{
			auto const exts = aref.extents();

			using std::get;

			BOOST_TEST( get<0>(exts) == multi::iextension(1, 3) );
			BOOST_TEST( get<1>(exts).first()  == 1 );
			BOOST_TEST( get<1>(exts).last () == 4 );
			BOOST_TEST( get<1>(exts) == multi::iextension(1, 4) );

			BOOST_TEST( exts == decltype(exts)(multi::iextension(1, 3), multi::iextension(1, 4)) );
		}
		{
			using std::get;  // workaround no prior declaration in function call with explicit template arguments is a C++20 extension [-Wc++20-extensions]

			BOOST_TEST( get<0>(aref.extents()) == multi::iextension(1, 3) );

			BOOST_TEST( get<1>(aref.extents()).first()  == 1 );
			BOOST_TEST( get<1>(aref.extents()).last () == 4 );

			BOOST_TEST( get<1>(aref.extents()) == multi::iextension(1, 4) );
			BOOST_TEST( aref.extents() == decltype(aref.extents())(multi::iextension(1, 3), multi::iextension(1, 4)) );
		}
		{
			using std::get;  // workaround no prior declaration in function call with explicit template arguments is a C++20 extension [-Wc++20-extensions]

			auto ss = aref.sizes();
			BOOST_TEST( get<0>(ss) == 2 );
			BOOST_TEST( get<1>(ss) == 3 );
			BOOST_TEST( ss == decltype(ss)(2, 3) );
		}
		{
			auto [nn, mm] = aref.sizes();
			BOOST_TEST( nn == 2 );
			BOOST_TEST( mm == 3 );
		}
		{
			using std::get;

			auto const ss = aref.sizes();
			BOOST_TEST( get<0>(ss) == 2 );
			BOOST_TEST( get<1>(ss) == 3 );
			BOOST_TEST( ss == decltype(ss)(2, 3) );
		}
		{
			using std::get;  // workaround no prior declaration in function call with explicit template arguments is a C++20 extension [-Wc++20-extensions]

			BOOST_TEST( get<0>(aref.sizes()) == 2 );
			BOOST_TEST( get<1>(aref.sizes()) == 3 );
			BOOST_TEST( aref.sizes() == decltype(aref.sizes())(2, 3) );
		}
		{
			auto const ss = aref.sizes();

			using std::get;  // workaround no prior declaration in function call with explicit template arguments is a C++20 extension [-Wc++20-extensions]

			BOOST_TEST( get<0>(ss) == 2 );
			BOOST_TEST( get<1>(ss) == 3 );
			BOOST_TEST( ss == decltype(ss)(2, 3) );
		}
		{
			using std::get;
			BOOST_TEST( get<0>(aref.sizes()) == 2 );
			BOOST_TEST( get<1>(aref.sizes()) == 3 );
			BOOST_TEST( aref.sizes() == decltype(aref.sizes())(2, 3) );
		}

		{
			using std::get;  // workaround no prior declaration in function call with explicit template arguments is a C++20 extension [-Wc++20-extensions]

			auto const ss = aref.sizes();
			BOOST_TEST( get<0>(ss) == 2 );
			BOOST_TEST( get<1>(ss) == 3 );
			BOOST_TEST( ss == decltype(ss)(2, 3) );
		}
		{
			using std::get;  // workaround no prior declaration in function call with explicit template arguments is a C++20 extension [-Wc++20-extensions]

			BOOST_TEST( get<0>(aref.sizes()) == 2 );
			BOOST_TEST( get<1>(aref.sizes()) == 3 );
			BOOST_TEST( aref.sizes() == decltype(aref.sizes())(2, 3) );
		}

		BOOST_TEST( &aref[1][1] == vec.data() );
	}

	// BOOST_AUTO_TEST_CASE(array_2D_with_offset)
	{
		multi::array<double, 2> const arr({multi::iextension(1, 3), multi::iextension(2, 5)}, 1.2);

		BOOST_TEST( arr.extent().first()  == 1 );
		BOOST_TEST( arr.extent().last () == 3 );
	}

	// BOOST_AUTO_TEST_CASE(array_ref_1D)
	{
		// clang-format off
		// NOLINTNEXTLINE(fuchsia-default-arguments-calls)
		std::array<std::string, 5> arr = {{"a", "b", "c", "d", "e"}};
		// clang-format on
		multi::array_ref<std::string, 1>&& mar = *&multi::array_ref<std::string, 1>{arr};
		// multi::Array<std::string(&)[1]> mar = *multi::Array<std::string(*)[1]>(&a);

		BOOST_TEST(  mar.extent().first() == 0 );
		BOOST_TEST(  mar.extent().last()  == 5 );

		auto&& mar1 = mar.reindexed(1);

		BOOST_TEST( mar1.extent().size() == mar.extent().size() );

		BOOST_TEST(  mar1.extent().first() == 1 );

		// BOOST_TEST( mar1.extent() == extent(mar1) );
		// BOOST_TEST(  extent(mar1).first() == 1 );  // TODO(correaa) failing with nvhpc 22.7 (but not later)

		BOOST_TEST(  mar1.extent().first() == 1 );
		BOOST_TEST(  mar1.extent().last()  == 6 );
		// BOOST_TEST( *extent(mar1).begin() == 1 );
		BOOST_TEST( *mar1.extent().begin() == 1 );

		BOOST_TEST( mar1.size() == mar.size() );
		BOOST_TEST( mar1.layout().extent().first() == 1 );

		BOOST_TEST( mar1.extent().first() == 1 );
		// BOOST_TEST( extent(mar1).first() == 1 );
		BOOST_TEST( &mar1[1]     == &arr[0] );     // NOLINT(readability-container-data-pointer) test access
		BOOST_TEST(  mar1.base() == &arr[0] );  // NOLINT(readability-container-data-pointer) test access
		BOOST_TEST(  mar1.base() ==  arr.data() );
	}

	// BOOST_AUTO_TEST_CASE(array_ref_original_tests_carray)
	{
		// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
		int darr[4][5] = {
			{1, 2},
			{2, 3},
		};
		multi::array_ref<int, 2>             ref(&darr[0][0], {4, 5});
		multi::array_ref<int, 2, int const*> cref(&darr[0][0], {4, 5});
		multi::array_ref<int const, 2>       crefc(&darr[0][0], {4, 5});
		multi::array_cref<int, 2>            ref2(&darr[0][0], {4, 5});

		BOOST_TEST( &ref[1][2] == &cref [1][2] );
		BOOST_TEST( &ref[1][2] == &crefc[1][2] );
		BOOST_TEST( &ref[1][2] == & ref2[1][2] );

		ref[1][1] = 20;
		BOOST_TEST( ref[1][1] == 20 );  // cppcheck-suppress knownConditionTrueFalse;

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif

		// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
		int darr2[4][5] = {
			{10, 00},
			{20, 30},
		};

		darr2[1][0] = 20;
		BOOST_TEST( darr2[1][0] == 20 );  // cppcheck-suppress knownConditionTrueFalse;

		// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
		auto const& dd = std::as_const(darr2);

		BOOST_TEST( &(dd[1][2]) == &(darr2[1][2]) );
		BOOST_TEST(( & ref[1].static_array_cast<int, int const*>()[1] == &ref[1][1] ));
		BOOST_TEST(( &multi::static_array_cast<int, int const*>(ref[1])[1] == &ref[1][1] ));

#ifdef __clang__
#pragma clang diagnostic pop
#endif
	}

	// BOOST_AUTO_TEST_CASE(array_ref_cast_carray)
	{
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif

		// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
		double darr[2][2] = {
			{1.0, 2.0},
			{2.0, 3.0},
		};
		multi::array_ref<double, 2> ref(&darr[0][0], {2, 2});

		// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
		auto&& other_darr = static_cast<double (&)[2][2]>(ref);

		// NOLINTNEXTLINE(hicpp-use-auto,modernize-use-auto,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
		double (&other_darr2)[2][2] = static_cast<double (&)[2][2]>(ref);
		double (&other_darr3)[2][2](ref);  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type

		BOOST_TEST( &ref        [1][0] == &darr[1][0] );
		BOOST_TEST( &other_darr [1][0] == &darr[1][0] );
		BOOST_TEST( &other_darr2[1][0] == &darr[1][0] );
		BOOST_TEST( &other_darr3[1][0] == &darr[1][0] );

#ifdef __clang__
#pragma clang diagnostic pop
#endif

		// TODO(correaa) adapt this test to Boost Test Lightweight
		// Homebrew GCC-13 terminates rather than having the expected exception caught.
		// MSVC 17.10 fails to compile
		// #if !(defined(__GNUC__) && __GNUC__ >= 5 && defined(__APPLE__)) && !defined(_MSC_VER)
		//  BOOST_REQUIRE_THROW(
		//      ([&] {
		//          double(&other_darr4)[3][3](ref);  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
		//          other_darr4[1][1] += 1.0;
		//      }()),
		//      std::bad_cast
		//  );
		// #endif
	}

	// BOOST_AUTO_TEST_CASE(array_ref_original_tests_const_carray)
	{
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif

		// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
		double const d2D[4][5] = {
			{1.0, 2.0},
			{2.0, 3.0},
		};
		multi::array_ref<double, 2, double const*> d2Rce(&d2D[0][0], {4, 5});

		BOOST_TEST( &d2Rce[2][3] == &d2D[2][3] );
		BOOST_TEST( d2Rce.size() == 4 );
		BOOST_TEST( d2Rce.num_elements() == 20 );

#ifdef __clang__
#pragma clang diagnostic pop
#endif
	}

	// BOOST_AUTO_TEST_CASE(array_ref_original_tests_const_carray_string)
	{
		using std::get;  // workaround no prior declaration in function call with explicit template arguments is a C++20 extension [-Wc++20-extensions]

		// NOLINTBEGIN(fuchsia-default-arguments-calls) std::string ctor
		// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
		std::string const dc3D[4][2][3] = {
			{{"A0a", "A0b", "A0c"}, {"A1a", "A1b", "A1c"}},
			{{"B0a", "B0b", "B0c"}, {"B1a", "B1b", "B1c"}},
			{{"C0a", "C0b", "C0c"}, {"C1a", "C1b", "C1c"}},
			{{"D0a", "D0b", "D0c"}, {"D1a", "D1b", "D1c"}},
		};
		// NOLINTEND(fuchsia-default-arguments-calls) std::string ctor

		multi::array_cref<std::string, 3> cref(&dc3D[0][0][0], {4, 2, 3});
		BOOST_TEST( cref.num_elements() == 24 && cref[2][1][1] == "C1b" );

		auto const& A2 = cref.sliced(0, 3).rotated()[1].sliced(0, 2).unrotated();
		BOOST_TEST( std::decay_t<decltype(A2)>::dimensionality == 2 && A2.num_elements() == 6 );
		// BOOST_TEST( multi::rank<std::decay_t<decltype(A2)>>{} == 2 && A2.num_elements() == 6 );

		BOOST_TEST( get<0>(sizes(A2)) == 3 && get<1>(sizes(A2)) == 2 );

		auto const& A3 = cref({0, 3}, 1, {0, 2});
		BOOST_TEST( std::decay_t<decltype(A3)>::dimensionality == 2 && A3.num_elements() == 6 );

		BOOST_TEST( A2.layout()[2][1] == &A2[2][1] - A2.base() );
		BOOST_TEST( A2.rotated().layout()[1][2] == &A2.rotated()[1][2] - A2.rotated().base() );
	}

	// BOOST_AUTO_TEST_CASE(array_ref_sizes_assingment)
	{
		multi::array_cref<int, 3> const cref(nullptr, {4, 2, 3});
		{
			auto [sizes1, sizes2, sizes3] = cref.sizes();
			BOOST_TEST( sizes1 == 4 );
			BOOST_TEST( sizes2 == 2 );
			BOOST_TEST( sizes3 == 3 );
		}
		{
			using std::get;  // workaround no prior declaration in function call with explicit template arguments is a C++20 extension [-Wc++20-extensions]

			auto sizes1 = get<0>(cref.sizes());
			auto sizes2 = get<1>(cref.sizes());
			auto sizes3 = get<2>(cref.sizes());
			BOOST_TEST( sizes1 == 4 );
			BOOST_TEST( sizes2 == 2 );
			BOOST_TEST( sizes3 == 3 );
		}
		// {
		// 	multi::size_t sizes1;  // NOLINT(cppcoreguidelines-init-variables)
		// 	multi::size_t sizes2;  // NOLINT(cppcoreguidelines-init-variables)
		// 	multi::size_t sizes3;  // NOLINT(cppcoreguidelines-init-variables)
		// 	std::tie(sizes1, sizes2, sizes3) = cref.sizes();

		// 	BOOST_TEST( sizes1 == 4 );
		// 	BOOST_TEST( sizes2 == 2 );
		// 	BOOST_TEST( sizes3 == 3 );
		// }
		{
			multi::ssize_t sizes1;  // NOLINT(cppcoreguidelines-init-variables)
			multi::ssize_t sizes2;  // NOLINT(cppcoreguidelines-init-variables)
			multi::ssize_t sizes3;  // NOLINT(cppcoreguidelines-init-variables)
			multi::tie(sizes1, sizes2, sizes3) = cref.sizes();

			BOOST_TEST( sizes1 == 4 );
			BOOST_TEST( sizes2 == 2 );
			BOOST_TEST( sizes3 == 3 );
		}
		{
			auto const [sizes1, sizes2, sizes3] = cref.sizes();

			BOOST_TEST( sizes1 == 4 );
			BOOST_TEST( sizes2 == 2 );
			BOOST_TEST( sizes3 == 3 );
		}
		// {
		//  // NOLINTNEXTLINE(runtime/int)
		//  long sizes1;  // NOLINT(google-runtime-int,cppcoreguidelines-init-variables) test bad idiom
		//  // NOLINTNEXTLINE(runtime/int)
		//  long sizes2;  // NOLINT(google-runtime-int,cppcoreguidelines-init-variables) test bad idiom
		//  // NOLINTNEXTLINE(runtime/int)
		//  long sizes3;  // NOLINT(google-runtime-int,cppcoreguidelines-init-variables) test bad idiom

		//  multi::tie(sizes1, sizes2, sizes3) = static_cast<multi::tuple<long, long, long>>(cref.sizes());

		//  BOOST_TEST( sizes1 == 4L );
		//  BOOST_TEST( sizes2 == 2L );
		//  BOOST_TEST( sizes3 == 3L );
		// }
		{
			// NOLINTNEXTLINE(runtime/int)
			long long sizes1;  // NOLINT(google-runtime-int,cppcoreguidelines-init-variables) test bad idiom
			// NOLINTNEXTLINE(runtime/int)
			long long sizes2;  // NOLINT(google-runtime-int,cppcoreguidelines-init-variables) test bad idiom
			// NOLINTNEXTLINE(runtime/int)
			long long sizes3;  // NOLINT(google-runtime-int,cppcoreguidelines-init-variables) test bad idiom
			// std::tie(sizes1, sizes2, sizes3) = cref.sizes();
			multi::tie(sizes1, sizes2, sizes3) = cref.sizes();

			BOOST_TEST( sizes1 == 4 );
			BOOST_TEST( sizes2 == 2 );
			BOOST_TEST( sizes3 == 3 );
		}
		{
			int64_t sizes1;  // NOLINT(cppcoreguidelines-init-variables)
			int64_t sizes2;  // NOLINT(cppcoreguidelines-init-variables)
			int64_t sizes3;  // NOLINT(cppcoreguidelines-init-variables)
			// std::tie(sizes1, sizes2, sizes3) = cref.sizes();
			multi::tie(sizes1, sizes2, sizes3) = cref.sizes();

			BOOST_TEST( sizes1 == 4 );
			BOOST_TEST( sizes2 == 2 );
			BOOST_TEST( sizes3 == 3 );
		}
	}

	// BOOST_AUTO_TEST_CASE(array_ref_rebuild_2D) {
	//  // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
	//  int d2D[4][5] = {
	//      {10, 20},
	//      {20, 30},
	//  };
	//  multi::array_ref<int, 2> d2R(&d2D[0][0], {4, 5});

	//  auto&& d2B     = d2R();
	//  auto&& d2B_ref = multi::ref(d2B.begin(), d2B.end());

	//  BOOST_TEST(  d2B[0][0]    ==  d2B_ref[0][0] );
	//  BOOST_TEST( &d2B[0][0]    == &d2B_ref[0][0] );

	//  BOOST_TEST(  d2B.base()   ==  d2B_ref.base() );
	//  BOOST_TEST(  d2B.layout() ==  d2B_ref.layout() );

	//  BOOST_TEST( &d2R() == &multi::ref(d2B.begin(), d2B.end()) );
	// }

	// BOOST_AUTO_TEST_CASE(array_ref_rebuild_1D) {
	// double d1D[5] = {1.0, 2.0, 3.0, 4.0, 5.0};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type

	// multi::array_ref<double, 1> d1R(&d1D[0], {5});

	// auto&& d1B     = d1R();
	// auto&& d1B_ref = multi::ref(d1B.begin(), d1B.end());

	// BOOST_TEST( d1B.base()   == d1B_ref.base() );
	// BOOST_TEST( d1B.layout() == d1B_ref.layout() );
	// BOOST_TEST( &d1R() == &multi::ref(d1B.begin(), d1B.end()) );
	//}

	// BOOST_AUTO_TEST_CASE(array_ref_move_assigment_2D)
	{
		{
			multi::array<double, 2> arr({5, 4});
			std::iota(arr.elements().begin(), arr.elements().end(), 0.);

			multi::array<double, 2> arr2({5, 4});
			std::iota(arr2.elements().begin(), arr2.elements().end(), 10.);

			auto&& Aref = multi::array_ref<double, 2>({5, 4}, arr.data_elements());
			auto&& Bref = multi::array_ref<double, 2>({5, 4}, arr2.data_elements());

			Bref = Aref;

			BOOST_TEST( arr2 == arr );
		}
		{
			multi::array<double, 2> arr({5, 4});
			std::iota(arr.elements().begin(), arr.elements().end(), 0.0);

			multi::array<double, 2> arr2({5, 4});
			std::iota(arr2.elements().begin(), arr2.elements().end(), 10.0);

			auto&& ref2 = multi::array_ref<double, 2>({5, 4}, arr2.data_elements());

			ref2 = multi::array_ref<double, 2>({5, 4}, arr.data_elements());

			BOOST_TEST( arr2 == arr );
		}
		{
			multi::array<double, 2> arr({5, 4});
			std::iota(arr.elements().begin(), arr.elements().end(), 0.0);

			multi::array<double, 2> arr2({5, 4});
			std::iota(arr2.elements().begin(), arr2.elements().end(), 10.0);

			auto&& ref  = multi::array_ref<double, 2>({5, 4}, arr.data_elements());
			auto&& ref2 = multi::array_ref<double, 2>({5, 4}, arr2.data_elements());

			ref2 = std::move(ref);

			BOOST_TEST( arr2 == arr );
		}
	}

	// BOOST_AUTO_TEST_CASE(array_ref_conversion_1D)
	{
		multi::array<int, 1> arr({5}, int{});
		BOOST_TEST( arr.size() == 5 );
		std::iota(arr.elements().begin(), arr.elements().end(), 0);

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif

		{
			auto const& carr = static_cast<int (&)[5]>(arr);  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
			BOOST_TEST( &carr[3] == &arr[3] );

			BOOST_TEST(f1d5(static_cast<int(&)[5]>(arr)) == 1 );  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
		}
		{
			int (&carr)[5](arr);  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
			BOOST_TEST( &carr[3] == &arr[3] );
		}

#ifdef __clang__
#pragma clang diagnostic pop
#endif
	}

	// BOOST_AUTO_TEST_CASE(array_ref_conversion_2D)
	{
		multi::array<int, 2> arr({5, 4});
		std::iota(arr.elements().begin(), arr.elements().end(), 0);

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif

		{
			auto const& carr = static_cast<int (&)[5][4]>(arr);  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
			BOOST_TEST( &carr[3][2] == &arr[3][2] );

			f2d54(static_cast<int (&)[5][4]>(arr));  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
		}
		{
			int (&carr)[5][4](arr);  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
			BOOST_TEST( &carr[3][2] == &arr[3][2] );
		}

#ifdef __clang__
#pragma clang diagnostic pop
#endif
	}

#ifndef _MSC_VER
	// BOOST_AUTO_TEST_CASE(as_span)
	{
#ifdef BOOST_MULTI_HAS_SPAN
		auto print_me0 = [](std::span<int> rng) {
			std::cout << "rng.size(): " << rng.size() << '\n';                                          // (4)
			std::for_each(rng.begin(), rng.end(), [](auto const& elem) { std::cout << elem << ' '; });  // NOLINT(modernize-use-ranges)
			std::cout << "\n\n";
		};
#endif

		auto print_me1 = [](multi::array_ref<int, 1> const& rng) -> void {
			std::cout << "rng.size(): " << rng.size() << '\n';                                                  // (4)
			std::for_each(rng.begin(), rng.end(), [](auto const& elem) -> void { std::cout << elem << ' '; });  // NOLINT(llvm-use-ranges,modernize-use-ranges)
			std::cout << "\n\n";
		};

		auto print_me2 = [](auto const ptr) {
			std::cout << "ptr->size(): " << ptr->size() << '\n';  // (4)
			std::for_each(ptr->begin(), ptr->end(), [](auto const& elem) { std::cout << elem << ' '; });
			std::cout << "\n\n";
		};

#ifdef BOOST_MULTI_HAS_SPAN
		{
			int arr[] = {1, 2, 3, 4};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy arrays
			print_me0(arr);

			// vvv this fails in certain versions of clang (14?ss)
			// std::vector vec = {1, 2, 3, 4, 5};  // NOLINT(fuchsia-default-arguments-calls)
			// print_me0(vec);

			// clang-format off
		std::array<int, 6> arr2 = {{1, 2, 3, 4, 5, 6}};
			// clang-format on

			print_me0(arr2);
		}
#endif
		{
			int arr[] = {1, 2, 3, 4};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test c-arrays

			print_me1(multi::array_ref<int, 1>{arr});
			print_me1(arr);

			std::vector<int> vec = {1, 2, 3, 4, 5};  // NOLINT(fuchsia-default-arguments-calls)

			print_me1(*&multi::array_ref<int, 1>({5}, vec.data()));

			// clang-format off
			std::array<int, 6> arr2 = {{1, 2, 3, 4, 5, 6}};
			// clang-format on

			print_me1(arr2);
			print_me1(*&multi::array_ref<int, 1>({6}, arr2.data()));

			multi::dynamic_array<int, 1> marr(
				// #ifdef _MSC_VER  // problems with MSVC 14.3 c++17
				multi::extents_t<1>
				// #endif
				{10},
				99
			);

			print_me1(*&multi::array_ref<int, 1>(10, marr.data_elements()));

			auto const& alias = marr;

			marr = alias;
			BOOST_TEST(marr[5] == 99);

			marr = alias();
			BOOST_TEST(marr[5] == 99);
		}
#if defined(BOOST_MULTI_HAS_SPAN) && !defined(__NVCC__)
#ifdef __cpp_lib_span
#if !defined(__clang__) || (__clang_major__ > 14)
#if defined(__clang__) || !defined(__GNUC__) || __GNUC__ > 11
		{
			std::vector<int> vec = {1, 2, 3};

			multi::array_ref<int, 1> aref(std::span{vec});
			BOOST_TEST( &aref[1] == &vec[1] );
		}
#endif
		{
			std::vector<int> vec = {1, 2, 3};

			multi::array_ref<int, 1> aref(std::span<int>{vec});
			BOOST_TEST( &aref[1] == &vec[1] );
		}
#endif
#endif
#endif
		{
			int arr[] = {1, 2, 3, 4};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test c-arrays
			print_me2(&multi::array_ref<int, 1>{arr});

			std::vector<int> vec = {1, 2, 3, 4, 5};  // NOLINT(fuchsia-default-arguments-calls)
			print_me2(&multi::array_ref<int, 1>(5, vec.data()));

			// clang-format off
			std::array<int, 6> arr2 = {{1, 2, 3, 4, 5, 6}};
			// clang-format on

			//  print_me2(&arr2);  // this crashes clang-tidy
			print_me2(&multi::array_ref<int, 1>({6}, arr2.data()));

			//  multi::dynamic_array<int, 1> marr({10}, 99);
			//  print_me2(&marr);  // TODO(correaa) make this work
		}
	}
#endif

	// BOOST_AUTO_TEST_CASE(diagonal)
	{
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif

		// NOLINTNEXTLINE(hicpp-avoid-c-arrays, modernize-avoid-c-arrays, cppcoreguidelines-avoid-c-arrays): test
		int arr[4][3] = {
			{  0,  10,  20},
			{ 50,  10,  70},
			{100, 110,  20},
			{990, 990, 999},
		};

		// NOLINTNEXTLINE(hicpp-avoid-c-arrays, modernize-avoid-c-arrays, cppcoreguidelines-avoid-c-arrays): special type
		multi::array_ref<int, 2> mar = *&multi::array_ref<int, 2>(arr);

		BOOST_TEST( &mar({0, 3}, {0, 3}).diagonal()[0] == &arr[0][0] );
		BOOST_TEST( &mar({0, 3}, {0, 3}).diagonal()[1] == &arr[1][1] );
		BOOST_TEST( &mar({0, 3}, {0, 3}).diagonal()[2] == &arr[2][2] );

#ifdef __clang__
#pragma clang diagnostic pop
#endif

		auto sum = 0;

		// NOLINTNEXTLINE(altera-unroll-loops) test for-range loop
		for(auto const& aii : mar.diagonal()) {
			sum += aii;  // cppcheck-suppress useStlAlgorithm;
		}
		BOOST_TEST( sum == mar[0][0] + mar[1][1] + mar[2][2]);
	}

	// BOOST_AUTO_TEST_CASE(function_passing)
	{
		multi::array<int, 2>      arr({3, 3});
		multi::array_ref<int, 2>& arrR = arr;

		arrR[0][0] = 21;

		arr.reextent({5, 5});

		BOOST_TEST( arr [0][0] == 21 );  // cppcheck-suppress knownConditionTrueFalse;
		BOOST_TEST(&arrR[0][0] == &arr[0][0]);
	}

	// BOOST_AUTO_TEST_CASE(function_passing_2)
	{
		multi::Array<int, 2>                   arr({3, 3});
		[[maybe_unused]] multi::Array<int&, 2> arrR = arr;

		arrR[0][0] = 51;

		[[maybe_unused]] multi::Array<int const&, 2> arrCR = arr;

		BOOST_TEST(&arrCR[0][0] == &arrR[0][0]);

		[[maybe_unused]] multi::Array<int const&, 2> arrCR2 = arrCR;

		arr.reextent({5, 5});
		BOOST_TEST( arrR[0][0] == 51 );  // cppcheck-suppress knownConditionTrueFalse;
		BOOST_TEST(&arrR[0][0] == &arr[0][0]);
	}

	// BOOST_AUTO_TEST_CASE(function_passing_3)
	{
		multi::array<int, 2> const arr({3, 3}, 10);

		BOOST_TEST( trace_array_deduce     (arr) == 30 );
		BOOST_TEST( trace_array_deduce<int>(arr) == 30 );

		multi::array<int, 2> const arr_paren_copy{arr()};
		BOOST_TEST( arr_paren_copy.size() == 3 );

		BOOST_TEST(( trace_separate_ref                         (arr) == 30 ));
		BOOST_TEST(( trace_separate_sub                         (arr) == 30 ));
	}

#if __cplusplus > 202002L || (defined(_MSVC_LANG) && _MSVC_LANG > 202002L)
	// BOOST_AUTO_TEST_CASE(function_passing_3_lambdas)
	{
		auto buffer = std::make_unique<int[]>(9);
		std::fill_n(buffer.get(), 9, 1);

		multi::array<int, 2> const     arr({3, 3}, 1);
		multi::array_ref<int, 2> const aref(buffer.get(), {3, 3});

		auto const& asub = arr({0, 3}, {0, 3});

		auto deduce_array = []<class Arr>(Arr const& a) {
			return std::accumulate(a.diagonal().begin(), a.diagonal().end(), typename Arr::element{0});
		};  // NOLINT(readability/braces) bug in cpplint 1.6.1

		BOOST_TEST( deduce_array(arr ) == 3 );
		BOOST_TEST( deduce_array(aref) == 3 );
		BOOST_TEST( deduce_array(asub) == 3 );

		auto deduce_element = []<class T>(multi::array<T, 2> const& a) {
			return std::accumulate(a.diagonal().begin(), a.diagonal().end(), T{0});
		};  // NOLINT(readability/braces) bug in cpplint 1.6.1

		BOOST_TEST( deduce_element(arr) == 3 );
		// BOOST_TEST( deduce_element(aref) == 30 );
		// BOOST_TEST( deduce_element(asub) == 30 );

		auto deduce_element_ref = []<class T>(multi::array_ref<T, 2> const& a) {
			return std::accumulate(a.diagonal().begin(), a.diagonal().end(), T{0});
		};  // NOLINT(readability/braces) cpplint 1.6.1 gets confused

		BOOST_TEST( deduce_element_ref(arr) == 3 );
		BOOST_TEST( deduce_element_ref(aref) == 3 );
		// BOOST_TEST( deduce_element_ref(asub) == 3 );

		auto deduce_element_sub = []<class T, class Ptr>(multi::const_subarray<T, 2, Ptr> const& a) {
			return std::accumulate(a.diagonal().begin(), a.diagonal().end(), T{0});
		};  // NOLINT(readability/braces) cpplint 1.6.1 gets confused

		BOOST_TEST( deduce_element_sub(arr) == 3 );
		BOOST_TEST( deduce_element_sub(aref) == 3 );
		BOOST_TEST( deduce_element_sub(asub) == 3 );
	}
#endif

	// BOOST_AUTO_TEST_CASE(function_passing_4) {
	//  multi::array<int, 2> arr({3, 3}, 10);

	//  BOOST_TEST( mut_trace_array_deduce     (arr) == 30 );
	//  BOOST_TEST( mut_trace_array_deduce<int>(arr) == 30 );

	//  BOOST_TEST(  mut_trace_generic                       (arr) == 30  );
	//  BOOST_TEST(( mut_trace_generic<multi::array<int, 2> >(arr) == 30 ));
	// }

	// BOOST_AUTO_TEST_CASE(array_fill_constructor)
	{
		multi::array<int, 2> arr(3, multi::array<int, 1>{10, 20, 30, 40});

		BOOST_TEST( arr[0][1] == 20 );
		BOOST_TEST( arr[1][1] == 20 );
	}

	// BOOST_AUTO_TEST_CASE(array_fill_constructor_1D)
	{
		using std::get;  // workaround no prior declaration in function call with explicit template arguments is a C++20 extension [-Wc++20-extensions]

		multi::array<int, 1> arr(3, 10);

		BOOST_TEST( arr[0] == 10 );
		BOOST_TEST( arr[1] == 10 );

		BOOST_TEST( get<0>(arr().sizes()) ==  3 );
		BOOST_TEST( get<0>(sizes(arr())) ==  3 );

		// BOOST_TEST( get<1>(sizes(arr())) == 10 );

		// BOOST_TEST( get<0>(sizes(arr)) ==  3 );
		// BOOST_TEST( get<1>(sizes(arr)) == 10 );
	}

	// BOOST_AUTO_TEST_CASE(array_fill_constructor_2D)
	{
		using std::get;  // workaround no prior declaration in function call with explicit template arguments is a C++20 extension [-Wc++20-extensions]

		multi::array<int, 2> arr({3, 4}, 10);

		BOOST_TEST( get<0>(arr().sizes()) ==  3 );
		BOOST_TEST( get<0>(sizes(arr())) ==  3 );

		BOOST_TEST( get<1>(arr().sizes()) ==  4 );
		BOOST_TEST( get<1>(sizes(arr())) ==  4 );

		// BOOST_TEST( get<1>(sizes(arr())) == 10 );

		// BOOST_TEST( get<0>(sizes(arr)) ==  3 );
		// BOOST_TEST( get<1>(sizes(arr)) == 10 );
	}

	return boost::report_errors();
}  // NOLINT(readability/fn_size)
