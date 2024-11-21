// Copyright 2019-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>
#include <boost/multi/detail/static_allocator.hpp>  // TODO(correaa) export IWYU

#include <algorithm>  // for transform, is_sorted
#include <array>      // for array, operator==
#include <cassert>    // for _LIBCPP_VERSION  // IWYU pragma: keep
#include <cstddef>    // for __GLIBCXX__, size_t
#include <iterator>   // for size, back_insert...
#include <memory>     // for make_unique, uniq...
#ifdef BOOST_MULTI_HAS_MEMORY_RESOURCE
	#include <memory_resource>  // for monotonic_buffer_...
#endif
#include <new>      // for operator new  // NOLINT(misc-include-cleaner)
#include <string>   // for basic_string, string
#include <utility>  // for move, forward
#include <vector>   // for vector, allocator

namespace multi = boost::multi;

template<class T, multi::dimensionality_type D, std::size_t Capacity = 4UL * 4UL>
using small_array = multi::static_array<T, D, multi::detail::static_allocator<T, Capacity>>;
// https://godbolt.org/z/d8ozWahna

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
BOOST_AUTO_TEST_CASE(static_array_allocator) {
	multi::array<int, 2> const                            ma({2, 3}, 99);
	multi::static_array<int, 2, std::allocator<int>> const sma(ma(), std::allocator<int>{});
	BOOST_TEST( sma == ma );
}

BOOST_AUTO_TEST_CASE(empty_stride) {
	multi::array<double, 2> ma;
	BOOST_TEST(ma.size() == 0);
	BOOST_TEST(ma.stride() != 0);
	BOOST_TEST(size(ma) == 0);

	multi::array<double, 2> ma0({0, 0}, 0.0);
	BOOST_TEST(ma0.size() == 0);
	BOOST_TEST(ma0.stride() != 0);
#ifndef _MSC_VER  // doesn't work with msvc 14.3 c++17 permissive mode
	BOOST_TEST(size(ma0) == 0);
#endif
}

BOOST_AUTO_TEST_CASE(std_vector_of_arrays_check_size) {
	multi::array<int, 2> const ma;
	BOOST_TEST( ma.size() == 0 );
	BOOST_TEST( ma.num_elements() == 0 );

	std::vector<multi::array<int, 2>> va(1);  // NOLINT(fuchsia-default-arguments-calls) vector

	BOOST_TEST( va[0].size() == 0 );
}

BOOST_AUTO_TEST_CASE(std_vector_of_arrays_manual_emplaceback_ctor) {
	std::vector<multi::array<int, 2>> va;

	va.emplace_back(multi::extensions_t<2>{3, 3}, 3);
	va.emplace_back(multi::extensions_t<2>{2, 2}, 2);
	va.emplace_back(multi::extensions_t<2>{1, 1}, 1);
	va.emplace_back(multi::extensions_t<2>{0, 0}, 0);
}

BOOST_AUTO_TEST_CASE(std_vector_of_arrays_manual_emplaceback) {
	std::vector<multi::array<int, 2>> va;

	va.emplace_back(multi::array<int, 2>({2, 2}, 2));
	va.emplace_back(multi::array<int, 2>({1, 1}, 1));
	va.emplace_back(multi::array<int, 2>({0, 0}, 0));
}

BOOST_AUTO_TEST_CASE(std_vector_of_arrays_manual_pushback) {
	std::vector<multi::array<int, 2>> va;

	va.push_back(multi::array<int, 2>({2, 2}, 2));
	va.push_back(multi::array<int, 2>({1, 1}, 1));
	va.push_back(multi::array<int, 2>({0, 0}, 0));
}

BOOST_AUTO_TEST_CASE(std_vector_of_arrays) {
	std::vector<multi::array<int, 2>> va;

	std::transform(
		multi::iextension(3).begin(), multi::iextension(3).end(),
		std::back_inserter(va),
		[](auto idx) { return multi::array<int, 2>({idx, idx}, static_cast<int>(idx)); }
	);

#ifndef _MSC_VER  // doesn't work with msvc 14.3 c++17 permissive mode
	BOOST_TEST( size(va[0]) == 0 );
	BOOST_TEST( size(va[1]) == 1 );
	BOOST_TEST( size(va[2]) == 2 );
#endif

	BOOST_TEST( va[1] [0][0] == 1 );
	BOOST_TEST( va[2] [0][0] == 2 );

	using namespace std::string_literals;  // NOLINT(build/namespaces)

#ifndef _MSC_VER  // doesn't work with msvc 14.3 c++17 permissive mode
	// NOLINTNEXTLINE(fuchsia-default-arguments-calls)
	std::vector<multi::array<int, 2>> const wa = {
		multi::array<int, 2>({0, 0}, 0),
		multi::array<int, 2>({1, 1}, 1),
		multi::array<int, 2>({2, 2}, 2),
	};
#else
	// NOLINTNEXTLINE(fuchsia-default-arguments-calls)
	std::vector<multi::array<int, 2>> const wa = {
		multi::array<int, 2>(multi::extensions_t<2>(0, 0), 0),
		multi::array<int, 2>(multi::extensions_t<2>(1, 1), 1),
		multi::array<int, 2>(multi::extensions_t<2>(2, 2), 2),
	};
#endif

	BOOST_TEST( va.size() == wa.size() );
	BOOST_TEST( va == wa );

	std::vector<multi::array<int, 2>> ua(3, std::allocator<multi::array<double, 2>>{});

	auto iex = multi::iextension(static_cast<multi::size_type>(ua.size()));

	std::transform(
		begin(iex), end(iex),
		begin(ua),
		[](auto idx) { return multi::array<int, 2>({idx, idx}, static_cast<int>(idx)); }
	);
	BOOST_TEST( ua == va );
}

BOOST_AUTO_TEST_CASE(std_vector_of_arrays_with_string_instead_of_int) {
	// NOLINTBEGIN(fuchsia-default-arguments-calls)  // string uses default parameter
	std::vector<multi::array<std::string, 2>> va;
	std::transform(
		begin(multi::iextension(3)), end(multi::iextension(3)),
		std::back_inserter(va),
		[](auto idx) { return multi::array<std::string, 2>({idx, idx}, std::to_string(idx)); }
	);

#ifndef _MSC_VER  // doesn't work with msvc 14.3 c++17 permissive mode
	BOOST_TEST( size(va[0]) == 0 );
	BOOST_TEST( size(va[1]) == 1 );
	BOOST_TEST( size(va[2]) == 2 );
#endif
	using namespace std::string_literals;  // NOLINT(build/namespaces)

	BOOST_TEST( va[1] [0][0] == "1"s );  // NOLINT(misc-include-cleaner) bug in clang-tidy 18
	BOOST_TEST( va[2] [0][0] == "2"s );

#ifndef _MSC_VER  // doesn't work with msvc 14.3 c++17 permissive mode
	std::vector<multi::array<std::string, 2>> const wa = {
		multi::array<std::string, 2>({0, 0}, "0"s),
		multi::array<std::string, 2>({1, 1}, "1"s),
		multi::array<std::string, 2>({2, 2}, "2"s),
	};
#else
	std::vector<multi::array<std::string, 2>> const wa = {
		multi::array<std::string, 2>(multi::extensions_t<2>(0, 0), "0"s),
		multi::array<std::string, 2>(multi::extensions_t<2>(1, 1), "1"s),
		multi::array<std::string, 2>(multi::extensions_t<2>(2, 2), "2"s),
	};
#endif

#ifndef _MSC_VER  // doesn't work with msvc 14.3 c++17 permissive mode
	BOOST_TEST( size(va) == size(wa) );
#endif
	BOOST_TEST( va == wa );

	std::vector<multi::array<std::string, 2>> ua(3, std::allocator<multi::array<double, 2>>{});

	auto iex = multi::iextension(static_cast<multi::size_type>(ua.size()));

	std::transform(
		begin(iex), end(iex),
		begin(ua),
		[](auto idx) { return multi::array<std::string, 2>({idx, idx}, std::to_string(idx)); }
	);

	BOOST_TEST( ua == va );

	// NOLINTEND(fuchsia-default-arguments-calls)  // string uses default parameter
}

// TODO(correaa) make this code work with nvcc compiler (non device function called from device host through adl uninitialized_fill)
#if !(defined(__NVCC__) || defined(__HIP_PLATFORM_NVIDIA__) || defined(__HIP_PLATFORM_AMD__) || defined(__HIPCC__))
BOOST_AUTO_TEST_CASE(array1d_of_arrays2d) {
	multi::array<multi::array<std::string, 2>, 1> arr(multi::extensions_t<1>(multi::iextension{10}), multi::array<std::string, 2>{});
	BOOST_TEST( size(arr) == 10 );

	std::transform(
		begin(extension(arr)), end(extension(arr)), begin(arr),
		[](auto idx) { return multi::array<std::string, 2>({idx, idx}, std::to_string(idx)); }
	);

	BOOST_TEST( size(arr[0]) == 0 );
	BOOST_TEST( size(arr[1]) == 1 );
	BOOST_TEST( size(arr[8]) == 8 );

	using namespace std::string_literals;  // NOLINT(build/namespaces)
	BOOST_TEST( arr[8][4][4] == "8"s );
}

BOOST_AUTO_TEST_CASE(array_3d_of_array_2d) {
	multi::array<multi::array<int, 3>, 2> AA({10, 20}, multi::array<int, 3>{});
	std::transform(AA.extension().begin(), AA.extension().end(), AA.begin(), AA.begin(), [](auto idx, auto&& row) -> decltype(row) {
		std::transform(row.extension().begin(), row.extension().end(), row.begin(), [idx](auto jdx) {
			return multi::array<int, 3>({idx + jdx, idx + jdx, idx + jdx}, 99);
		});
		return std::forward<decltype(row)>(row);
	});

	// BOOST_TEST( AA[9][19].size() == 9 + 19 );

	// BOOST_TEST( std::size(AA[9][19]) == 9 + 19 );  // doesn't work on nvhpc 22.11
	// BOOST_TEST( size(AA[9][19]) == 9 + 19 );

	// BOOST_TEST( AA[9][19][1][1][1] == 99 );
}

BOOST_AUTO_TEST_CASE(array_3d_of_array_2d_no_init) {
	multi::array<multi::array<int, 3>, 2> AA({10, 20});
	std::transform(extension(AA).begin(), extension(AA).end(), AA.begin(), AA.begin(), [](auto idx, auto&& row) -> decltype(row) {
		std::transform(extension(row).begin(), extension(row).end(), row.begin(), [idx](auto jdx) {
			return multi::array<int, 3>({idx + jdx, idx + jdx, idx + jdx}, 99);
		});
		return std::forward<decltype(row)>(row);
	});

	BOOST_TEST( AA[9][19].size() == 9 + 19 );
	// BOOST_TEST( std::size(AA[9][19]) == 9 + 19 );  // doesn't work on nvhpc 22.11
	BOOST_TEST( size(AA[9][19]) == 9 + 19 );

	BOOST_TEST( AA[9][19][1][1][1] == 99 );
}
#endif

BOOST_AUTO_TEST_CASE(const_elements) {
	auto ptr = std::make_unique<int const>(2);
	// ok, can't assign  //  *ptr = 3.0;
	BOOST_TEST( *ptr == 2 );
}

#ifdef BOOST_MULTI_HAS_MEMORY_RESOURCE
BOOST_AUTO_TEST_CASE(pmr) {
	std::array<char, 13> buffer = {
		{'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C'}
    };

	std::pmr::monotonic_buffer_resource pool{std::data(buffer), std::size(buffer)};

	multi::array<char, 2, std::pmr::polymorphic_allocator<char>> Aarr({2, 2}, 'x', &pool);
	Aarr[0][0] = 'x';
	Aarr[0][1] = 'y';
	Aarr[1][0] = 'z';
	Aarr[1][1] = '&';

	multi::array<char, 2, std::pmr::polymorphic_allocator<char>> Barr({3, 2}, 'o', &pool);

	BOOST_TEST(( buffer != std::array<char, 13>{{'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C'}} ));

	#if defined(__GLIBCXX__)
	BOOST_TEST(( buffer == std::array<char, 13>{{'x', 'y', 'z', '&', 'o', 'o', 'o', 'o', 'o', 'o', 'A', 'B', 'C'}} ));
	#endif
	#if defined(_LIBCPP_VERSION)
	BOOST_TEST(( buffer == std::array<char, 13>{{'0', '1', '2', 'o', 'o', 'o', 'o', 'o', 'o', 'x', 'y', 'z', '&'}} ));
	#endif

	BOOST_TEST(Aarr[0][0] == 'x');
	BOOST_TEST(Barr[0][0] == 'o');
}

BOOST_AUTO_TEST_CASE(pmr2) {
	// clang-format off
	std::array<char, 13> buffer = {{'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X'}};
	// clang-format on

	std::pmr::monotonic_buffer_resource pool{std::data(buffer), std::size(buffer)};

	#ifndef _MSC_VER
	multi::pmr::array<char, 2> Aarr({2, 2}, 'a', &pool);
	multi::pmr::array<char, 2> Barr({3, 2}, 'b', &pool);
	#else
	multi::pmr::array<char, 2> Aarr(multi::extensions_t<2>{2, 2}, 'a', &pool);
	multi::pmr::array<char, 2> Barr(multi::extensions_t<2>{3, 2}, 'b', &pool);
	#endif

	#if defined(__GLIBCXX__)
	BOOST_TEST(( buffer == std::array<char, 13>{{'a', 'a', 'a', 'a', 'b', 'b', 'b', 'b', 'b', 'b', 'X', 'X', 'X'}} ));
	#endif
	#if defined(_LIBCPP_VERSION)
	BOOST_TEST(( buffer == std::array<char, 13>{{'X', 'X', 'X', 'b', 'b', 'b', 'b', 'b', 'b', 'a', 'a', 'a', 'a'}} ));
	#endif

	BOOST_TEST(Aarr[0][0] == 'a');
	BOOST_TEST(Barr[0][0] == 'b');
}

BOOST_AUTO_TEST_CASE(pmr_double_uninitialized) {
	std::array<int, 12> buffer{
		{4, 5, 6, 7, 8, 9, 10, 11, 996, 997, 998, 999}
    };

	std::pmr::monotonic_buffer_resource pool{static_cast<void*>(std::data(buffer)), 12 * sizeof(int)};

	multi::pmr::array<int, 2> Aarr({2, 2}, &pool);

	BOOST_TEST( buffer[0] == 4 );
	BOOST_TEST( buffer[1] == 5 );

	#if defined(__GLIBCXX__)
	BOOST_TEST(Aarr[0][0] == 4);
	#endif
	#if defined(_LIBCPP_VERSION)
	BOOST_TEST(Aarr[0][0] == 996);
	#endif
}
#endif

BOOST_AUTO_TEST_CASE(static_allocator) {
	using T = int;
	multi::detail::static_allocator<T, 32> sa{};

	auto* pp = sa.allocate(10);

	new (std::next(pp, 8)) T{42};

	BOOST_TEST( *std::next(pp, 8) == 42 );
	// (pp + 8)->~double();
	sa.deallocate(pp, 10);
}

#if defined(__cpp_constexpr) && (__cpp_constexpr > 202306L)
constexpr auto f() {
	std::vector<int> v = {1, 2, 3};
	return v.size();
}

BOOST_AUTO_TEST_CASE(constexpr_allocator_vector) {
	static_assert(f() == 3);
	BOOST_TEST( f() == 3 );
}

constexpr auto g() {
	multi::array<int, 2> arr = {
		{4, 5, 6},
		{1, 2, 3},
		{7, 8, 9},
	};
	std::sort(arr.begin(), arr.end());
	for(auto it = arr.diagonal().begin(); it != arr.diagonal().end(); ++it) {
		*it += 5;
	}
	auto ret = arr[1][1];
	return ret;
}

BOOST_AUTO_TEST_CASE(constexpr_allocator) {
	constexpr auto gg = g();
	static_assert(gg == 10);
	BOOST_TEST( gg == 10 );
}
#endif

#if !defined(_MSC_VER)  // static allocator does not work with MSVC implementation pf vector
BOOST_AUTO_TEST_CASE(static_allocator_on_vector_int) {
	std::vector<int, multi::detail::static_allocator<int, 32>> vv(10, 42);  // NOLINT(fuchsia-default-arguments-calls)
	BOOST_TEST( vv[3] == 42 );

	// auto ww = vv;
	// BOOST_TEST( ww[3] == 42 );

	// ww[3] = 51;
	// BOOST_TEST( ww[3] == 51 );
	// BOOST_TEST( vv[3] == 42 );

	// auto xx = std::move(ww);
	// BOOST_TEST( ww.empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved)
	// BOOST_TEST( vv[3] == 42 );
	// BOOST_TEST( xx[3] == 51 );

	// {
	//  std::vector<std::vector<int, multi::detail::static_allocator<int, 32>>> const VV = {vv, xx, vv};  // NOLINT(fuchsia-default-arguments-calls)
	//  BOOST_TEST( VV.size() == 3 );
	// }
}

BOOST_AUTO_TEST_CASE(static_allocator_on_vector_string) {
	std::string const cat = "catcatcatcatcatcatcatcatcatcatcatcatcatcatcatcatcatcatcatcatcatcatcatcat";  // NOLINT(fuchsia-default-arguments-calls)
	std::string const dog = "dogdogdogdogdogdogdogdogdogdogdogdogdogdogdogdogdogdogdogdogdogdogdogdog";  // NOLINT(fuchsia-default-arguments-calls)

	std::vector<std::string, multi::detail::static_allocator<std::string, 32>> vv(10, cat);  // NOLINT(fuchsia-default-arguments-calls)
	BOOST_TEST( vv[3] == cat );

	auto ww = vv;
	BOOST_TEST( ww[3] == cat );

	ww[3] = dog;
	BOOST_TEST( ww[3] == dog );
	BOOST_TEST( vv[3] == cat );

	auto xx = std::move(ww);
	BOOST_TEST( vv[3] == cat );
	BOOST_TEST( xx[3] == dog );
	BOOST_TEST( ww.empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved)

	// vv.resize(15);

	// swap(xx, vv);
	// BOOST_TEST( vv[3] == dog );
	// BOOST_TEST( xx[3] == cat );

	{
		std::vector<std::vector<std::string, multi::detail::static_allocator<std::string, 32>>> const VV = {vv, xx, vv};  // NOLINT(fuchsia-default-arguments-calls)
		BOOST_TEST( VV.size() == 3 );
		// swap(VV[0], VV[1]);
		// std::sort(VV.begin(), VV.end());
		// BOOST_TEST( std::is_sorted(VV.begin(), VV.end()) );
		// VV.resize(10, xx);
		// std::sort(VV.begin(), VV.end());
		// BOOST_TEST( std::is_sorted(VV.begin(), VV.end()) );
	}
}
#endif

#if !defined(_MSC_VER) || (_MSC_VER > 193030706)  // TODO(correaa) doesn't work on MSVC 14.3 in c++17 mode
BOOST_AUTO_TEST_CASE(small_array_int) {
	small_array<int, 2, 4UL * 4UL> vv({4, 4}, 42);

	BOOST_TEST( vv[3][3] == 42 );

	auto ww = vv;

	BOOST_TEST( ww[3][3] == 42 );
	BOOST_TEST( ww.base() != vv.base() );

	auto const* wwb = ww.base();
	auto const* vvb = vv.base();

	ww[3][3] = 51;

	BOOST_TEST( ww[3][3] == 51 );
	BOOST_TEST( vv[3][3] == 42 );

	swap(ww, vv);

	BOOST_TEST( vv[3][3] == 51 );
	BOOST_TEST( ww[3][3] == 42 );

	BOOST_TEST( ww.base() == wwb );
	BOOST_TEST( vv.base() == vvb );

	auto xx{std::move(ww)};

	BOOST_TEST( vv[3][3] == 51 );
	BOOST_TEST( xx[3][3] == 42 );
	// BOOST_TEST( ww[3][3] == 42 );
	BOOST_TEST( xx.base() != vv.base() );
	// BOOST_TEST( ww.empty() );

	small_array<int, 2, 4UL * 4UL> yy({4, 4});
	yy = vv;

	BOOST_TEST( yy == vv );

	yy = std::move(vv);
	BOOST_TEST( vv.size() == 4 );  // NOLINT(clang-analyzer-cplusplus.Move,bugprone-use-after-move,hicpp-invalid-access-moved)

	{
		std::vector<small_array<int, 2, 4UL * 4UL>> VV = {vv, xx, vv};  // NOLINT(fuchsia-default-arguments-calls)
		BOOST_TEST( VV.size() == 3 );
		swap(VV[0], VV[1]);
		std::sort(VV.begin(), VV.end());
		BOOST_TEST( std::is_sorted(VV.begin(), VV.end()) );
		VV.resize(10, xx);
		std::sort(VV.begin(), VV.end());
		BOOST_TEST( std::is_sorted(VV.begin(), VV.end()) );
	}
}
#endif

BOOST_AUTO_TEST_CASE(props_of_static_allocator) {
	{
		std::vector<int> vv(20, 11);  // NOLINT(fuchsia-default-arguments-calls)
		std::vector<int> ww = vv;
		BOOST_TEST( ww == vv );

		ww = vv;
		BOOST_TEST( ww == vv );

		ww = std::move(vv);
		BOOST_TEST( vv.size() == 0 );  // NOLINT(readability-container-size-empty,bugprone-use-after-move,hicpp-invalid-access-moved,clang-analyzer-cplusplus.Move)

		std::vector<int> xx(20, 22);  // NOLINT(fuchsia-default-arguments-calls)
		swap(ww, xx);
		BOOST_TEST( ww == std::vector<int>(20, 22) );  // NOLINT(fuchsia-default-arguments-calls)
	}
#if !defined(_MSC_VER)  // static_allocator doesn't work with MSVC implementation of vector
	{
		std::vector<int, multi::detail::static_allocator<int, 32>> vv(20, 11);  // NOLINT(fuchsia-default-arguments-calls)
		std::vector<int, multi::detail::static_allocator<int, 32>> ww = vv;
		BOOST_TEST( ww == vv );

		ww = vv;
		BOOST_TEST( ww == vv );

		ww = std::move(vv);
		BOOST_TEST( vv.size() == 0 );  // NOLINT(readability-container-size-empty,bugprone-use-after-move,hicpp-invalid-access-moved,clang-analyzer-cplusplus.Move)

		std::vector<int, multi::detail::static_allocator<int, 32>> xx(20, 22);  // NOLINT(fuchsia-default-arguments-calls)
		swap(ww, xx);
		BOOST_TEST(( ww == std::vector<int, multi::detail::static_allocator<int, 32>>(20, 22) ));  // NOLINT(fuchsia-default-arguments-calls)
	}
#endif
}

return boost::report_errors();}
