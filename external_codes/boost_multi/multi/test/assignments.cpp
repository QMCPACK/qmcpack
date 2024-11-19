// Copyright 2019-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// #if defined(__clang__)
//  #pragma clang diagnostic push
//  #pragma clang diagnostic ignored "-Wunknown-warning-option"
//  #pragma clang diagnostic ignored "-Wconversion"
//     #pragma clang diagnostic ignored "-Wextra-semi-stmt"
//  #pragma clang diagnostic ignored "-Wold-style-cast"
//  #pragma clang diagnostic ignored "-Wsign-conversion"
//     #pragma clang diagnostic ignored "-Wswitch-default"
//  #pragma clang diagnostic ignored "-Wundef"
// #elif defined(__GNUC__)
//  #pragma GCC diagnostic push
//  #if (__GNUC__ > 7)
//      #pragma GCC diagnostic ignored "-Wcast-function-type"
//  #endif
//  #pragma GCC diagnostic ignored "-Wconversion"
//  #pragma GCC diagnostic ignored "-Wold-style-cast"
//  #pragma GCC diagnostic ignored "-Wsign-conversion"
//  #pragma GCC diagnostic ignored "-Wundef"
// #elif defined(_MSC_VER)
//  #pragma warning(push)
//  #pragma warning(disable : 4244)
// #endif

// #ifndef BOOST_TEST_MODULE
//  #define BOOST_TEST_MAIN
// #endif

// #include <boost/test/included/unit_test.hpp>

// #if defined(__clang__)
//  #pragma clang diagnostic pop
// #elif defined(__GNUC__)
//  #pragma GCC diagnostic pop
// #elif defined(_MSC_VER)
//  #pragma warning(pop)
// #endif

#include <boost/multi/array.hpp>  // for array, apply, operator==, layout_t

#include <algorithm>  // for fill
#include <complex>    // for complex
#include <cstddef>    // for size_t
#include <iterator>   // for size
#include <memory>     // for std::allocator  // IWYU pragma: keep
// IWYU pragma: no_include <type_traits>  // for decay_t
#include <utility>  // for move
#include <vector>   // for vector, allocator

namespace multi = boost::multi;

namespace {

constexpr auto make_ref(int* ptr) {
	return multi::array_ref<int, 2>(ptr, { 5, 7 });
}

}  // namespace

template<class T, class Allocator>
auto eye(multi::extensions_t<2> exts, Allocator alloc) {
	multi::array<T, 2, Allocator> ret(exts, 0, alloc);
	std::fill(ret.diagonal().begin(), ret.diagonal().end(), 1);
	return ret;
}

template<class T>
auto eye(multi::extensions_t<2> exts) { return eye<T>(exts, std::allocator<T>{}); }

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE)

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
BOOST_AUTO_TEST_CASE(equality_1D) {
	multi::array<int, 1> arr  = { 10, 20, 30 };
	multi::array<int, 1> arr2 = { 10, 20, 30 };

	BOOST_TEST(    arr == arr2  );
	BOOST_TEST( !(arr != arr2) );

	BOOST_TEST(    arr() == arr2()  );
	BOOST_TEST( !(arr() != arr2()) );
}

BOOST_AUTO_TEST_CASE(equality_2D) {
	multi::array<int, 2> arr = {
		{ 10, 20, 30 },
		{ 40, 50, 60 },
	};
	multi::array<int, 2> arr2 = {
		{ 10, 20, 30 },
		{ 40, 50, 60 },
	};

	BOOST_TEST( arr == arr2 );
	BOOST_TEST( !(arr != arr2) );

	BOOST_TEST( arr() == arr2() );
	BOOST_TEST( !(arr() != arr2()) );

	BOOST_TEST( arr[0] == arr2[0] );
	BOOST_TEST( !(arr[0] != arr2[0]) );
}

BOOST_AUTO_TEST_CASE(multi_copy_move) {
	multi::array<double, 2> arr({ 3, 3 }, 0.0);
	multi::array<double, 2> arr2 = arr;
	BOOST_TEST( arr == arr2 );

	auto* arr_data = arr.data_elements();

	multi::array<double, 2> arr3 = std::move(arr);

	BOOST_TEST( arr3.data_elements() == arr_data );

	multi::array<double, 2> const arr4(std::move(arr2));
	BOOST_TEST( size(arr4) == 3 );
}

BOOST_AUTO_TEST_CASE(range_assignment) {
	{
		auto const ext = multi::make_extension_t(10L);

		multi::array<multi::size_t, 1> vec(ext.begin(), ext.end());

		BOOST_TEST( ext.size() == vec.size() );
		BOOST_TEST( vec[1] == 1L );
	}
	{
		multi::array<multi::size_t, 1> vec(multi::extensions_t<1>{ multi::iextension{ 10 } });

		auto const ext = extension(vec);

		vec.assign(ext.begin(), ext.end());
		BOOST_TEST( vec[1] == 1 );
	}
}

BOOST_AUTO_TEST_CASE(rearranged_assignment) {
	multi::array<int, 4> const tmp(
#ifdef _MSC_VER  // problem with 14.3 c++17
		multi::extensions_t<4>
#endif
		{ 14, 14, 7, 4 }
	);

	auto const ext5 = multi::extensions_t<5>{ 2, 14, 14, 7, 2 };

	[[maybe_unused]] auto const ext52 = ext5;

	[[maybe_unused]] multi::array<int, 5> const src_test(ext5);

	multi::array<int, 5> src(
#ifdef _MSC_VER  // problem with 14.3 c++17
		multi::extensions_t<5>
#endif
		{ 2, 14, 14, 7, 2 }
	);

	src[0][1][2][3][1] = 99;

	// BOOST_TEST( tmp.unrotated().partitioned(2).transposed().rotated().extensions() == src.extensions() );
	// BOOST_TEST( extensions(tmp.unrotated().partitioned(2).transposed().rotated()) == extensions(src) );
}

BOOST_AUTO_TEST_CASE(rearranged_assignment_resize) {
	multi::array<double, 2> const arrA({ 4, 5 });
	multi::array<double, 2>       arrB({ 2, 3 });

	arrB = arrA;
	BOOST_TEST( arrB.size() == 4 );
}

#ifndef _MSVER  // TODO(correaa) fix
// seems to produce a deterministic divide by zero
// Assertion failed: stride_ != 0, file D:\a\boost-multi\boost-root\boost/multi/detail/layout.hpp, line 767
// D:\a\boost-multi\boost-root\boost\multi\detail\layout.hpp(770) : error C2220: the following warning is treated as an error
// D:\a\boost-multi\boost-root\boost\multi\detail\layout.hpp(770) : warning C4723: potential divide by 0
// D:\a\boost-multi\boost-root\boost\multi\detail\layout.hpp(770) : warning C4723: potential divide by 0

BOOST_AUTO_TEST_CASE(rvalue_assignments) {
	using complex = std::complex<double>;

	std::vector<double> const vec1(200, 99.0);  // NOLINT(fuchsia-default-arguments-calls)
	std::vector<complex>      vec2(200);        // NOLINT(fuchsia-default-arguments-calls)

	auto linear1 = [&] { return multi::array_cptr<double, 1>(vec1.data(), 200); };
	auto linear2 = [&] { return multi::array_ptr<complex, 1>(vec2.data(), 200); };

	*linear2() = *linear1();
}
#endif

BOOST_AUTO_TEST_CASE(assignments) {
	{
		std::vector<int> vec(static_cast<std::size_t>(5 * 7), 99);  // NOLINT(fuchsia-default-arguments-calls)

		constexpr int val = 33;

		multi::array<int, 2> const arr({ 5, 7 }, val);
		multi::array_ref<int, 2>(vec.data(), { 5, 7 }) = arr();  // arr() is a subarray

		BOOST_TEST( vec[9] == val );
		BOOST_TEST( !vec.empty() );
		BOOST_TEST( !is_empty(arr) );
	}
	{
		std::vector<int> vec(5 * 7L, 99);  // NOLINT(fuchsia-default-arguments-calls)
		std::vector<int> wec(5 * 7L, 33);  // NOLINT(fuchsia-default-arguments-calls)

		multi::array_ptr<int, 2> const Bp(wec.data(), { 5, 7 });
		make_ref(vec.data()) = *Bp;

		auto&& mref = make_ref(vec.data());
		// mref        = (*Bp).sliced(0, 5);
		mref        = Bp->sliced(0, 5);

		// make_ref(vec.data()) = (*Bp).sliced(0, 5);
		make_ref(vec.data()) = Bp->sliced(0, 5);

		BOOST_TEST( vec[9] == 33 );
	}
	{
		std::vector<int> vec(5 * 7L, 99);  // NOLINT(fuchsia-default-arguments-calls)
		std::vector<int> wec(5 * 7L, 33);  // NOLINT(fuchsia-default-arguments-calls)

		make_ref(vec.data()) = make_ref(wec.data());

		BOOST_TEST( vec[9] == 33 );
	}
}

BOOST_AUTO_TEST_CASE(assigment_temporary) {
	multi::array<int, 2> Id = eye<int>(multi::extensions_t<2>({ 3, 3 }));
	BOOST_TEST( Id == eye<double>({3, 3}) );
	BOOST_TEST( Id[1][1] == 1 );
	BOOST_TEST( Id[1][0] == 0 );
}

return boost::report_errors();}
