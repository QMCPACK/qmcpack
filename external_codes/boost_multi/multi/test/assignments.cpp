// Copyright 2019-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <complex>

#include <boost/multi/array.hpp>

// Suppress warnings from boost.test
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wundef"
#pragma clang diagnostic ignored "-Wconversion"
#pragma clang diagnostic ignored "-Wsign-conversion"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wundef"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#elif defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4244)
#endif

#ifndef BOOST_TEST_MODULE
#define BOOST_TEST_MAIN
#endif

#include <boost/test/unit_test.hpp>

// #if defined(__clang__)
// #pragma clang diagnostic pop
// #elif defined(__GNUC__)
// #pragma GCC diagnostic pop
// #elif defined(_MSC_VER)
// #pragma warning(pop)
// #endif

namespace multi = boost::multi;

namespace {

constexpr auto make_ref(int* ptr) {
	return multi::array_ref<int, 2>(ptr, {5, 7});
}

}  // namespace

BOOST_AUTO_TEST_CASE(equality_1D) {
	multi::array<double, 1> arr  = {1.0, 2.0, 3.0};
	multi::array<double, 1> arr2 = {1.0, 2.0, 3.0};

	BOOST_REQUIRE(    arr == arr2  );
	BOOST_REQUIRE( ! (arr != arr2) );

	BOOST_REQUIRE(    arr() == arr2()  );
	BOOST_REQUIRE( ! (arr() != arr2()) );
}

BOOST_AUTO_TEST_CASE(equality_2D) {
	multi::array<double, 2> arr = {
		{1.0, 2.0, 3.0},
		{4.0, 5.0, 6.0},
	};
	multi::array<double, 2> arr2 = {
		{1.0, 2.0, 3.0},
		{4.0, 5.0, 6.0},
	};
	BOOST_REQUIRE( arr == arr2 );
	BOOST_REQUIRE( ! (arr != arr2) );

	BOOST_REQUIRE( arr() == arr2() );
	BOOST_REQUIRE( ! (arr() != arr2()) );

	BOOST_REQUIRE( arr[0] == arr2[0] );
	BOOST_REQUIRE( ! (arr[0] != arr2[0]) );
}

BOOST_AUTO_TEST_CASE(multi_copy_move) {
	multi::array<double, 2> arr({3, 3}, 0.0);
	multi::array<double, 2> arr2 = arr;
	BOOST_REQUIRE( arr == arr2 );

	auto* arr_data = arr.data_elements();

	multi::array<double, 2> arr3 = std::move(arr);

	BOOST_REQUIRE( arr3.data_elements() == arr_data );

	multi::array<double, 2> const arr4(std::move(arr2));
	BOOST_REQUIRE( size(arr4) == 3 );
}

BOOST_AUTO_TEST_CASE(range_assignment) {
	{
		auto const ext = multi::make_extension_t(10L);

		multi::array<int, 1> vec(ext.begin(), ext.end());

		BOOST_REQUIRE( ext.size() == vec.size() );
		BOOST_REQUIRE( vec[1] = 10 );
	}
	{
		multi::array<int, 1> vec(multi::extensions_t<1>{multi::iextension{10}});

		auto const ext = extension(vec);

		vec.assign(ext.begin(), ext.end());
		BOOST_REQUIRE( vec[1] == 1 );
	}
}

BOOST_AUTO_TEST_CASE(rearranged_assignment) {
	multi::array<int, 4> tmp(
#ifdef _MSC_VER  // problem with 14.3 c++17
		multi::extensions_t<4>
#endif
		{14, 14, 7, 4}
	);

	auto const ext5 = multi::extensions_t<5>{2, 14, 14, 7, 2};

	[[maybe_unused]] auto const ext52 = ext5;

	[[maybe_unused]] multi::array<int, 5> const src_test(ext5);

	multi::array<int, 5> src(
#ifdef _MSC_VER  // problem with 14.3 c++17
		multi::extensions_t<5>
#endif
		{2, 14, 14, 7, 2}
	);

	src[0][1][2][3][1] = 99;

	BOOST_REQUIRE( extensions(tmp.unrotated().partitioned(2).transposed().rotated()) == extensions(src) );
}

BOOST_AUTO_TEST_CASE(rearranged_assignment_resize) {
	multi::array<double, 2> const arrA({4, 5});
	multi::array<double, 2>       arrB({2, 3});

	arrB = arrA;
	BOOST_REQUIRE( arrB.size() == 4 );
}

BOOST_AUTO_TEST_CASE(rvalue_assignments) {
	using complex = std::complex<double>;

	std::vector<double> const vec1(200, 99.0);  // NOLINT(fuchsia-default-arguments-calls)
	std::vector<complex>      vec2(200);  // NOLINT(fuchsia-default-arguments-calls)

	auto linear1 = [&] { return multi::array_cptr<double, 1>(vec1.data(), 200); };
	auto linear2 = [&] { return multi::array_ptr<complex, 1>(vec2.data(), 200); };

	*linear2() = *linear1();
}

BOOST_AUTO_TEST_CASE(assignments) {
	{
		std::vector<int> vec(static_cast<std::size_t>(5 * 7), 99);  // NOLINT(fuchsia-default-arguments-calls)

		constexpr int val = 33;

		multi::array<int, 2> const arr({5, 7}, val);
		multi::array_ref<int, 2>(vec.data(), {5, 7}) = arr();  // arr() is a subarray

		BOOST_REQUIRE( vec[9] == val );
		BOOST_REQUIRE( ! vec.empty() );
		BOOST_REQUIRE( ! is_empty(arr) );
	}
	{
		std::vector<int> vec(5 * 7L, 99);  // NOLINT(fuchsia-default-arguments-calls)
		std::vector<int> wec(5 * 7L, 33);  // NOLINT(fuchsia-default-arguments-calls)

		multi::array_ptr<int, 2> const Bp(wec.data(), {5, 7});
		make_ref(vec.data()) = *Bp;

		auto&& mref = make_ref(vec.data());
		mref        = Bp->sliced(0, 5);

		make_ref(vec.data()) = Bp->sliced(0, 5);

		BOOST_REQUIRE( vec[9] == 33 );
	}
	{
		std::vector<int> vec(5 * 7L, 99);  // NOLINT(fuchsia-default-arguments-calls)
		std::vector<int> wec(5 * 7L, 33);  // NOLINT(fuchsia-default-arguments-calls)

		make_ref(vec.data()) = make_ref(wec.data());

		BOOST_REQUIRE( vec[9] == 33 );
	}
}

template<class T, class Allocator>
auto eye(multi::extensions_t<2> exts, Allocator alloc) {
	multi::array<T, 2, Allocator> ret(exts, 0, alloc);
	ret.diagonal().fill(1);
	return ret;
}

template<class T>
auto eye(multi::extensions_t<2> exts) { return eye<T>(exts, std::allocator<T>{}); }

BOOST_AUTO_TEST_CASE(assigment_temporary) {
	multi::array<int, 2> Id = eye<int>(multi::extensions_t<2>({3, 3}));
	BOOST_REQUIRE( Id == eye<double>({3, 3}) );
	BOOST_REQUIRE( Id[1][1] == 1 );
	BOOST_REQUIRE( Id[1][0] == 0 );
}
