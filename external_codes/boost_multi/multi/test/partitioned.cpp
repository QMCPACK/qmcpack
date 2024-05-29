// Copyright 2018-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <array>

// Suppress warnings from boost.test
#if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wold-style-cast"
#  pragma clang diagnostic ignored "-Wundef"
#  pragma clang diagnostic ignored "-Wconversion"
#  pragma clang diagnostic ignored "-Wsign-conversion"
#  pragma clang diagnostic ignored "-Wfloat-equal"
#elif defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wold-style-cast"
#  pragma GCC diagnostic ignored "-Wundef"
#  pragma GCC diagnostic ignored "-Wconversion"
#  pragma GCC diagnostic ignored "-Wsign-conversion"
#  pragma GCC diagnostic ignored "-Wfloat-equal"
#endif

#ifndef BOOST_TEST_MODULE
#  define BOOST_TEST_MAIN
#endif

#include <boost/test/unit_test.hpp>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(array_partitioned_1d) {
	multi::array<double, 1> A1 = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};

	auto&& A2_ref = A1.partitioned(2);

	static_assert(std::decay<decltype(A2_ref)>::type::rank{} == decltype(A1)::rank{} + 1);
	static_assert(std::decay_t<decltype(A2_ref)>::rank_v == decltype(A1)::rank_v + 1);

	BOOST_REQUIRE( size(A2_ref   ) == 2 );
	BOOST_REQUIRE( size(A2_ref[0]) == 3 );

	BOOST_REQUIRE( &A2_ref[1][0] == &A1[3] );

	BOOST_REQUIRE(( A2_ref == multi::array<double, 2>{ {0, 1, 2}, {3, 4, 5} } ));
}

BOOST_AUTO_TEST_CASE(array_partitioned_2d) {
	multi::array<double, 2> A2 = {
		{ 0.0,  1.0,  2.0,  3.0,  4.0,  5.0},
		{ 6.0,  7.0,  8.0,  9.0, 10.0, 11.0},

		{12.0, 13.0, 14.0, 15.0, 16.0, 17.0},
		{18.0, 19.0, 20.0, 21.0, 22.0, 23.0},
	};
	auto&& A3_ref = A2.partitioned(2);

	static_assert(std::decay_t<decltype(A3_ref)>::rank{} == decltype(A2)::rank{} + 1);
	static_assert(std::decay_t<decltype(A3_ref)>::rank_v == decltype(A2)::rank_v + 1);

	BOOST_REQUIRE( num_elements(A3_ref) == num_elements(A2) );
	BOOST_REQUIRE( size(A3_ref)==2 );
	BOOST_REQUIRE( size(A3_ref[0])==2 );
	BOOST_REQUIRE( size(A3_ref[0][0])==6 );
	BOOST_REQUIRE( &A3_ref[1][1][0] == &A2[3][0] );
}

BOOST_AUTO_TEST_CASE(array_partitioned) {
	using namespace std::string_literals;  // NOLINT(build/namespaces) for ""s

	multi::array<std::string, 2> A2 = {
		{"s0P0"s, "s1P0"s},
		{"s0P1"s, "s1P1"s},
		{"s0P2"s, "s1P2"s},
		{"s0P3"s, "s1P3"s},
		{"s0P4"s, "s1P4"s},
		{"s0P5"s, "s1P5"s},
	};

	BOOST_REQUIRE(  size(A2) == 6 );

	BOOST_REQUIRE( std::get<0>(A2.sizes()) == 6 );
	BOOST_REQUIRE( std::get<1>(A2.sizes()) == 2 );

	BOOST_REQUIRE(( A2.sizes() == decltype(A2.sizes()){6, 2} ));

	BOOST_REQUIRE( std::get<0>(sizes(A2)) == 6 );
	BOOST_REQUIRE( std::get<1>(sizes(A2)) == 2 );

	BOOST_REQUIRE( size(A2.partitioned(3)) == 3 );

	static_assert(decltype(A2.partitioned(3))::rank{} == 3);
	static_assert(decltype(A2.partitioned(3))::rank::value == 3);
	static_assert(decltype(A2.partitioned(3))::rank_v == 3);

	BOOST_REQUIRE(( sizes(A2.partitioned(3)) == decltype(sizes(A2.partitioned(3))){3, 2, 2} ));

	BOOST_REQUIRE( std::get<0>(sizes(A2.partitioned(3))) == 3 );
	BOOST_REQUIRE( std::get<1>(sizes(A2.partitioned(3))) == 2 );
	BOOST_REQUIRE( std::get<2>(sizes(A2.partitioned(3))) == 2 );

	BOOST_REQUIRE( size(A2.partitioned(1)) == 1 );

	static_assert(decltype(A2.partitioned(1))::rank{} == 3);
	static_assert(decltype(A2.partitioned(1))::rank::value == 3);
	static_assert(decltype(A2.partitioned(1))::rank_v == 3);

	BOOST_REQUIRE( &A2.partitioned(1).rotated()[3][1][0] == &A2[3][1] );
}

template<class Ref> class propagate_const;

template<class T> class propagate_const<T&> {
	T& r_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)

 public:
	explicit propagate_const(T& other) : r_{other} {}
	propagate_const(propagate_const const&) = delete;
	propagate_const(propagate_const&&)      = delete;

	auto operator=(propagate_const const&) -> propagate_const&     = default;
	auto operator=(propagate_const&&) noexcept -> propagate_const& = default;

	auto operator=(T const& other) -> propagate_const& {
		r_ = other;
		return *this;
	}

	~propagate_const() noexcept = default;

	explicit operator T const&() const noexcept { return r_; }
	explicit operator T&() noexcept { return r_; }
};

template<class T> class propagate_const<T const&> {
	T const& r_;  // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members)

 public:
	explicit propagate_const(T const& other) : r_{other} {}
	auto operator=(T const& other) -> propagate_const& = delete;
	explicit operator T const&() const noexcept { return r_; }
};

BOOST_AUTO_TEST_CASE(array_encoded_subarray) {
	// arr[walker][encoded_property]  // 7 walkers
	multi::array<int, 2> arr = {
		{990, 990, 1000, 001,  10,  11,  20,  21, 990},
		{990, 990,  100, 101, 110, 111, 120, 121, 990},
		{990, 990,  200, 201, 210, 211, 220, 221, 990},
		{990, 990,  300, 301, 310, 311, 320, 321, 990},
		{990, 990,  400, 401, 410, 411, 420, 421, 990},
		{990, 990,  500, 501, 510, 511, 520, 521, 990},
		{990, 990,  600, 601, 610, 611, 620, 621, 990},
	};

	multi::iextension const encoded_3x2_range = {2, 8};

	auto&& arrRPU = arr.rotated()(encoded_3x2_range).partitioned(3).unrotated();

	static_assert(decltype(+arrRPU)::rank::value == 3);
	static_assert(decltype(+arrRPU)::rank{} == 3);
	static_assert(decltype(+arrRPU)::rank_v == 3);

	BOOST_REQUIRE(( sizes(arrRPU) == decltype(sizes(arrRPU)){7, 3, 2} ));
	BOOST_REQUIRE( arrRPU[4].num_elements() == 3*2L );

	BOOST_REQUIRE( &arrRPU[4][1][0] == &arr[4][4] );
	BOOST_REQUIRE( arrRPU[4][1][0] == 410 );

	BOOST_REQUIRE((
		arrRPU[4] == multi::array<double, 2>{
			{400, 401},
			{410, 411},
			{420, 421},
		}
	));

	arrRPU[4][1][0] = 11110;
	BOOST_REQUIRE( arr[4][4] == 11110 );

	class walker_ref {
		using raw_source_reference = decltype(std::declval<multi::array<int, 2>&>()[0]);
		using internal_array_type  = decltype(std::declval<raw_source_reference>()({2, 8}).partitioned(3));

	 public:  // NOLINT(whitespace/indent) bug in cpplint
		propagate_const<int&> prop1;  // NOLINT(misc-non-private-member-variables-in-classes)
		propagate_const<int&> prop2;  // NOLINT(misc-non-private-member-variables-in-classes)
		internal_array_type slater_array;  // NOLINT(misc-non-private-member-variables-in-classes)
		propagate_const<int&> prop3;  // NOLINT(misc-non-private-member-variables-in-classes)

		explicit walker_ref(raw_source_reference&& row) : prop1{row[0]}, prop2{row[1]}, slater_array{row({2, 8}).partitioned(3)}, prop3{std::move(row)[8]} {}
	};

	auto&& wr = walker_ref(arr[5]);

	wr.prop1 = 88;

	BOOST_REQUIRE( wr.slater_array[2][1] == 521 );

	wr.slater_array[2][1] = 99990;
}

BOOST_AUTO_TEST_CASE(array_partitioned_add_to_last) {
	multi::array<double, 3> arr = {
		{
			{  0.0,  1.0,  2.0,  3.0,  4.0,  5.0},
			{  6.0,  7.0,  8.0,  9.0, 10.0, 11.0},
			{ 12.0, 13.0, 14.0, 15.0, 16.0, 17.0},
			{ 18.0, 19.0, 20.0, 21.0, 22.0, 23.0},
		},
		{
			{  0.0,  1.0,  2.0,  3.0,  4.0,  5.0},
			{  6.0,  7.0,  8.0,  9.0, 10.0, 11.0},
			{ 12.0, 13.0, 14.0, 15.0, 16.0, 17.0},
			{ 18.0, 19.0, 20.0, 21.0, 22.0, 23.0},
		}
	};

	auto strides = std::apply([](auto... strds) { return std::array<std::ptrdiff_t, sizeof...(strds)>{{strds...}}; }, arr.layout().strides());
	// auto strides = std::apply([](auto... strds) { return std::array<std::ptrdiff_t, sizeof...(strds)>{{strds...}}; }, arr.strides());

	BOOST_REQUIRE( std::is_sorted(strides.rbegin(), strides.rend()) && arr.num_elements() == arr.nelems() );  // contiguous c-ordering

#ifndef _MSC_VER  // problem with MSVC 14.3 c++17
	auto&& A4 = arr.reinterpret_array_cast<double>(1);

	BOOST_REQUIRE(( arr.extensions() == decltype(arr.extensions()){2, 4, 6} ));
	BOOST_REQUIRE(( A4.extensions() == decltype(A4.extensions()){2, 4, 6, 1} ));

//  BOOST_REQUIRE( A4.is_flattable() );
//  BOOST_REQUIRE( A4.flatted().is_flattable() );

	BOOST_REQUIRE( &A4[1][2][3][0] == &arr[1][2][3] );
#endif
}

BOOST_AUTO_TEST_CASE(array_partitioned_vs_chunked_1D) {
	multi::array<double, 1> arr = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0};
	BOOST_REQUIRE( size(arr.partitioned(3)) == 3 );
	BOOST_REQUIRE( arr.partitioned(3)[1] == decltype(+arr.partitioned(3)[1])({4.0, 5.0, 6.0, 7.0}) );
	BOOST_REQUIRE( &arr.partitioned(3)[1][2] == &arr[6] );

	BOOST_REQUIRE( size(arr.chunked(3)) == 4 );
	BOOST_REQUIRE( arr.chunked(3)[1] == decltype(+arr.chunked(3)[1])({3.0, 4.0, 5.0}) );
	BOOST_REQUIRE( &arr.chunked(3)[1][2] == &arr[5] );
}

BOOST_AUTO_TEST_CASE(array_partitioned_vs_chunked_2D) {
	multi::array<double, 2> arr({100, 53});
	BOOST_REQUIRE( size(arr.partitioned(20)) == 20 );
	BOOST_REQUIRE( &arr.partitioned(20)[1][2] == &arr[7] );

	BOOST_REQUIRE( size(arr.chunked(5)) == 20 );
	BOOST_REQUIRE( &arr.chunked(5)[1][2] == &arr[7] );
}
