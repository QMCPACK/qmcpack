// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo A. Correa 2018-2021

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi partitioned operation"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(array_partitioned_1d) {
	multi::array<double, 1>	A1 = {0, 1, 2, 3, 4, 5};
	auto&& A2_ref = A1.partitioned(2);

	static_assert( std::decay<decltype(A2_ref)>::type::rank {} == decltype(A1)::rank {} + 1 , "!");
	static_assert( std::decay_t<decltype(A2_ref)>::rank_v == decltype(A1)::rank_v +1 , "!");

	BOOST_REQUIRE( size(A2_ref)==2 );
	BOOST_REQUIRE( size(A2_ref[0])==3 );
	BOOST_REQUIRE( &A2_ref[1][0] == &A1[3] );

	BOOST_REQUIRE(( A2_ref == multi::array<double, 2>{ {0, 1, 2}, {3, 4, 5} } ));
}

BOOST_AUTO_TEST_CASE(array_partitioned_2d) {
	multi::array<double, 2>	A2 = {
		{  0,  1,  2,  3,  4,  5},
		{  6,  7,  8,  9, 10, 11},

		{ 12, 13, 14, 15, 16, 17},
		{ 18, 19, 20, 21, 22, 23},
	};
	auto&& A3_ref = A2.partitioned(2);
	static_assert( std::decay_t<decltype(A3_ref)>::rank_v == decltype(A2)::rank_v + 1 , "!");
	BOOST_REQUIRE( num_elements(A3_ref) == num_elements(A2) );
	BOOST_REQUIRE( size(A3_ref)==2 );
	BOOST_REQUIRE( size(A3_ref[0])==2 );
	BOOST_REQUIRE( size(A3_ref[0][0])==6 );
	BOOST_REQUIRE( &A3_ref[1][1][0] == &A2[3][0] );
}

BOOST_AUTO_TEST_CASE(array_partitioned) {
	multi::array<std::string, 2> A2 = {
		{  "s0P0",  "s1P0"},
		{  "s0P1",  "s1P1"},
		{  "s0P2",  "s1P2"},
		{  "s0P3",  "s1P3"},
		{  "s0P4",  "s1P4"},
		{  "s0P5",  "s1P5"},
	};

	BOOST_REQUIRE(  size(A2) == 6 );
	BOOST_REQUIRE(( sizes(A2) == decltype(sizes(A2)){6, 2} ));

	BOOST_REQUIRE( size(A2.partitioned(3)) == 3 );
	static_assert( decltype(A2.partitioned(3))::rank_v == 3 , "!");

	BOOST_REQUIRE(( sizes(A2.partitioned(3)) == decltype(sizes(A2.partitioned(3))){3, 2, 2} ));

	BOOST_REQUIRE( size(A2.partitioned(1)) == 1 );
	static_assert( decltype(A2.partitioned(1))::rank_v == 3 , "!");
	BOOST_REQUIRE( &A2.partitioned(1).rotated()[3][1][0] == &A2[3][1] );
}

template<class Ref> class propagate_const;

template<class T> class propagate_const<T&>{
	T& r_;

 public:
	explicit propagate_const(T& other) : r_{other} {}
	propagate_const(propagate_const const&) = delete;
	propagate_const(propagate_const&&) = delete;
	// NOLINTNEXTLINE(fuchsia-trailing-return,-warnings-as-errors): reference adaptor
	auto operator=(propagate_const const&) -> propagate_const& = default;
	// NOLINTNEXTLINE(fuchsia-trailing-return,-warnings-as-errors): reference adaptor
	auto operator=(propagate_const&&) noexcept -> propagate_const& = default;
	// NOLINTNEXTLINE(fuchsia-trailing-return,-warnings-as-errors): reference adaptor
	auto operator=(T const& other) -> propagate_const& {r_ = other; return *this;}
	~propagate_const() noexcept = default;
	explicit operator T const&() const noexcept {return r_;}
	explicit operator T      &()       noexcept {return r_;}
};

template<class T> class propagate_const<T const&>{
	T const& r_;

 public:
	explicit propagate_const(T const& other) : r_{other} {}
	// NOLINTNEXTLINE(fuchsia-trailing-return,-warnings-as-errors): reference adaptor
	auto operator=(T const& other) -> propagate_const& = delete;
	explicit operator T const&() const noexcept {return r_;}
};

BOOST_AUTO_TEST_CASE(array_encoded_subarray) {
	multi::array<double, 2> A = { // A[walker][encoded_property] // 7 walkers
		{99., 99., 0.00, 0.01, 0.10, 0.11, 0.20, 0.21, 99.},
		{99., 99., 1.00, 1.01, 1.10, 1.11, 1.20, 1.21, 99.},
		{99., 99., 2.00, 2.01, 2.10, 2.11, 2.20, 2.21, 99.},
		{99., 99., 3.00, 3.01, 3.10, 3.11, 3.20, 3.21, 99.},
		{99., 99., 4.00, 4.01, 4.10, 4.11, 4.20, 4.21, 99.},
		{99., 99., 5.00, 5.01, 5.10, 5.11, 5.20, 5.21, 99.},
		{99., 99., 6.00, 6.01, 6.10, 6.11, 6.20, 6.21, 99.},
	};

	multi::iextension const encoded_3x2_range = {2, 8};
	auto&& B = A.rotated()(encoded_3x2_range).partitioned(3).unrotated();

	static_assert( decltype(+B)::rank_v == 3 , "!");
	BOOST_REQUIRE(( sizes(B) == decltype(sizes(B)){7, 3, 2} ));
	BOOST_REQUIRE( B[4].num_elements() == 3*2L );

	BOOST_REQUIRE( &B[4][1][0] == &A[4][4] );
	BOOST_REQUIRE( B[4][1][0] == 4.10 );

	BOOST_REQUIRE((
		B[4] == multi::array<double, 2>{
			{4.00, 4.01},
			{4.10, 4.11},
			{4.20, 4.21},
		}
	));

	B[4][1][0] = 1111.;
	BOOST_REQUIRE( A[4][4] == 1111. );

	class walker_ref{
		using raw_source_reference = decltype(std::declval<multi::array<double, 2>&>()[0]);
		using internal_array_type = decltype(std::declval<raw_source_reference>()({2, 8}).partitioned(3));
	public:
		propagate_const<double&> prop1;   // NOLINT(misc-non-private-member-variables-in-classes)
		propagate_const<double&> prop2;   // NOLINT(misc-non-private-member-variables-in-classes)
		internal_array_type slater_array; // NOLINT(misc-non-private-member-variables-in-classes)
		propagate_const<double&> prop3;   // NOLINT(misc-non-private-member-variables-in-classes)
		explicit walker_ref(raw_source_reference&& row) : prop1{row[0]}, prop2{row[1]}, slater_array{row({2, 8}).partitioned(3)}, prop3{row[8]}{}
	};

	auto&& wr = walker_ref(A[5]);
	wr.prop1 = 88;
	BOOST_REQUIRE( wr.slater_array[2][1] == 5.21 );

	wr.slater_array[2][1] = 9999.;
}

BOOST_AUTO_TEST_CASE(array_partitioned_add_to_last) {
	multi::array<double, 3>	A3 = {
		{
			{  0.,  1.,  2.,  3.,  4.,  5.},
			{  6.,  7.,  8.,  9., 10., 11.},
			{ 12., 13., 14., 15., 16., 17.},
			{ 18., 19., 20., 21., 22., 23.},
		},
		{
			{  0.,  1.,  2.,  3.,  4.,  5.},
			{  6.,  7.,  8.,  9., 10., 11.},
			{ 12., 13., 14., 15., 16., 17.},
			{ 18., 19., 20., 21., 22., 23.},
		}
	};

#if(__cplusplus >= 201703L)
	auto strides = std::apply([](auto... e){return std::array<std::ptrdiff_t, sizeof...(e)>{e...};}, A3.strides());

	BOOST_REQUIRE( std::is_sorted(strides.rbegin(), strides.rend()) and A3.num_elements() == A3.nelems() ); // contiguous c-ordering
#endif

	auto&& A4 = A3.reinterpret_array_cast<double>(1);

	BOOST_REQUIRE(( A3.extensions() == decltype(A3.extensions()){2, 4, 6} ));
	BOOST_REQUIRE(( A4.extensions() == decltype(A4.extensions()){2, 4, 6, 1} ));

	BOOST_REQUIRE( A4.is_flattable() );
	BOOST_REQUIRE( A4.flatted().is_flattable() );

	BOOST_REQUIRE( &A4[1][2][3][0] == &A3[1][2][3] );
}

BOOST_AUTO_TEST_CASE(array_partitioned_vs_chunked_1D) {
	multi::array<double, 1> A = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11.};
	BOOST_REQUIRE( size(A.partitioned(3)) == 3 );
	BOOST_REQUIRE( A.partitioned(3)[1] == decltype(+A.partitioned(3)[1])({4., 5., 6., 7.}) );
	BOOST_REQUIRE( &A.partitioned(3)[1][2] == &A[6] );

	BOOST_REQUIRE( size(A.chunked(3)) == 4 );
	BOOST_REQUIRE( A.chunked(3)[1] == decltype(+A.chunked(3)[1])({3., 4., 5.}) );
	BOOST_REQUIRE( &A.chunked(3)[1][2] == &A[5] );
}

BOOST_AUTO_TEST_CASE(array_partitioned_vs_chunked_2D) {
	multi::array<double, 2> A({100, 53});
	BOOST_REQUIRE( size(A.partitioned(20)) == 20 );
	BOOST_REQUIRE( &A.partitioned(20)[1][2] == &A[7] );

	BOOST_REQUIRE( size(A.chunked(5)) == 20 );
	BOOST_REQUIRE( &A.chunked(5)[1][2] == &A[7] );
}

