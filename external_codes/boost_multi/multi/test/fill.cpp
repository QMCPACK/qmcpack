// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi fill"
#include<boost/test/unit_test.hpp>

#include "../include/multi/array.hpp"

#include<algorithm>  // for transform
#include<limits>
#include<random>
#include<type_traits> // enable_if_t

// from Howard Hinnart hash
static constexpr auto fnv1a(void const* key, std::size_t len, std::size_t hash) noexcept {  // NOLINT(bugprone-easily-swappable-parameters)
	auto const *first = static_cast<unsigned char const*>(key);
	unsigned char const* const last = first + len;  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic): low level
	for(; first < last; ++first) {  // NOLINT(altera-id-dependent-backward-branch,cppcoreguidelines-pro-bounds-pointer-arithmetic): low level
		hash = (hash ^ *first) * 1099511628211U;  // prime
	}
	return hash;
}

// static constexpr auto fnv1a(void const* key, std::size_t len) noexcept {
// 	return fnv1a(key, len, 14695981039346656037U);
// }

class fnv1a_t {
	std::size_t h = 14695981039346656037U;  // offset

 public:
	using result_type = std::size_t;
	static constexpr auto min() {return std::numeric_limits<result_type>::min();}
	static constexpr auto max() {return std::numeric_limits<result_type>::max();}
	void operator()(void const* key, std::size_t len) noexcept {h = fnv1a(key, len, h);}
	template<class T, std::enable_if_t<std::is_fundamental_v<T>, int> = 0>
	auto operator()(T const& value) noexcept -> decltype(auto) {operator()(&value, sizeof(value)); return *this;}
//  result_type operator()() && noexcept{return h;}
	auto operator()() const& noexcept {return h;}
//  explicit operator result_type() && noexcept {return h;}
	explicit operator result_type() const& noexcept {return h;}
};

BOOST_AUTO_TEST_CASE(fill_1d) {
	namespace multi = boost::multi;
	{
		multi::array<multi::index, 1> d1D(multi::extensions_t<1>{multi::iextension{10}});
		static_assert( std::is_same_v<std::iterator_traits<decltype(begin(d1D))>::value_type, multi::index>, "!");

		using std::copy;
		copy(begin(extension(d1D)), end(extension(d1D)), begin(d1D));
		BOOST_REQUIRE( d1D[0] == 0 );
		BOOST_REQUIRE( d1D[1] == 1 );
		BOOST_REQUIRE( d1D[9] == 9 );

		d1D.assign(extension(d1D));
		BOOST_REQUIRE( d1D[0] == 0 );
		BOOST_REQUIRE( d1D[1] == 1 );
		BOOST_REQUIRE( d1D[9] == 9 );
	}
	{
		multi::array<multi::index, 1> d1D(begin(multi::index_extension(10)), end(multi::index_extension(10)));
		BOOST_REQUIRE( size(d1D) == 10 );
		BOOST_REQUIRE( d1D[0] == 0 );
		BOOST_REQUIRE( d1D[1] == 1 );
		BOOST_REQUIRE( d1D[9] == 9 );
	}
	{
		multi::array<multi::index, 1> d1D(multi::extensions_t<1>{multi::iextension{10}});
		BOOST_REQUIRE( size(d1D) == 10 );

		d1D.assign(begin(extension(d1D)), end(extension(d1D)));
		BOOST_REQUIRE( d1D[0] == 0 );
		BOOST_REQUIRE( d1D[1] == 1 );
		BOOST_REQUIRE( d1D[9] == 9 );
	}
	{
		multi::array<multi::index, 1> d1D(multi::extensions_t<1>{multi::iextension{10}});
		d1D.assign(extension(d1D));
		BOOST_REQUIRE( d1D[0] == 0 );
		BOOST_REQUIRE( d1D[1] == 1 );
		BOOST_REQUIRE( d1D[9] == 9 );
	}
}

BOOST_AUTO_TEST_CASE(fill_member) {
	namespace multi = boost::multi;
	multi::array<double, 1> d1D = {1., 2., 3., 4.};
	d1D.fill(42.);

	multi::array<double, 2> d2D = {
		{150., 16., 17., 18., 19.},
		{  5.,  5.,  5.,  5.,  5.},
		{100., 11., 12., 13., 14.},
		{ 50.,  6.,  7.,  8.,  9.}
	};

	BOOST_REQUIRE(   d2D.elements().size()  == d2D.num_elements()  );
	BOOST_REQUIRE(   d2D.elements().base()  == d2D.base()          );
	BOOST_REQUIRE(   d2D.elements()[3]      == 18.                 );
	BOOST_REQUIRE( &*d2D.elements().begin() == d2D.data_elements() );
	BOOST_REQUIRE( &*d2D.elements().end()   == d2D.data_elements() + d2D.num_elements() );
//	std::fill( d2D.elements().begin(), d2D.elements().end() , 99. );
//	multi::adl_fill_n( d2D.elements().begin(), d2D.elements().size(), 99. );
	d2D.elements().fill(99.);

	BOOST_REQUIRE( d2D[1][1] == 99. );
}

BOOST_AUTO_TEST_CASE(fill) {
	std::random_device randdev;

	namespace multi = boost::multi;

	multi::array<double, 2> d2D = {
		{150., 16., 17., 18., 19.},
		{  5.,  5.,  5.,  5.,  5.},
		{100., 11., 12., 13., 14.},
		{ 50.,  6.,  7.,  8.,  9.}
	};
	using std::all_of;
	BOOST_REQUIRE( all_of(begin(d2D[1]), end(d2D[1]), [](auto const& elem) {return elem == 5.;}) );

	using std::fill;
	fill(d2D[1].begin(), d2D[1].end(), 8.);

	BOOST_REQUIRE( all_of(begin(d2D[1]), end(d2D[1]), [](auto const& elem) {return elem == 8.;}) );

	fill(begin(rotated(d2D)[1]), end(rotated(d2D)[1]), 8.);
	BOOST_REQUIRE( all_of(begin(rotated(d2D)[1]), end(rotated(d2D)[1]), [](auto&& elem) {return elem == 8.;}) );

	fill(begin((d2D.rotated())[1]), end((d2D.rotated())[1]), 8.);
	BOOST_REQUIRE( all_of(begin((d2D.rotated())[1]), end((d2D.rotated())[1]), [](auto&& elem) {return elem == 8.;}) );

	auto rand = [gauss = std::normal_distribution<>{}, gen = std::mt19937{randdev()}]() mutable {return gauss(gen);};
	multi::array<double, 2> r2D({5, 5});
	std::for_each(begin(r2D), end(r2D), [&](auto&& elem) {std::generate(begin(elem), end(elem), rand);});
}

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(fill_1D) {
	multi::array<double, 1> arr = {1., 2., 3.};
	multi::array<double, 2> arr2({10, 3});

	std::fill( begin(arr2), end(arr2), arr );

	BOOST_REQUIRE( arr2[0] == arr );
	BOOST_REQUIRE( arr2[1] == arr );
	// ...
	BOOST_REQUIRE( arr2[9] == arr );
}

#define FWD(a) std::forward<decltype(a)>(a)

template<class BinaryOp, class Column, class Array, class Out>
auto broadcast(BinaryOp op, Column const& col, Array const& in, Out&& out) -> Out&& {  // NOLINT(readability-identifier-length) clang-tidy 14 bug
	std::transform(
		begin(~in), end(~in), begin(~out), begin(~out),
		[acol = (~col)[0], &op](auto const& Acol, auto&& Bcol) {
			std::transform(begin(Acol), end(Acol), begin(acol), begin(Bcol), op);
			return FWD(Bcol);
		}
	);

	return std::forward<Out>(out);
}

BOOST_AUTO_TEST_CASE (julia_broadcast, *boost::unit_test::tolerance(0.00001) ) {
	multi::array<double, 2> col = {
		{0.1},
		{0.2}
	};
	multi::array<double, 2> arr = {
		{1.10813, 1.72068, 1.15387},
		{1.36851, 1.66401, 1.47846}
	};
	{  // "broadcast"
		multi::array<double, 2> arr2(extensions(arr));
		broadcast(std::plus<>{}, col, arr, arr2);

		BOOST_TEST( arr2[0][0] == 1.20813 ); BOOST_TEST( arr2[0][1] == 1.82068 ); BOOST_TEST( arr2[0][2] == 1.25387 );
		BOOST_TEST( arr2[1][0] == 1.56851 ); BOOST_TEST( arr2[1][1] == 1.86401 ); BOOST_TEST( arr2[1][2] == 1.67846 );
	}
	{  // inefficient: replicate the vector before summing elementwise
		multi::array<double, 2> ax3({2, 3});

		std::fill( begin(~ax3), end(~ax3), (~col)[0] );
		BOOST_TEST( ax3[0][0] == 0.1 ); BOOST_TEST( ax3[0][1] == 0.1 ); BOOST_TEST( ax3[0][2] == 0.1 );
		BOOST_TEST( ax3[1][0] == 0.2 ); BOOST_TEST( ax3[1][1] == 0.2 ); BOOST_TEST( ax3[1][2] == 0.2 );

		multi::array<double, 2> Ap(extensions(arr));
		std::transform(begin(arr.elements()), end(arr.elements()), begin(ax3.elements()), begin(Ap.elements()), std::plus<>{});

		BOOST_TEST( Ap[0][0] == 1.20813 ); BOOST_TEST( Ap[0][1] == 1.82068 ); BOOST_TEST( Ap[0][2] == 1.25387 );
		BOOST_TEST( Ap[1][0] == 1.56851 ); BOOST_TEST( Ap[1][1] == 1.86401 ); BOOST_TEST( Ap[1][2] == 1.67846 );
	}
}
