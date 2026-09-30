// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>
#include <boost/multi/detail/outer.hpp>

#include <boost/core/lightweight_test.hpp>  // IWYU pragma: keep

#include <iterator>  // IWYU pragma: keep  // for incrementable
#include <tuple>     // IWYU pragma: keep
// IWYU pragma: no_include <type_traits>  // for integral_constant

namespace multi = boost::multi;

auto main() -> int {  // NOLINT(bugprone-exception-escape,readability-function-cognitive-complexity)
	{
		multi::array<int, 1> const arr1d(3);

		auto const x1d = multi::detail::outer_t(arr1d.extent());

		BOOST_TEST( x1d.size() == 3 );

		auto const y1d = multi::detail::outer_t(3);
		BOOST_TEST( y1d.size() == 3 );
	}
	{
		multi::detail::outer_t const x2d(4, 3);
		BOOST_TEST( x2d.size() == 4 );
		auto [x0, x1] = x2d;

		BOOST_TEST( x0.size() == 4 );
		BOOST_TEST( x1.size() == 3 );

		using std::get;
		BOOST_TEST( x0 == get<0>(x2d) );
		BOOST_TEST( x1 == get<1>(x2d) );

		auto it = x2d.begin();
		++it;

#if defined(__cpp_concepts) && (__cpp_concepts >= 201907L) && (!defined(__clang__) || __clang_major__ > 10)
		static_assert(std::incrementable<decltype(it)>);
#endif
		*it;
		x2d.begin()[1];
		// BOOST_TEST( *it == xs[1] );

		auto x2d_it  = x2d.begin();
		auto x2d_it2 = x2d_it + 2;
		auto x2d_it3 = x2d_it2 + 1;
		BOOST_TEST( x2d_it3 == x2d.begin() + 3 );
	}
	{
		// outer_t behaves like a (lazy) array of coordinate tuples
		multi::detail::outer_t const x2d(2, 3);  // 2 x 3 structured cartesian product

		BOOST_TEST( x2d.num_elements() == 6 );

		auto const els = x2d.elements();  // lazy random-access range of (i, j) tuples
		BOOST_TEST( els.size() == 6 );

		using std::get;
		// row-major order: (0,0)(0,1)(0,2)(1,0)(1,1)(1,2)
		BOOST_TEST(( els[0] == std::make_tuple(0, 0) ));
		BOOST_TEST(( els[1] == std::make_tuple(0, 1) ));
		BOOST_TEST(( els[3] == std::make_tuple(1, 0) ));
		BOOST_TEST(( els[5] == std::make_tuple(1, 2) ));

		// random-access iteration
		auto it = els.begin();
		BOOST_TEST(( *it == std::make_tuple(0, 0) ));
		BOOST_TEST(( it[4] == std::make_tuple(1, 1) ));
		BOOST_TEST( els.end() - els.begin() == 6 );
		BOOST_TEST(( *(els.begin() + 5) == std::make_tuple(1, 2) ));

		// operator-(iterator, iterator): distance between two interior positions (probes outer.hpp:151)
		// both operands have nonzero index, so `-`->`+` and operand-swap mutations are caught
		BOOST_TEST( (els.begin() + 5) - (els.begin() + 2) ==  3 );
		BOOST_TEST( (els.begin() + 2) - (els.begin() + 5) == -3 );

		// sizes() reports per-dimension lengths
		BOOST_TEST(( x2d.sizes() == std::make_tuple(2, 3) ));
	}
	{
		multi::array<int, 3> const arr({2, 3, 5});
		auto [is, js, ks] = arr.extents();
		multi::array<int, 3> const brr({ks, js, is}, multi::uninitialized_elements);
	}
	{
		multi::array<int, 1> const arr({5});
		using std::get;
		auto is = get<0>(arr.extents());

		multi::array<int, 1> const brr(multi::extents_t<1>{is});
	}
	{
		multi::array<int, 1> const arr({5});
		auto [is] = arr.extents();

		multi::array<int, 1> const brr(multi::extents_t<1>{is});

		BOOST_TEST( arr.extents() == brr.extents() );
	}
	{
		multi::array<int, 2> const arr({2, 3});

		BOOST_TEST( arr.size() == 2 );
		BOOST_TEST( arr.transposed().size() == 3 );

		auto [is, js] = arr.extents();

		BOOST_TEST( is.size() == 2 );
		BOOST_TEST( js.size() == 3 );

		multi::array<int, 2> const brr(multi::extents_t<2>{js, is});

		BOOST_TEST( brr.size() == 3 );

		BOOST_TEST( brr.size() == arr.transposed().size() );
		BOOST_TEST( brr.extents() == arr.transposed().extents() );
	}
	{
		multi::array<int, 2> const arr({2, 3});

		BOOST_TEST( arr.size() == 2 );
		BOOST_TEST( arr.transposed().size() == 3 );

		auto [is, js] = arr.extents();

		BOOST_TEST( is.size() == 2 );
		BOOST_TEST( js.size() == 3 );

		// braced `{js, is}` would read as iota-rows (extension_t -> array<int,1>); use extents_t<2> for extents
		multi::array<int, 2> const brr(multi::extents_t<2>{js, is});

		BOOST_TEST( brr.size() == 3 );

		BOOST_TEST( brr.size() == arr.transposed().size() );
		BOOST_TEST( brr.extents() == arr.transposed().extents() );
	}
	{
		multi::array<int, 3> const arr({2, 3, 5});
		auto [is, js, ks] = arr.extents();
		multi::array<int, 3> const brr(multi::extents_t<3>{ks, js, is});

		BOOST_TEST( brr.extents() == arr.rotated().transposed().extents() );
	}
	{
		multi::array<int, 3> const arr({2, 3, 5});
		auto [is, js, ks] = arr.extents();
		multi::array<int, 3> const brr({ks, js, is});

		BOOST_TEST( brr.extents() == arr.rotated().transposed().extents() );
	}

	return boost::report_errors();
}
