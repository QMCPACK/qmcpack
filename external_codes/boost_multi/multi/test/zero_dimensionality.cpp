// Copyright 2019-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <complex>

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

BOOST_AUTO_TEST_CASE(zero_dimensionality_part1) {
	{
		std::vector<int> v1 = {10, 20, 30};  // NOLINT(fuchsia-default-arguments-calls)

		multi::array_ref<int, 1> m1(v1.data(), multi::extensions_t<1>{multi::iextension{3}});
		BOOST_REQUIRE( size(m1) == 3 );
		BOOST_REQUIRE( &m1[1] == &v1[1] );
		BOOST_REQUIRE( num_elements(m1) == 3 );

		multi::array_ref<int, 0> m0(v1.data(), {});
//      BOOST_REQUIRE(( &m0 == multi::array_ptr<double, 0>(v1.data(), {}) ));
		BOOST_REQUIRE( data_elements(m0) == v1.data() );
		BOOST_REQUIRE( num_elements(m0) == 1 );

		m0 = 51;
		BOOST_REQUIRE( v1[0] == 51 );

		int const& doub = std::move(m0);
		BOOST_REQUIRE( doub == 51 );
	}
	{
		multi::static_array<double, 0> a0 = multi::static_array<double, 0>{45.0};  // TODO(correaa) this might trigger a compiler crash with g++ 7.5 because of operator&() && overloads
		BOOST_REQUIRE( num_elements(a0) == 1 );
		BOOST_REQUIRE( a0 == 45.0 );

		a0 = multi::static_array<double, 0>{60.0};
		BOOST_REQUIRE( a0 == 60.0 );
	}
	{
		std::allocator<double> const alloc;
		multi::static_array<double, 0> a0(45.0, alloc);
		BOOST_REQUIRE( num_elements(a0) == 1 );
		BOOST_REQUIRE( a0 == 45.0 );

		a0 = multi::static_array<double, 0>{60.0};
		BOOST_REQUIRE( a0 == 60.0 );
	}
}

BOOST_AUTO_TEST_CASE(zero_dimensionality_part2) {
	{
		multi::array<std::complex<double>, 2> const arr(
		#ifdef _MSC_VER  // problem with 14.3 c++17
			multi::extensions_t<2>
		#endif
			{1, 2},
			std::allocator<std::complex<double>>{}
		);
		BOOST_REQUIRE( arr.size() == 1 );
	}
	{
		double doub = 2.0;
		multi::array_ref<double, 0> arr(doub);
		double const& the_doub = static_cast<double&>(arr);
		BOOST_REQUIRE(  the_doub ==  doub );
		BOOST_REQUIRE( &the_doub == &doub );
	}
	{
		double doub = 2.0;
		double dd{multi::array_ref<double, 0>(&doub, {})};

		BOOST_REQUIRE( dd == doub );

		multi::array_ptr<double, 1> const ap1(&doub, multi::extensions_t<1>{{0, 1}});
		BOOST_REQUIRE( ap1->base() == &doub );
		BOOST_REQUIRE( (*ap1).base() == &doub );

		multi::array_ptr<double, 0> const ap0(&doub, {});

		BOOST_REQUIRE(( ap0 == multi::array_ptr<double, 0>(&doub, {}) ));
		BOOST_REQUIRE(( ap0 != multi::array_ptr<double, 0>(&dd, {}) ));
		BOOST_REQUIRE( ap0->base() == &doub );
		BOOST_REQUIRE( (*ap0).base() == &doub );

		multi::array_ptr<double, 0> const ap0dd{&dd};
		BOOST_REQUIRE( ap0dd != ap0 );
		BOOST_REQUIRE( *ap0 == *ap0dd );
		double d3 = 3.14159265358979323846;
		BOOST_REQUIRE(( *multi::array_ptr<double, 0>(&d3, {}) == 3.14159265358979323846 ));
	}
}
