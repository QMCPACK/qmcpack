// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo Correa 2019-2021

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi zero dimensionality"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

#include<complex>

namespace multi = boost::multi;

template<class T> void what(T&&) = delete;

BOOST_AUTO_TEST_CASE(zero_dimensionality) {
	{
		std::vector<double> v1 = {1., 2., 3.};

		multi::array_ref<double, 1> m1(v1.data(), multi::extensions_t<1>{multi::iextension{3}});
		BOOST_REQUIRE( size(m1) == 3 );
		BOOST_REQUIRE( &m1[1] == &v1[1] );
		BOOST_REQUIRE( num_elements(m1) == 3 );

		multi::array_ref<double, 0> m0(v1.data(), {});
//		BOOST_REQUIRE(( &m0 == multi::array_ptr<double, 0>(v1.data(), {}) ));
		BOOST_REQUIRE( data_elements(m0) == v1.data() );
		BOOST_REQUIRE( num_elements(m0) == 1 );

		m0 = 5.1;
		BOOST_REQUIRE( v1[0] == 5.1 );

		double const& d = std::move(m0);
		BOOST_REQUIRE( d == 5.1 );
	}
	{
		multi::static_array<double, 0> a0 = multi::static_array<double, 0>{45.}; // TODO(correaa) this might trigger a compiler crash with g++ 7.5 because of operator&() && overloads
		BOOST_REQUIRE( num_elements(a0) == 1 );
		BOOST_REQUIRE( a0 == 45. );

		a0 = multi::static_array<double, 0>{60.};
		BOOST_REQUIRE( a0 == 60. );
	}
	{
		std::allocator<double> alloc;
		multi::static_array<double, 0> a0(45., alloc);
		BOOST_REQUIRE( num_elements(a0) == 1 );
		BOOST_REQUIRE( a0 == 45. );

		a0 = multi::static_array<double, 0>{60.};
		BOOST_REQUIRE( a0 == 60. );
	}
	 {
		multi::array<std::complex<double>, 2> a({1, 2}, std::allocator<std::complex<double>>{});
		BOOST_REQUIRE( size(a) == 1 );
	}
	 {
		double d = 2.;
	//	multi::array_ref<double, 0> ar0(&d, {});
	//	double dd{ar0};
		double dd{multi::array_ref<double, 0>(&d, {})};

		BOOST_REQUIRE( dd == d );

		multi::array_ptr<double, 1> ap1(&d, 1);
		BOOST_REQUIRE( ap1->base() == &d );
		BOOST_REQUIRE( (*ap1).base() == &d );

		multi::array_ptr<double, 0> ap0(&d, {});

		BOOST_REQUIRE(( ap0 == multi::array_ptr<double, 0>(&d, {}) ));
		BOOST_REQUIRE(( ap0 != multi::array_ptr<double, 0>(&dd, {}) ));
		BOOST_REQUIRE( ap0->base() == &d );
		BOOST_REQUIRE( (*ap0).base() == &d );

		multi::array_ptr<double, 0> ap0dd{&dd};
		BOOST_REQUIRE( ap0dd != ap0 );
		BOOST_REQUIRE( *ap0 == *ap0dd );
		double d3 = 3.;
		BOOST_REQUIRE(( *multi::array_ptr<double, 0>(&d3, {}) == 3. ));

		#if defined(__cpp_deduction_guides)
		BOOST_REQUIRE(( *multi::array_ptr {&d3, multi::extensions_t<0>{}} == 3. ));
		BOOST_REQUIRE((  multi::array_ptr {&d3, multi::extensions_t<0>{}} == multi::array_ptr<double, 0>(&d3, {}) ));
		#endif
	}
}

