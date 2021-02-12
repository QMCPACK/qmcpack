#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $CXXFLAGS $0 -o $0.$X -lboost_unit_test_framework&&$0.$X&&rm $0.$X;exit
#endif
// © Alfredo Correa 2019-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi zero dimensionality"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include<iostream>

#include "../array.hpp"
//#include "../adaptors/cuda.hpp"

#include<complex>

namespace multi = boost::multi;

template<class T> void what(T&&) = delete;

BOOST_AUTO_TEST_CASE(zero_dimensionality){
	{
		std::vector<double> v1 = {1., 2., 3.};

		multi::array_ref<double, 1> m1(v1.data(), 3);
		BOOST_REQUIRE( size(m1) == 3 );
		BOOST_REQUIRE( &m1[1] == &v1[1] );
		BOOST_REQUIRE( num_elements(m1) == 3 );

		multi::array_ref<double, 0> m0(v1.data());
		BOOST_REQUIRE( &m0 == v1.data() );
		BOOST_REQUIRE( data_elements(m0) == v1.data() );
		BOOST_REQUIRE( num_elements(m0) == 1 );

		m0 = 5.1;
		BOOST_REQUIRE( v1[0] == 5.1 );

		double const& d = std::move(m0);
		BOOST_REQUIRE( d == 5.1 );
	}
	{
		multi::static_array<double, 0> a0 = 45.;
		BOOST_REQUIRE( num_elements(a0) == 1 );
		BOOST_REQUIRE( a0 == 45. );

		a0 = 60.;
		BOOST_REQUIRE( a0 == 60. );
	}
	{
		std::allocator<double> alloc;
		multi::static_array<double, 0> a0(45., alloc);
		BOOST_REQUIRE( num_elements(a0) == 1 );
		BOOST_REQUIRE( a0 == 45. );

		a0 = 60.;
		BOOST_REQUIRE( a0 == 60. );
	}
	{
		multi::array<std::complex<double>, 2> a({1, 2}, std::allocator<std::complex<double>>{});
		assert( size(a) == 1 );
	}
	{
		multi::array<std::complex<double>, 0> a = std::complex<double>{1., 2.};
		assert( num_elements(a) == 1 );
	}
	{
		double d = 2.;
		multi::array_ref<double, 0> ar0(&d, {});
		double dd = ar0;
		BOOST_REQUIRE( dd == d );

		multi::array_ptr<double, 1> ap1(&d, 1);
		BOOST_REQUIRE( ap1->base() == &d );
		BOOST_REQUIRE( (*ap1).base() == &d );
		
		multi::array_ptr<double, 0> ap0 = &d;

		BOOST_REQUIRE( ap0 == &d );
		BOOST_REQUIRE( ap0 != &dd );
		BOOST_REQUIRE( ap0->base() == &d );
		BOOST_REQUIRE( (*ap0).base() == &d );

		multi::array_ptr<double, 0> ap0dd = &dd;
		BOOST_REQUIRE( ap0dd != ap0 );
		BOOST_REQUIRE( *ap0 == *ap0dd );
		double d3 = 3.;
		BOOST_REQUIRE(( *multi::array_ptr<double, 0>{&d3} == 3. ));
		BOOST_REQUIRE(( &multi::array_ref<double, 0>{&d3} == multi::array_ptr<double, 0>{&d3} ));

		#if defined(__cpp_deduction_guides)
		BOOST_REQUIRE(( *multi::array_ptr{&d3} == 3. ));
		BOOST_REQUIRE(( multi::array_ptr{&d3} == multi::array_ptr<double, 0>(&d3, {}) ));
		#endif
	}
#if 0
	{
		multi::cuda::array<double, 0> a0; a0 = 45.;
		multi::cuda::array<double, 0> b0; b0 = 45.;
		BOOST_REQUIRE( a0 == b0 );
	}
	{
		multi::cuda::managed::array<double, 0> a0; a0 = 45.;
		multi::cuda::managed::array<double, 0> b0; b0 = 45.;
		BOOST_REQUIRE( a0 == b0 );
	}
#endif
}

