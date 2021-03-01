#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;-*-
$CXX $0 -o $0.$X -lboost_unit_test_framework&&$0.$X&&rm $0.$X;exit
#endif
// © Alfredo A. Correa 2018-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi overload resolution"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

#include<complex>

namespace multi = boost::multi;

std::string what_is(multi::array<             double , 2>&){return "real"   ;}
std::string what_is(multi::array<std::complex<double>, 2>&){return "complex";}

BOOST_AUTO_TEST_CASE(multi_array_range_section){
	multi::array<             double , 2> real_A({10, 20});
	multi::array<std::complex<double>, 2> cplx_A({10, 20});

	std::string real_str    = what_is(real_A);
	std::string complex_str = what_is(cplx_A);

	BOOST_REQUIRE( real_str    == "real"    );
	BOOST_REQUIRE( complex_str == "complex" );
}

