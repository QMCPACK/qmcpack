// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo A. Correa 2021-2022

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi index range"
#include<boost/test/unit_test.hpp>

#include "multi/array_ref.hpp"

#include <boost/iterator/transform_iterator.hpp>

#include <boost/serialization/nvp.hpp>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

#include<numeric>  // for accumulate

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(xml_serialization_index_range) {
	std::stringstream ss;
	multi::range<std::ptrdiff_t> const rg{5, 10};
	{
	    boost::archive::xml_oarchive oa{ss};
		oa<< ::boost::serialization::make_nvp("rg", rg);
	}
	{
		boost::archive::xml_iarchive ia{ss};
		multi::range<std::ptrdiff_t> rg2;
		ia>> ::boost::serialization::make_nvp("rg2", rg2);
		BOOST_REQUIRE( rg == rg2 );
	}
}

BOOST_AUTO_TEST_CASE(multi_range) {
#if defined(__cpp_deduction_guides) and __cpp_deduction_guides and not defined(__NVCC__)
	BOOST_REQUIRE(( multi::range{5, 5}.empty() ));
#else
	BOOST_REQUIRE(( multi::range<std::ptrdiff_t>{5, 5}.empty() ));
#endif
{
	auto r = multi::range<std::ptrdiff_t>{5, 10};
	std::vector<double> v(r.begin(), r.end());
	BOOST_REQUIRE( v[1] == 6 );
}
{
	auto r = multi::range<std::ptrdiff_t>{5, 10};
	auto f = [](auto x) {return x+1;};
	std::vector<double> v(
		boost::make_transform_iterator(r.begin(), f),
		boost::make_transform_iterator(r.end()  , f)
	);
	BOOST_REQUIRE( v[1] == 7 );
}
}

BOOST_AUTO_TEST_CASE(multi_range_in_constexpr) {
	BOOST_REQUIRE( multi::extension_t<int>{5} == 5 ); // this is not a constexpr in cuda 10
	BOOST_REQUIRE(( multi::extension_t<int>{5, 12}.contains(10) ));

	multi::range<int> rr{5, 12};
	BOOST_REQUIRE( rr.contains(6) );
	BOOST_REQUIRE( not rr.contains(12) );
}

BOOST_AUTO_TEST_CASE(multi_range2) {
	multi::index_extension x(10);

	BOOST_REQUIRE( *begin(x) == 0 );
	BOOST_REQUIRE( size(x) == 10 );
	BOOST_REQUIRE( x[0] == 0 );
	BOOST_REQUIRE( x[1] == 1 );
	BOOST_REQUIRE( x[9] == 9 );

	auto b = begin(x);
	BOOST_REQUIRE( b[0] == x[0] );
	BOOST_REQUIRE( b[1] == x[1] );

//	static_assert( ranges::forward_iterator< std::decay_t<decltype(b)> > , "!");

	BOOST_REQUIRE( std::accumulate( begin(x), end(x), 0) == 0 + 1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9 );

//  #if(__cpp_structured_bindings >= 201606)
//   {
//  	multi::iextensions<3> ies({{0, 3}, {0, 4}, {0, 5}});
//  	BOOST_REQUIRE( std::get<1>(ies).size() == 4 );
//  	auto [is, js, ks] = ies;
//  	BOOST_REQUIRE( is.size() == 3 );
//  }
//  #endif
}
