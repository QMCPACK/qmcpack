// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo A. Correa 2021-2022

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi index range"
#include<boost/test/unit_test.hpp>

#include "multi/array_ref.hpp"

#include <boost/iterator/transform_iterator.hpp>

// #include <boost/serialization/nvp.hpp>

// #include <boost/archive/xml_iarchive.hpp>
// #include <boost/archive/xml_oarchive.hpp>

#include<numeric>  // for accumulate

namespace multi = boost::multi;

#if 0
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
#endif

BOOST_AUTO_TEST_CASE(multi_range) {
#if defined(__cpp_deduction_guides) and __cpp_deduction_guides and not defined(__NVCC__)
	BOOST_REQUIRE(( multi::range{5, 5}.empty() ));
#else
	BOOST_REQUIRE(( multi::range<std::ptrdiff_t>{5, 5}.empty() ));
#endif
{
	auto drng = multi::range<std::ptrdiff_t>{5, 10};
	std::vector<double> vec(drng.begin(), drng.end());
	BOOST_REQUIRE( vec[1] == 6 );
}
{
	auto drng = multi::range<std::ptrdiff_t>{5, 10};
	auto fun = [](auto idx) {return idx + 1;};
	std::vector<double> vec(
		boost::make_transform_iterator(drng.begin(), fun),
		boost::make_transform_iterator(drng.end()  , fun)
	);
	BOOST_REQUIRE( vec[1] == 7 );
}
}

BOOST_AUTO_TEST_CASE(multi_range_in_constexpr) {
	BOOST_REQUIRE( multi::extension_t<int>{5} == 5 ); // this is not a constexpr in cuda 10
	BOOST_REQUIRE(( multi::extension_t<int>{5, 12}.contains(10) ));

	multi::range<int> irng{5, 12};

	BOOST_REQUIRE( irng.contains(6) );
	BOOST_REQUIRE( not irng.contains(12) );

	BOOST_REQUIRE( * irng.begin()      ==  5 );
	BOOST_REQUIRE( *(irng.begin() + 1) ==  6 );

	BOOST_REQUIRE(   irng.first()      ==  5 );
	BOOST_REQUIRE(   irng.last()       == 12 );

	BOOST_REQUIRE(   irng.front()      ==  5 );
	BOOST_REQUIRE(   irng.back ()      == 11 );

	std::vector<int> vec = {5, 6, 7, 8, 9, 10, 11};

	assert( std::equal( irng.begin(), irng.end(), vec.begin(), vec.end() ) );

	auto sum = 0;
	for(auto elem : irng) {
		sum += elem;
	}
	BOOST_REQUIRE( sum == 5 + 6 + 7 + 8 + 9 + 10 + 11 );

}

BOOST_AUTO_TEST_CASE(multi_range2) {
	multi::index_extension iex(10);

	BOOST_REQUIRE( *begin(iex) == 0 );
	BOOST_REQUIRE( size(iex) == 10 );
	BOOST_REQUIRE( iex[0] == 0 );
	BOOST_REQUIRE( iex[1] == 1 );
	BOOST_REQUIRE( iex[9] == 9 );

	auto const xbeg = begin(iex);
	BOOST_REQUIRE( xbeg[0] == iex[0] );
	BOOST_REQUIRE( xbeg[1] == iex[1] );

	BOOST_REQUIRE( std::accumulate( begin(iex), end(iex), 0) == 0 + 1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9 );

//  #if(__cpp_structured_bindings >= 201606)
	{
		multi::iextensions<3> ies({{0, 3}, {0, 4}, {0, 5}});
		BOOST_REQUIRE( std::get<1>(ies).size() == 4 );
		auto [eyes, jays, kays] = ies;
		BOOST_REQUIRE( eyes.size() == 3 );
		BOOST_REQUIRE( jays.size() == 4 );
		BOOST_REQUIRE( kays.size() == 5 );
	}
//  #endif
}
