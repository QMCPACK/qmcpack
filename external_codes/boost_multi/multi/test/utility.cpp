#ifdef COMPILATION// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4-*-
$CXX $CXXFLAGS $0 -o $0.$X -lboost_unit_test_framework&&$0.$X&&rm $0.$X;exit
#endif
// © Alfredo A. Correa 2018-2020

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi utility"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"


#include<boost/archive/xml_iarchive.hpp>
#include<boost/archive/xml_oarchive.hpp>

#include<boost/multi_array.hpp>

#include<boost/serialization/binary_object.hpp>

#include "../adaptors/serialization/xml_archive.hpp"

#include<fstream>
#include<numeric> // iota

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(std_array_extensions_3d){
	std::array<std::array<std::array<double, 5>, 4>, 3> arr = {};

	static_assert( std::is_same<typename multi::array_traits<decltype(arr)>::element, double>{}, "!" );

	using multi::dimensionality;
	BOOST_REQUIRE( dimensionality(arr) == 3 );

	using multi::extension;
	BOOST_REQUIRE( extension(arr) == 3 );

	using multi::extensions;
	BOOST_REQUIRE(( extensions(arr) == decltype(extensions(arr)){3, 4, 5} ));

	using multi::data_elements;
	BOOST_REQUIRE( data_elements(arr) == &arr[0][0][0] );

	using multi::num_elements;
	BOOST_REQUIRE( num_elements(arr) == 60 );
	
	multi::array<double, 3> marr({3, 4, 5});
	using multi::layout;
	BOOST_REQUIRE( layout(arr) == layout(marr) );
	
	BOOST_REQUIRE( extensions(arr) == extensions(marr) );
}

BOOST_AUTO_TEST_CASE(std_array_extensions_2d){
	std::array<std::array<double, 4>, 3> arr = {};

	static_assert( std::is_same<typename multi::array_traits<decltype(arr)>::element, double>{}, "!" );

	using multi::dimensionality;
	BOOST_REQUIRE( dimensionality(arr) == 2 );
	
	using multi::extension;
	BOOST_REQUIRE( extension(arr) == 3 );

	using multi::extensions;
	BOOST_REQUIRE(( extensions(arr) == decltype(extensions(arr)){3, 4} ));

	using multi::data_elements;
	BOOST_REQUIRE( data_elements(arr) == &arr[0][0] );

	using multi::num_elements;
	BOOST_REQUIRE( num_elements(arr) == 12 );
	
	multi::array<double, 2> marr({3, 4});
	using multi::layout;
	BOOST_REQUIRE( layout(arr) == layout(marr) );
	
	BOOST_REQUIRE( extensions(arr) == extensions(marr) );
}

BOOST_AUTO_TEST_CASE(std_array_extensions_1d){
	std::array<double, 4> arr = {};

	static_assert( std::is_same<typename multi::array_traits<decltype(arr)>::element, double>{}, "!" );

	using multi::dimensionality;
	BOOST_REQUIRE( dimensionality(arr) == 1 );

	using multi::extension;
	BOOST_REQUIRE( extension(arr) == 4 );

	using multi::extensions;
	BOOST_REQUIRE(( extensions(arr) == decltype(extensions(arr)){4} ));

	using multi::data_elements;
	BOOST_REQUIRE( data_elements(arr) == &arr[0] );

	using multi::num_elements;
	BOOST_REQUIRE( num_elements(arr) == 4 );
}

BOOST_AUTO_TEST_CASE(test_utility_1d){

	std::array<double, 10> carr = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9.};
	multi::array_ref<double, 1> marr(&carr[0], 10);
//	boost::multi_array_ref<double, 1> Marr(&carr[0], boost::extents[10]);
	std::vector<double> varr(10); std::iota(begin(varr), end(varr), 0);
	std::array<double, 10> aarr{}; std::iota(begin(aarr), end(aarr), 0);

	BOOST_REQUIRE( size(marr) == 10 );
	using multi::size;
//	BOOST_REQUIRE( multi::size(varr) == size(marr) );
//	BOOST_REQUIRE( size(Marr) == size(marr) );
	BOOST_REQUIRE( size(carr) == size(marr) );
	BOOST_REQUIRE( static_cast<multi::size_type>(size(aarr)) == size(marr) );

	BOOST_REQUIRE( carr[7] == marr[7] );
//	BOOST_REQUIRE( Marr[7] == marr[7] );
	BOOST_REQUIRE( varr[7] == marr[7] );
	BOOST_REQUIRE( aarr[7] == marr[7] );

	BOOST_REQUIRE( &carr[7] == &marr[7] );
//	BOOST_REQUIRE( &Marr[7] == &marr[7] );

	using multi::num_elements;
	BOOST_REQUIRE( num_elements(carr) == num_elements(marr) );
//	BOOST_REQUIRE( num_elements(Marr) == num_elements(marr) );
	BOOST_REQUIRE( num_elements(varr) == num_elements(marr) );
	BOOST_REQUIRE( num_elements(aarr) == num_elements(aarr) );

//	static_assert( multi::has_data<decltype(Marr)>{}, "!");
//	static_assert( multi::has_num_elements<decltype(Marr)>{}, "!");
//	static_assert( not multi::has_data_elements<decltype(Marr)>{}, "!");

	using multi::data_elements;
	BOOST_REQUIRE( carr.data() == data_elements(marr) );
//	BOOST_REQUIRE( data_elements(Marr) == data_elements(marr) );
//	BOOST_REQUIRE( data_elements(varr) != data_elements(marr) ); // TODO: compat with std::vector


	BOOST_REQUIRE( *begin(varr) == *begin(marr) );
	using std::begin;
//	BOOST_REQUIRE( *begin(Marr) == *begin(marr) );
	BOOST_REQUIRE( *begin(carr) == *begin(marr) );

//	BOOST_TEST( *(rend(varr)-1) == *(rend(marr)-1) );
	using std::rend;
//	BOOST_REQUIRE( *(rend(Marr)-1) == *(rend(marr)-1) );
//	BOOST_REQUIRE( *(rend(carr)-1) == *(rend(marr)-1) );

	BOOST_REQUIRE( *(end(varr)-1) == *(end(marr)-1) );
	using std::end;
//	BOOST_REQUIRE( *(end(Marr)-1) == *(end(marr)-1) );
	BOOST_REQUIRE( *(end(carr)-1) == *(end(marr)-1) );

//	using std::equal;
//	BOOST_REQUIRE( equal(begin(varr), end(varr), begin(marr), end(marr)) );
//	BOOST_REQUIRE( equal(begin(Marr), end(Marr), begin(marr), end(marr)) );
//	BOOST_REQUIRE( equal(begin(carr), end(carr), begin(marr), end(marr)) );

}

BOOST_AUTO_TEST_CASE(test_utility_2d){

	std::array<std::array<double, 10>, 3> carr{
		{
			{ 0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.},
			{10., 11., 12., 13., 14., 15., 16., 17., 18., 19.},
			{20., 21., 22., 23., 24., 25., 26., 27., 28., 29.},
		}
	};
	multi::array_ref<double, 2> marr(&carr[0][0], {3, 10});
//	boost::multi_array_ref<double, 2> Marr(&carr[0][0], boost::extents[3][10]);

	using multi::size;
	BOOST_REQUIRE( size(carr) == size(marr) );
 //	BOOST_REQUIRE( size(Marr) == size(marr) );

	BOOST_REQUIRE( carr[1][7] == marr[1][7] );
//	BOOST_REQUIRE( Marr[1][7] == marr[1][7] );

	BOOST_REQUIRE( &carr[1][7] == &marr[1][7] );
//	BOOST_REQUIRE( &Marr[1][7] == &marr[1][7] );

	using multi::num_elements;
	BOOST_REQUIRE( num_elements(carr) == num_elements(marr) );
//	BOOST_REQUIRE( num_elements(Marr) == num_elements(marr) );

//	static_assert( multi::has_data<decltype(Marr)>{}, "!"); // TODO make array_traits 
//	static_assert( multi::has_num_elements<decltype(Marr)>{}, "!");
//	static_assert( not multi::has_data_elements<decltype(Marr)>{}, "!");

	using multi::data_elements;
	BOOST_REQUIRE( data_elements(carr) == data_elements(marr) );
//	BOOST_REQUIRE( data_elements(Marr) == data_elements(marr) );

}

