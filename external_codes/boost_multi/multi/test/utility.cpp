// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo A. Correa 2018-2021

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi utility"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"

//#include<boost/archive/xml_iarchive.hpp>
//#include<boost/archive/xml_oarchive.hpp>

//#include<boost/serialization/binary_object.hpp>

//#include "../adaptors/serialization/xml_archive.hpp"

#include<fstream>
#include<numeric> // iota

namespace multi = boost::multi;

//BOOST_AUTO_TEST_CASE(c_array_2d){
//	double const A[2][2] = {{1., 2.}, {3., 4.}};
//	BOOST_REQUIRE( std::data(A) == &A[0][0] );
//}

// TODO(correaa) add test for reinterpret_pointer_cast

BOOST_AUTO_TEST_CASE(std_array_extensions_3d) {
	std::array<std::array<std::array<double, 5>, 4>, 3> arr = {};

	static_assert(std::is_same<typename multi::array_traits<decltype(arr)>::element, double>{}, "!");

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

BOOST_AUTO_TEST_CASE(std_array_extensions_2d) {
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

BOOST_AUTO_TEST_CASE(std_array_extensions_1d) {
	std::array<double, 4> arr = {};

	static_assert( std::is_same<typename multi::array_traits<decltype(arr)>::element, double>{}, "!" );

	using multi::dimensionality;
	BOOST_REQUIRE( dimensionality(arr) == 1 );

	using multi::extension;
	BOOST_REQUIRE( extension(arr) == 4 );

	using multi::extensions;
	BOOST_REQUIRE(( extensions(arr) == decltype(extensions(arr)){multi::iextension{4}} ));

	using multi::data_elements;
	BOOST_REQUIRE( data_elements(arr) == &arr[0] );

	using multi::num_elements;
	BOOST_REQUIRE( num_elements(arr) == 4 );
}

BOOST_AUTO_TEST_CASE(test_utility_1d) {
	std::array<double, 10> carr = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9.};
	multi::array_ref<double, 1> marr(&carr[0], {multi::iextension{10}});
//	boost::multi_array_ref<double, 1> Marr(&carr[0], boost::extents[10]);
	std::vector<double> varr(10); std::iota(begin(varr), end(varr), 0);
	std::array<double, 10> aarr{}; std::iota(begin(aarr), end(aarr), 0);

	BOOST_REQUIRE( size(marr) == 10 );

	BOOST_REQUIRE( static_cast<multi::size_t>(carr.size()) == size(marr) );
	BOOST_REQUIRE( static_cast<multi::size_t>(aarr.size()) == size(marr) );

	BOOST_REQUIRE( carr[7] == marr[7] );
	BOOST_REQUIRE( varr[7] == marr[7] );
	BOOST_REQUIRE( aarr[7] == marr[7] );

	BOOST_REQUIRE( &carr[7] == &marr[7] );

	using multi::num_elements;
	BOOST_REQUIRE( num_elements(carr) == num_elements(marr) );
	BOOST_REQUIRE( num_elements(varr) == num_elements(marr) );
	BOOST_REQUIRE( num_elements(aarr) == num_elements(aarr) );

	using multi::data_elements;
	BOOST_REQUIRE( carr.data() == data_elements(marr) );

	BOOST_REQUIRE( *begin(varr) == *begin(marr) );

	using std::begin;
	BOOST_REQUIRE( *begin(carr) == *begin(marr) );

	using std::rend;
	BOOST_REQUIRE( *(end(varr)-1) == *(end(marr)-1) );

	using std::end;
	BOOST_REQUIRE( *(end(carr)-1) == *(end(marr)-1) );
}

BOOST_AUTO_TEST_CASE(test_utility_2d) {
	std::array<std::array<double, 10>, 3> carr{
		{
			{ 0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.},
			{10., 11., 12., 13., 14., 15., 16., 17., 18., 19.},
			{20., 21., 22., 23., 24., 25., 26., 27., 28., 29.},
		}
	};
	multi::array_ref<double, 2> marr(&carr[0][0], {3, 10});

	BOOST_REQUIRE( static_cast<multi::size_t>(carr.size()) == size(marr) );

	BOOST_REQUIRE( carr[1][7] == marr[1][7] );

	BOOST_REQUIRE( &carr[1][7] == &marr[1][7] );

	using multi::num_elements;
	BOOST_REQUIRE( num_elements(carr) == num_elements(marr) );

	using multi::data_elements;
	BOOST_REQUIRE( data_elements(carr) == data_elements(marr) );
}

BOOST_AUTO_TEST_CASE(multi_utility_test) {
	static_assert( std::is_same<std::iterator_traits<double const*>::value_type, double>{}, "!");

	using multi::corigin;
	using multi::dimensionality;
	using multi::extension;
	using multi::extensions;
//  using multi::origin;
	using multi::size;
	using multi::sizes;
	using multi::num_elements;
{
	double A[4] = {1., 2., 3., 4.};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : test legacy types
	BOOST_REQUIRE( dimensionality(A) == 1 );
	BOOST_REQUIRE( extension(A).first() == 0 );
	BOOST_REQUIRE( extension(A).last() == 4 );

	BOOST_REQUIRE( size(A) == 4 );
	BOOST_REQUIRE( std::get<0>(sizes(A)) == size(A) );
	using multi::get_allocator;

	static_assert(std::is_same<decltype(get_allocator(A)), std::allocator<double> >{}, "!");

	using std::addressof;

	using multi::data_elements;
	static_assert( std::is_same<decltype(data_elements(A)), double*>{} , "!");
//	BOOST_REQUIRE( data(A) == addressof(A[0]) );
	BOOST_REQUIRE( data_elements(A) == addressof(A[0]) );
} {
	double A[2][3] = {{1., 2., 3.}, {4., 5., 6.}};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : test legacy types
	BOOST_REQUIRE( dimensionality(A) == 2 );
	BOOST_REQUIRE( extension(A).first() == 0 );
	BOOST_REQUIRE( extension(A).last() == 2 );
//	int a = extensions(A);

//  BOOST_REQUIRE( origin(A) == &A[0][0] );
//	*origin(A) = 99.;
	A[0][0] = 99.;

	BOOST_REQUIRE( A[0][0] == 99. );
	BOOST_REQUIRE( corigin(A) == &A[0][0] );
	BOOST_REQUIRE( size(A) == 2 );
	BOOST_REQUIRE( std::get<0>(sizes(A)) == size(A) );
	BOOST_REQUIRE( num_elements(A) == 6 );

	static_assert( num_elements(A) == 6 , "!" );
}
}

