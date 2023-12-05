// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi reinterpret array"
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

#include<complex>
#include<numeric>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(multi_reinterpret_array_cast_struct_to_dimension) {
	struct vec3 {
		double x, y, z;
	};
	multi::array<vec3, 1> arr(multi::extensions_t<1>{multi::iextension{100}});
	arr[8] = {1., 2., 3.};
	BOOST_REQUIRE( arr[8].y == 2. );

	BOOST_REQUIRE( arr.reinterpret_array_cast<double>(3)[8][1] == arr[8].y );

	multi::array<double, 2> A2D = arr.reinterpret_array_cast<double>(3);
	BOOST_REQUIRE( dimensionality(A2D) == dimensionality(arr) + 1 );
	BOOST_REQUIRE( size(A2D) == size(arr) );
	BOOST_REQUIRE(  A2D[8][1] ==  arr[8].y );
	BOOST_REQUIRE( &A2D[8][1] != &arr[8].y );

	BOOST_REQUIRE( & arr[8].x == & arr.reinterpret_array_cast<double>(3)[8][0] );
	BOOST_REQUIRE( & arr[8].y == & arr.reinterpret_array_cast<double>(3)[8][1] );
	BOOST_REQUIRE( & arr[8].z == & arr.reinterpret_array_cast<double>(3)[8][2] );
}

BOOST_AUTO_TEST_CASE(multi_reinterpret_array_cast_complex_to_real_extra_dimension) {
	using complex = std::complex<double>;
	multi::array<complex, 1> arr(multi::extensions_t<1>{multi::iextension{100}}, complex{1., 2.});
	BOOST_REQUIRE(  size(arr) == 100 );

	BOOST_REQUIRE( real(arr[0]) == 1. );
	BOOST_REQUIRE( imag(arr[0]) == 2. );

	BOOST_REQUIRE(( arr[0] == complex{1., 2.} ));

	multi::array<double, 1> arr2 = arr.reinterpret_array_cast<double>();
	BOOST_REQUIRE( dimensionality(arr2) == dimensionality(arr) );
	BOOST_REQUIRE( arr2[0] == 1 and arr2[1] == 1 );

	multi::array<double, 2> arr3 = arr.reinterpret_array_cast<double>(2);

	BOOST_REQUIRE(( sizes(arr3)==decltype(sizes(arr3)){100, 2} ));
	BOOST_REQUIRE( arr3[5][0] == real(arr[5]) );
	BOOST_REQUIRE( arr3[5][1] == imag(arr[5]) );
}

BOOST_AUTO_TEST_CASE(multi_reinterpret_array_cast_tuple_as_extra_dimension) {
	using vector3 = std::array<double, 3>;
//	using vector3 = std::tuple<double, double, double>; // for tuples reinterpret_array_cast is implementation dependent!!

	vector3 v3d;
	// NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast, cppcoreguidelines-avoid-c-arrays, hicpp-avoid-c-arrays, modernize-avoid-c-arrays): test
	BOOST_REQUIRE( &reinterpret_cast<double(&)[3]>(v3d)[1] == &std::get<1>(v3d) );
	{
		multi::array<vector3, 1> arr(multi::extensions_t<1>{multi::iextension{10}});
		BOOST_REQUIRE( &arr.reinterpret_array_cast<double>(3)[2][1] == &std::get<1>(arr[2]) );
	}
	{
		multi::array<vector3, 2> arr({10, 20});
		BOOST_REQUIRE( &arr.reinterpret_array_cast<double>(3)[5][7][2] == &std::get<2>(arr[5][7]) );
	}
	{
		multi::array<vector3, 2> const arr({4, 5}, vector3{{1., 2., 3.}} );

		BOOST_REQUIRE( dimensionality(arr.reinterpret_array_cast<double>(3)) == 3 );
		BOOST_REQUIRE( arr.reinterpret_array_cast<double>(3).num_elements() == arr.num_elements()*3 );
		BOOST_REQUIRE( arr.reinterpret_array_cast<double>(3).size() == 4 );
		BOOST_REQUIRE( arr.reinterpret_array_cast<double>(3)[0].size() == 5 );
		BOOST_REQUIRE( arr.reinterpret_array_cast<double>(3)[0][0].size() == 3 );
		BOOST_REQUIRE( &arr.reinterpret_array_cast<double>(3)[2][3][0] == &std::get<0>(arr[2][3]) );
		BOOST_REQUIRE( &arr.reinterpret_array_cast<double>(3)[2][3][1] == &std::get<1>(arr[2][3]) );
		BOOST_REQUIRE( &arr.reinterpret_array_cast<double>(3)[2][3][2] == &std::get<2>(arr[2][3]) );

		multi::array<double, 3> const arr2 = arr.reinterpret_array_cast<double>(3);
		BOOST_REQUIRE( arr2[2][3][0] == std::get<0>(arr[2][3]) );
		BOOST_REQUIRE( arr2[2][3][1] == std::get<1>(arr[2][3]) );
		BOOST_REQUIRE( arr2[2][3][2] == std::get<2>(arr[2][3]) );

		auto arr3 = +arr.reinterpret_array_cast<double>(3);
		BOOST_REQUIRE( arr3 == arr2 );
	}
}

template<class T> struct complex_dummy{T real; T imag;};

BOOST_AUTO_TEST_CASE(multi_reinterpret_array_cast) {
{
	std::complex<double> cee{1, 2};
	auto *ptr = reinterpret_cast<complex_dummy<double>*>(&cee); // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
	ptr->real = 11;
	BOOST_REQUIRE(real(cee)==11);
}
{
	multi::array<std::complex<double>, 1> arr(multi::extensions_t<1>{multi::iextension{10}});
	std::iota( begin(arr), end(arr), 1.);
	BOOST_REQUIRE( arr[8] == 9. );
	auto&& arr2 = arr.reinterpret_array_cast<complex_dummy<double>>();
	arr2[8].real = 1000.;
	BOOST_REQUIRE( arr[8] == 1000. );
}
}

BOOST_AUTO_TEST_CASE(multi_reinterpret_array_cast_realcomplex) {
	using complex = std::complex<double>;
{
	complex cee{1, 2};
	auto *conjd_cee = reinterpret_cast<std::array<double, 2>*>(&cee); // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
	(*conjd_cee)[0] = 11;
	BOOST_REQUIRE( conjd_cee );
	BOOST_REQUIRE(real(cee)==11);
}
{
	complex cee{1, 2};
	// NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast, cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays): test purposes
	auto *ceePC = reinterpret_cast<double(*)[2]>(&cee);
	(*ceePC)[0] = 11;
	BOOST_REQUIRE( ceePC );
	BOOST_REQUIRE(real(cee)==11);
}
{
	multi::array<complex, 1> arr(multi::extensions_t<1>{multi::iextension{10}});
	auto&& arr2 = arr.reinterpret_array_cast<double>(2);
	arr2[8][0] = 1000.;
	arr2[8][1] = 2000.;
	BOOST_REQUIRE( arr[8] == std::complex<double>(1000., 2000.) );
}
}

BOOST_AUTO_TEST_CASE(multi_reinterpret_array_cast_pair_to_complex) {
	using complex = std::complex<double>;
	using pair    = std::pair<double, double>;
	multi::array<complex, 2> arr({10, 10}, complex{3., 4.});

	multi::array<complex, 2> const& Aconst = arr;
	auto&& A_block = Aconst({0, 5}, {0, 5});

	auto const& Apair_block = A_block.template reinterpret_array_cast<pair const>(); // const is important // cppcheck 1.90 needs `template` to avoid internal bug
	BOOST_REQUIRE( &Apair_block[1][2] == static_cast<void*>(&arr[1][2]) );

	auto&& Adoubles_block = A_block.reinterpret_array_cast<double const>(2);
	BOOST_REQUIRE( &Adoubles_block[1][2][0] == static_cast<void*>(&arr[1][2]) );
}

BOOST_AUTO_TEST_CASE(multi_reinterpret_array_cast_pointer) {
	multi::array<double, 2> arr({10, 10}, 5.);

	auto&& Aconstcast = arr.reinterpret_array_cast<double, double const*>();
	BOOST_REQUIRE( &arr[0][0] == &Aconstcast[0][0] );
	static_assert( std::is_same<decltype(Aconstcast[1][2]), double const&>{}, "!" );
}
