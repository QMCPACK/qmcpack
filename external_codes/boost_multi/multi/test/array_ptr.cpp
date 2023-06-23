// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo A. Correa 2019-2022

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi array pointer"
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

namespace multi = boost::multi;

// NOLINTNEXTLINE(fuchsia-trailing-return): trailing return helps readability
template<class T> auto fwd_array(T&& array) -> T&& {return std::forward<T>(array);}

BOOST_AUTO_TEST_CASE(multi_array_ptr_equality) {
	multi::array<double, 2> arr = {
		{1., 2., 3.},
		{4., 5., 6.},
		{7., 8., 9.},
		{1., 2., 3.}
	};
	BOOST_REQUIRE(  arr[2] ==  arr[2] );
	BOOST_REQUIRE( &arr[2] == &arr[2] );
	BOOST_REQUIRE( &arr[2] == &fwd_array(arr[2]) );
	BOOST_REQUIRE( &fwd_array(arr[2]) == &arr[2] );

//	auto const& A2 = fwd_array(A[2]);
	auto const& carr2 = arr[2];
	BOOST_REQUIRE( carr2[0] == arr[2][0] );
	BOOST_REQUIRE( carr2.base() == arr[2].base() );
	BOOST_REQUIRE( &carr2 == &std::as_const(arr)[2] );
	BOOST_REQUIRE( &carr2 == &              arr [2] );

	auto const& ac2 = carr2; //fwd_array(A[2]);
	BOOST_REQUIRE( &ac2 == &std::as_const(arr)[2] );
	BOOST_REQUIRE( &ac2 == &              arr [2] );
}

BOOST_AUTO_TEST_CASE(multi_array_ptr) {
	{
		std::array<std::array<double, 5>, 4> arr {{
			{{ 0.,  1.,  2.,  3.,  4.}},
			{{ 5.,  6.,  7.,  8.,  9.}},
			{{10., 11., 12., 13., 14.}},
			{{15., 16., 17., 18., 19.}}
		}};

		multi::array_ptr<double, 2> arrP{&arr};

		BOOST_REQUIRE( arrP->extensions() == multi::extensions(arr) );
		BOOST_REQUIRE( extensions(*arrP) == multi::extensions(arr) );

		using multi::extensions;
		BOOST_REQUIRE( extensions(*arrP) == extensions(arr) );
		BOOST_REQUIRE( &arrP->operator[](1)[1] == &arr[1][1] );

		multi::array_ptr<double, 2> arrP2{&arr};
		BOOST_REQUIRE( arrP == arrP2 ); BOOST_REQUIRE( not (arrP != arrP2) );

		std::array<std::array<double, 5>, 4> arr2{};
		multi::array_ptr<double, 2> arr2P{&arr2};
		BOOST_REQUIRE( arr2P != arrP ); BOOST_REQUIRE( not (arr2P == arrP) );

		arr2P = arrP;
		BOOST_REQUIRE(  arrP ==  arr2P );
		BOOST_REQUIRE( *arrP == *arr2P );
		BOOST_REQUIRE(  arrP->operator==(*arrP) );

		auto&& arrR = *arrP;
		BOOST_REQUIRE( &arrR[1][1] == &arr[1][1] );
		BOOST_REQUIRE( arrR == *arrP );
		BOOST_REQUIRE( std::equal(arrR.begin(), arrR.end(), arrP->begin(), arrP->end()) );
		BOOST_REQUIRE( size(arrR) == arrP->size() );
	}
	{
		std::array<std::array<double, 5>, 4> arr = {{
			std::array<double, 5>{{ 0.,  1.,  2.,  3.,  4.}},
			std::array<double, 5>{{ 5.,  6.,  7.,  8.,  9.}},
			std::array<double, 5>{{10., 11., 12., 13., 14.}},
			std::array<double, 5>{{15., 16., 17., 18., 19.}}
		}};

		std::vector<multi::array_ptr<double, 1>> ptrs;
		ptrs.emplace_back(&arr[0][0], 5);  // NOLINT(readability-container-data-pointer) test access
		ptrs.emplace_back(&arr[2][0], 5);  // NOLINT(readability-container-data-pointer) test access
		ptrs.emplace_back(&arr[3][0], 5);  // NOLINT(readability-container-data-pointer) test access

		BOOST_REQUIRE( &(*ptrs[2])[4] == &arr[3][4]   );
		BOOST_REQUIRE(  (*ptrs[2])[4] == 19         );
		BOOST_REQUIRE(    ptrs[2]->operator[](4) == 19 );
	}
	{
		std::vector<double> v1(100, 3.);
		std::vector<double> const v2(100, 4.);
		multi::array_ptr<double, 2> v1P2D(v1.data(), {10, 10});
		multi::array_cptr<double, 2> v2P2D(v2.data(), {10, 10});

		*v1P2D = *v2P2D;
		v1P2D->operator=(*v2P2D);

		BOOST_REQUIRE( v1[8] == 4. );
	}
}

BOOST_AUTO_TEST_CASE(span_like) {
	std::vector<double> vec = {0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};

	using my_span = multi::array_ref<double, 1>;

	auto aP = & my_span{vec.data() + 2,{5}};                                         // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
	BOOST_REQUIRE( aP->size() == 5 );
	BOOST_REQUIRE( (*aP)[0] == 2. );

	auto const& aCRef = *aP;
	BOOST_REQUIRE(  aCRef.size() == 5 );

	BOOST_REQUIRE( &aCRef[0] == &vec[2] );
	BOOST_REQUIRE(  aCRef[0] == 2.    );

	auto&& aRef = *aP;
	aRef[0] = 99.;
	BOOST_REQUIRE( vec[2] == 99. );
}
