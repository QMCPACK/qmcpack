// Copyright 2019-2024 Alfredo A. Correa

#include <boost/test/unit_test.hpp>

#include <multi/array.hpp>

#include <array>

namespace multi = boost::multi;

// NOLINTNEXTLINE(fuchsia-trailing-return): trailing return helps readability
template<class T> auto fwd_array(T&& array) -> T&& { return std::forward<T>(array); }

BOOST_AUTO_TEST_CASE(multi_array_ptr_equality) {
	multi::array<double, 2> arr = {
		{1.0, 2.0, 3.0},
		{4.0, 5.0, 6.0},
		{7.0, 8.0, 9.0},
		{1.0, 2.0, 3.0},
	};
	BOOST_REQUIRE(  arr[2] ==  arr[2] );
	BOOST_REQUIRE( &arr[2] == &arr[2] );
	BOOST_REQUIRE( &arr[2] != &(arr[2]({0, 2})) );
	BOOST_REQUIRE( !( &arr[2] == &std::as_const(arr)[2]({0, 2})) );
	BOOST_REQUIRE( &arr[2] == &fwd_array(arr[2]) );
	BOOST_REQUIRE( &fwd_array(arr[2]) == &arr[2] );

	auto arr_ptr = &arr[2];
	BOOST_REQUIRE( arr_ptr == arr_ptr );

	auto& arr_ptr_ref = arr_ptr;
	arr_ptr = arr_ptr_ref;
	arr_ptr = std::move(arr_ptr_ref);

	auto arr_ptr2 = &std::as_const(arr)[2];
	BOOST_REQUIRE( arr_ptr == arr_ptr2 );
	BOOST_REQUIRE( arr_ptr2 == arr_ptr );
	BOOST_REQUIRE( !(arr_ptr != arr_ptr) );

	auto& arr_ptr2_ref = arr_ptr2;
	arr_ptr2 = arr_ptr2_ref;
	arr_ptr2_ref = arr_ptr2;

	auto const& carr2 = arr[2];
	BOOST_REQUIRE( carr2[0] == arr[2][0] );
	BOOST_REQUIRE( carr2.base() == arr[2].base() );
	BOOST_REQUIRE( &carr2 == &std::as_const(arr)[2] );
	BOOST_REQUIRE( &carr2 == &              arr [2] );

	auto const& ac2 = carr2;  // fwd_array(A[2]);
	BOOST_REQUIRE( &ac2 == &std::as_const(arr)[2] );
	BOOST_REQUIRE( &std::as_const(arr)[2] == &ac2 );
	BOOST_REQUIRE( &ac2 == &              arr [2] );
}

BOOST_AUTO_TEST_CASE(multi_array_ptr) {
	{
		// clang-format off
		std::array<std::array<double, 5>, 4> arr{
			{{{0.0, 1.0, 2.0, 3.0, 4.0}},
			 {{5.0, 6.0, 7.0, 8.0, 9.0}},
			 {{10.0, 11.0, 12.0, 13.0, 14.0}},
			 {{15.0, 16.0, 17.0, 18.0, 19.0}}},
		};
		// clang-format on

		multi::array_ptr<double, 2> const arrP{&arr};

		BOOST_REQUIRE( arrP->extensions() == multi::extensions(arr) );
		BOOST_REQUIRE( extensions(*arrP) == multi::extensions(arr) );

		using multi::extensions;
		BOOST_REQUIRE( extensions(*arrP) == extensions(arr) );
		BOOST_REQUIRE( &arrP->operator[](1)[1] == &arr[1][1] );

		multi::array_ptr<double, 2> const arrP2{&arr};
		BOOST_REQUIRE( arrP == arrP2 );
		BOOST_REQUIRE( not (arrP != arrP2) );

		std::array<std::array<double, 5>, 4> arr2{};
		multi::array_ptr<double, 2>          arr2P{&arr2};
		BOOST_REQUIRE( arr2P != arrP );
		BOOST_REQUIRE( not (arr2P == arrP) );

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
		std::array<std::array<double, 5>, 4> arr = {
			{std::array<double, 5>{{0.0, 1.0, 2.0, 3.0, 4.0}},
			 std::array<double, 5>{{5.0, 6.0, 7.0, 8.0, 9.0}},
			 std::array<double, 5>{{10.0, 11.0, 12.0, 13.0, 14.0}},
			 std::array<double, 5>{{15.0, 16.0, 17.0, 18.0, 19.0}}},
		};

		std::vector<multi::array_ptr<double, 1>> ptrs;
		ptrs.emplace_back(&arr[0][0], 5);  // NOLINT(readability-container-data-pointer) test access
		ptrs.emplace_back(arr[2].data(), 5);
		ptrs.emplace_back(&arr[3][0], 5);  // NOLINT(readability-container-data-pointer) test access

		BOOST_REQUIRE( &(*ptrs[2])[4] == &arr[3][4]   );
		BOOST_REQUIRE(  (*ptrs[2])[4] == 19         );
		BOOST_REQUIRE(    ptrs[2]->operator[](4) == 19 );
	}
	{
		std::vector<double>                v1(100, 3.0);  // testing std::vector of multi:array NOLINT(fuchsia-default-arguments-calls)
		std::vector<double> const          v2(100, 4.0);  // testing std::vector of multi:array NOLINT(fuchsia-default-arguments-calls)
		multi::array_ptr<double, 2> const  v1P2D(v1.data(), {10, 10});
		multi::array_cptr<double, 2> const v2P2D(v2.data(), {10, 10});

		*v1P2D = *v2P2D;
		v1P2D->operator=(*v2P2D);

		BOOST_REQUIRE( v1[8] == 4.0 );
	}
}

BOOST_AUTO_TEST_CASE(span_like) {
	std::vector<double> vec = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};  // testing std::vector of multi:array NOLINT(fuchsia-default-arguments-calls)

	using my_span = multi::array_ref<double, 1>;

	auto aP = &my_span{vec.data() + 2, {5}};  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
	BOOST_REQUIRE( aP->size() == 5 );
	BOOST_REQUIRE( (*aP)[0] == 2.0 );

	auto const& aCRef = *aP;
	BOOST_REQUIRE(  aCRef.size() == 5 );

	BOOST_REQUIRE( &aCRef[0] == &vec[2] );
	BOOST_REQUIRE(  aCRef[0] == 2.0     );

	auto&& aRef = *aP;
	aRef[0]     = 99.0;
	BOOST_REQUIRE( vec[2] == 99.0 );
}

BOOST_AUTO_TEST_CASE(multi_array_ptr_assignment) {
	multi::array<double, 2> arr = {
		{1.0, 2.0, 3.0},
		{4.0, 5.0, 6.0},
		{7.0, 8.0, 9.0},
		{1.0, 2.0, 3.0},
	};
	{
		auto rowP = &arr[2];

		rowP = *std::addressof(rowP);

		auto rowP2 = rowP;
		rowP2      = rowP;  // self assigment

		BOOST_REQUIRE( rowP == rowP2 );
		BOOST_REQUIRE( not(rowP != rowP2) );

		auto rowP0 = &arr[0];

		BOOST_REQUIRE( rowP0 != rowP2 );
		BOOST_REQUIRE( not(rowP0 == rowP2) );

		rowP2 = decltype(rowP2){nullptr};
		BOOST_REQUIRE( not rowP2 );

		auto rowP3 = std::exchange(rowP, nullptr);
		BOOST_REQUIRE( rowP3 == &arr[2] );
		BOOST_REQUIRE( rowP == nullptr );
		BOOST_REQUIRE( not rowP );
	}
	{
		auto rowP = &arr();

		rowP = *std::addressof(rowP);

		decltype(rowP) rowP2;
		rowP2 = rowP;

		BOOST_REQUIRE( rowP == rowP2 );

		rowP2 = decltype(rowP2){nullptr};
		BOOST_REQUIRE( not rowP2 );

		auto rowP3 = std::exchange(rowP, nullptr);
		BOOST_REQUIRE( rowP3 == &arr() );
		BOOST_REQUIRE( rowP == nullptr );
		BOOST_REQUIRE( not rowP );
	}
}
