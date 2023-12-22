// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2023 Alfredo A. Correa

// #define BOOST_TEST_MODULE "C++ Unit Tests for Multi allocators"  // NOLINT(cppcoreguidelines-macro-usage) title
#include <boost/test/unit_test.hpp>

#include <multi/array.hpp>

#if defined(__cpp_lib_memory_resource) and (__cpp_lib_memory_resource >= 201603L)
#include <memory_resource>  // for polymorphic memory resource, monotonic buffer
#endif

#include <vector>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(std_vector_of_arrays) {
	std::vector<multi::array<double, 2>> va;
	std::transform(
		begin(multi::iextension(3)), end(multi::iextension(3)),
		std::back_inserter(va),
		[](auto idx){return multi::array<double, 2>({idx, idx}, static_cast<double>(idx));}
	);

	BOOST_REQUIRE( size(va[0]) == 0 );
	BOOST_REQUIRE( size(va[1]) == 1 );
	BOOST_REQUIRE( size(va[2]) == 2 );
	BOOST_REQUIRE( va[1] [0][0] == 1 );
	BOOST_REQUIRE( va[2] [0][0] == 2 );

	std::vector<multi::array<double, 2>> const wa = {  // testing std::vector of multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)
		multi::array<double, 2>({0, 0}, 0.0),
		multi::array<double, 2>({1, 1}, 1.0),
		multi::array<double, 2>({2, 2}, 2.0),
	};
	BOOST_REQUIRE( size(va) == size(wa) );
	BOOST_REQUIRE( va == wa );

	std::vector<multi::array<double, 2>> ua(3, std::allocator<multi::array<double, 2>>{});
	auto iex = multi::iextension(static_cast<multi::size_type>(ua.size()));
	std::transform(
		begin(iex), end(iex),
		begin(ua),
		[](auto idx) {return multi::array<double, 2>({idx, idx}, static_cast<double>(idx));}
	);
	BOOST_REQUIRE( ua == va );
}

BOOST_AUTO_TEST_CASE(array1d_of_arrays2d) {
	multi::array<multi::array<double, 2>, 1> arr(multi::extensions_t<1>(multi::iextension{10}), multi::array<double, 2>{});
	BOOST_REQUIRE( size(arr) == 10 );

	std::transform(
		begin(extension(arr)), end(extension(arr)), begin(arr),
		[](auto idx) {return multi::array<double, 2>({idx, idx}, static_cast<double>(idx));}
	);

	BOOST_REQUIRE( size(arr[0]) == 0 );
	BOOST_REQUIRE( size(arr[1]) == 1 );
	BOOST_REQUIRE( size(arr[8]) == 8 );
	BOOST_REQUIRE( arr[8][4][4] == 8.0 );
}

BOOST_AUTO_TEST_CASE(array_3d_of_array_2d)  {
	multi::array<multi::array<double, 3>, 2> AA({10, 20}, multi::array<double, 3>{});
	std::transform(extension(AA).begin(), extension(AA).end(), AA.begin(), AA.begin(), [](auto idx, auto&& row) -> decltype(row) {
		std::transform(extension(row).begin(), extension(row).end(), row.begin(), [idx](auto jdx) {
			return multi::array<double, 3>({idx + jdx, idx + jdx, idx + jdx}, 99.0);
		});
		return std::forward<decltype(row)>(row);
	});

	BOOST_REQUIRE( size(AA[9][19]) == 9 + 19 );
	BOOST_REQUIRE( AA[9][19][1][1][1] == 99.0 );
}

BOOST_AUTO_TEST_CASE(array_3d_of_array_2d_no_init)  {
	multi::array<multi::array<double, 3>, 2> AA({10, 20});
	std::transform(extension(AA).begin(), extension(AA).end(), AA.begin(), AA.begin(), [](auto idx, auto&& row) -> decltype(row) {
		std::transform(extension(row).begin(), extension(row).end(), row.begin(), [idx](auto jdx) {
			return multi::array<double, 3>({idx + jdx, idx + jdx, idx + jdx}, 99.0);
		});
		return std::forward<decltype(row)>(row);
	});

	BOOST_REQUIRE( size(AA[9][19]) == 9 + 19 );
	BOOST_REQUIRE( AA[9][19][1][1][1] == 99. );
}


BOOST_AUTO_TEST_CASE(const_elements) {
	auto ptr = std::make_unique<double const>(2.0);
//  *ptr = 3.0;  // ok, can't assign
	BOOST_REQUIRE( *ptr == 2.0 );

//  multi::array<double const, 2, std::allocator<double>> arr({10, 10}, 99.0);
//
//  BOOST_REQUIRE( arr[1][2] == 99.0 );
}

#if defined(__cpp_lib_memory_resource) and (__cpp_lib_memory_resource >= 201603L)
BOOST_AUTO_TEST_CASE(pmr) {
	std::array<char, 13> buffer = {{'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X'}};
	std::pmr::monotonic_buffer_resource pool{std::data(buffer), std::size(buffer)};

	multi::array<char, 2, std::pmr::polymorphic_allocator<char>> Aarr({2, 2}, 'a', &pool);
	multi::array<char, 2, std::pmr::polymorphic_allocator<char>> Barr({3, 2}, 'b', &pool);

	BOOST_REQUIRE(( buffer == std::array<char, 13>{{'a', 'a', 'a', 'a', 'b', 'b', 'b', 'b', 'b', 'b', 'X', 'X', 'X'}} ));

	BOOST_REQUIRE(Aarr[0][0] == 'a');
	BOOST_REQUIRE(Barr[0][0] == 'b');
}
#endif

#if(MULTI_PROVIDES_PMR_ARRAY)
BOOST_AUTO_TEST_CASE(pmr2) {
	std::array<char, 13> buffer = {{'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X'}};
	std::pmr::monotonic_buffer_resource pool{std::data(buffer), std::size(buffer)};

	multi::pmr::array<char, 2> Aarr({2, 2}, 'a', &pool);
	multi::pmr::array<char, 2> Barr({3, 2}, 'b', &pool);

	BOOST_REQUIRE(( buffer == std::array<char, 13>{{'a', 'a', 'a', 'a', 'b', 'b', 'b', 'b', 'b', 'b', 'X', 'X', 'X'}} ));

	BOOST_REQUIRE(Aarr[0][0] == 'a');
	BOOST_REQUIRE(Barr[0][0] == 'b');
}

BOOST_AUTO_TEST_CASE(pmr_double_uninitialized) {
	std::array<double, 12> buffer = {{4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.00, 11.0,  999.9, 999.9, 999.9, 999.9}};
	std::pmr::monotonic_buffer_resource pool{static_cast<void*>(std::data(buffer)), 12*sizeof(double)};

	multi::pmr::array<double, 2> Aarr({2, 2}, &pool);

	BOOST_TEST( buffer[0] == 4.0 );
	BOOST_TEST( buffer[1] == 5.0 );

	BOOST_REQUIRE(Aarr[0][0] == 4.0);
}
#endif
