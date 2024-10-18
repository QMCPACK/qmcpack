// Copyright 2019-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <complex>

// Suppress warnings from boost.test
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wundef"
#pragma clang diagnostic ignored "-Wconversion"
#pragma clang diagnostic ignored "-Wsign-conversion"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wundef"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif

#ifndef BOOST_TEST_MODULE
#define BOOST_TEST_MAIN
#endif

#include <boost/test/unit_test.hpp>

#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

// workaround for libc++ and boost test
// namespace boost::unit_test::ut_detail {
//     auto normalize_test_case_name(const_string name) -> std::string {
//         return ( name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size() ));  // NOLINT(fuchsia-default-arguments-calls)
//     }
// }  // end namespace boost::unit_test::ut_detail

namespace multi = boost::multi;

template<class DynamicArray>  // e.g. std::vector or multi::array
void resize_copy_1(std::vector<int> const& source, DynamicArray& darr) {
	darr = DynamicArray(source);
}

template<class DynamicArray>  // e.g. std::vector or multi::array
void resize_copy_2(std::vector<int> const& source, DynamicArray& darr) {
	darr = DynamicArray(source.begin(), source.end());  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)
}

template<class DynamicArray>  // e.g. std::vector or multi::array
void resize_copy_3(std::vector<int> const& source, DynamicArray& darr) {
	darr = std::decay_t<decltype(darr)>(source.begin(), source.end());  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)
}

template<class It, class DynamicArray>  // e.g. std::vector or multi::array
void resize_copy_4(It first, It last, DynamicArray& darr) {
	darr = DynamicArray(first, last);  // or std::decay_t<decltype(da)>(source.begin(), source.end())  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)
}

template<class It, class DynamicArray>  // e.g. std::vector or multi::array
void resize_copy_5(It first, It last, DynamicArray& darr) {
	darr.assign(first, last);  // or std::decay_t<decltype(da)>(source.begin(), source.end())
}

// void resize_copy_6   ----> see below test_resize_copy_6

BOOST_AUTO_TEST_CASE(test_resize_copy_1) {
	std::vector<int> const source = {0, 1, 2, 3};  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)

	std::vector<int>     dest_v = {99, 99};  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)
	multi::array<int, 1> dest_a = {88, 88};

	BOOST_REQUIRE( dest_v.size() == 2 );
	BOOST_REQUIRE( dest_a.size() == 2 );

	resize_copy_1(source, dest_v);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3 );

	resize_copy_1(source, dest_a);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3 );
}

BOOST_AUTO_TEST_CASE(test_resize_copy_2) {
	std::vector<int> const source = {0, 1, 2, 3};  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)

	std::vector<int>     dest_v = {99, 99};  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)
	multi::array<int, 1> dest_a = {88, 88};

	BOOST_REQUIRE( dest_v.size() == 2 );
	BOOST_REQUIRE( dest_a.size() == 2 );

	resize_copy_2(source, dest_v);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3 );

	resize_copy_2(source, dest_a);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3 );
}

BOOST_AUTO_TEST_CASE(test_resize_copy_3) {
	std::vector<int> const source = {0, 10, 20, 30};  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)

	std::vector<int>     dest_v = {990, 990};  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)
	multi::array<int, 1> dest_a = {880, 880};

	BOOST_REQUIRE( dest_v.size() == 2 );
	BOOST_REQUIRE( dest_a.size() == 2 );

	resize_copy_3(source, dest_v);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 30 );

	resize_copy_3(source, dest_a);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 30 );
}

BOOST_AUTO_TEST_CASE(test_resize_copy_4) {
	std::vector<int> const source = {0, 10, 20, 30};  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)

	std::vector<int>     dest_v = {990, 990};  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)
	multi::array<int, 1> dest_a = {880, 880};

	BOOST_REQUIRE( dest_v.size() == 2 );
	BOOST_REQUIRE( dest_a.size() == 2 );

	resize_copy_4(source.begin(), source.end(), dest_v);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 30 );

	resize_copy_4(source.begin(), source.end(), dest_a);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 30 );
}

BOOST_AUTO_TEST_CASE(test_resize_copy_5) {
	std::vector<int> const source = {0, 10, 20, 30};  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)

	std::vector<int>     dest_v = {990, 990};  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)
	multi::array<int, 1> dest_a = {880, 880};

	BOOST_REQUIRE( dest_v.size() == 2 );
	BOOST_REQUIRE( dest_a.size() == 2 );

	resize_copy_5(source.begin(), source.end(), dest_v);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 30 );

	resize_copy_5(source.begin(), source.end(), dest_a);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 30 );
}

BOOST_AUTO_TEST_CASE(test_resize_copy_6) {
	std::vector<int> const source = {0, 10, 20, 30};  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)

	std::vector<int>     dest_v = {990, 990};  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)
	multi::array<int, 1> dest_a = {880, 880};

	BOOST_REQUIRE( dest_v.size() == 2 );
	BOOST_REQUIRE( dest_a.size() == 2 );

	{  // look same code as below
		dest_v = decltype(dest_v)(source);
	}

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 30 );

	{  // look same code as above
		dest_a = decltype(dest_a)(source);
	}

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 30 );
}

BOOST_AUTO_TEST_CASE(assign_equality) {
	{
		multi::array<int, 1> const AA = {10, 20, 30};
		std::vector<int> const     aa = {10, 20, 30};  // NOLINT(fuchsia-default-arguments-calls)

		BOOST_REQUIRE( std::equal(AA.begin(), AA.end(), aa.begin() ) );
	}
	{
		multi::array<int, 1> const AA = {10, 20, 30};
		std::vector<int> const     aa(AA.begin(), AA.end());  // NOLINT(fuchsia-default-arguments-calls)

		BOOST_REQUIRE( std::equal(AA.begin(), AA.end(), aa.begin() ) );
	}
	{
		multi::array<int, 1> const AA = {10, 20, 30};

		auto const aa(AA().operator std::vector<double>());

		BOOST_REQUIRE( std::equal(AA.begin(), AA.end(), aa.begin() ) );
	}
	// {
	//  multi::array<double, 1> const AA = {1.0, 2.0, 3.0};
	//  std::vector<double> const aa(AA);

	//  BOOST_REQUIRE( std::equal(AA.begin(), AA.end(), aa.begin() ) );
	// }
	{
		std::vector<int> const     aa = {10, 20, 30};  // NOLINT(fuchsia-default-arguments-calls)
		multi::array<int, 1> const AA(aa.begin(), aa.end());

		BOOST_REQUIRE( std::equal(AA.begin(), AA.end(), aa.begin() ) );
	}
	{
		std::vector<int> const     aa = {10, 20, 30};  // NOLINT(fuchsia-default-arguments-calls)
		multi::array<int, 1> const AA(aa);

		BOOST_REQUIRE( std::equal(AA.begin(), AA.end(), aa.begin() ) );
	}
}

BOOST_AUTO_TEST_CASE(construct_from_vector_2D) {
	{
		multi::array<int, 2> const AA = {
			{10, 20},
			{30, 40},
		};
		BOOST_REQUIRE( AA.num_elements() == 4 );

		std::vector<multi::array<double, 1>> const aa(AA.begin(), AA.end());  // NOLINT(fuchsia-default-arguments-calls)
	}
	{
		multi::array<int, 2> const AA = {
			{10, 20},
			{30, 40},
		};
		BOOST_REQUIRE( AA.num_elements() == 4 );

		auto const aa(AA().operator std::vector<std::vector<double>>());
	}
	{
		multi::array<int, 2> const AA = {
			{10, 20},
			{30, 40},
		};
		BOOST_REQUIRE( AA.num_elements() == 4 );

		auto const aa = AA.operator std::vector<std::vector<double>>();
	}
	{
		multi::array<int, 2> const AA = {
			{10, 20},
			{30, 40},
		};
		BOOST_REQUIRE( AA.num_elements() == 4 );

		auto const aa = static_cast<std::vector<std::vector<double>>>(AA);
	}
#if !defined(__circle_build__) || (__circle_build__ > 200)  // crashes circle 187-200 in docker
	{
		multi::array<double, 2> const AA = {
			{1.0, 2.0},
			{3.0, 4.0},
		};

		BOOST_REQUIRE( AA.num_elements() == 4 );

		std::vector<std::vector<double>> const aa(AA);

		BOOST_REQUIRE( aa.size() == 2 );
		// std::vector<std::vector<double>> const aaa = AA;  // doesn't compile, needs implicit conversion
	}
#endif
	{
		multi::array<double, 2> const AA = {
			{1.0, 2.0},
			{3.0, 4.0},
		};
		BOOST_REQUIRE( AA.num_elements() == 4 );
	}
#if !defined(__circle_build__) || (__circle_build__ > 200 )  // crashes circle 187-200 in docker
	{
		multi::array<double, 2> const AA = {
			{1.0, 2.0},
			{3.0, 4.0},
		};
		BOOST_REQUIRE( AA.num_elements() == 4 );

		std::vector<multi::array<double, 1>> const aa(AA);
		BOOST_REQUIRE( aa.size() == 2 );
	}
#endif
}
