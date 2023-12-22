// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2023 Alfredo A. Correa

// #define BOOST_TEST_MODULE "C++ Unit Tests for Multi legacy adaptor example"  // NOLINT(cppcoreguidelines-macro-usage) title
#include<boost/test/unit_test.hpp>

#include <multi/array.hpp>

#include<complex>

namespace multi = boost::multi;

template<class DynamicArray>  // e.g. std::vector or multi::array
void resize_copy_1(std::vector<double> const& source, DynamicArray& darr) {
	darr = DynamicArray(source);
}

template<class DynamicArray>  // e.g. std::vector or multi::array
void resize_copy_2(std::vector<double> const& source, DynamicArray& darr) {
	darr = DynamicArray(source.begin(), source.end());  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)
}

template<class DynamicArray>  // e.g. std::vector or multi::array
void resize_copy_3(std::vector<double> const& source, DynamicArray& darr) {
	darr = std::decay_t<decltype(darr)>(source.begin(), source.end());  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)
}

template<class It, class DynamicArray>   // e.g. std::vector or multi::array
void resize_copy_4(It first, It last, DynamicArray& darr) {
	darr = DynamicArray(first, last);  // or std::decay_t<decltype(da)>(source.begin(), source.end())  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)
}

template<class It, class DynamicArray>  // e.g. std::vector or multi::array
void resize_copy_5(It first, It last, DynamicArray& darr) {
	darr.assign(first, last);  // or std::decay_t<decltype(da)>(source.begin(), source.end())
}

// void resize_copy_6   ----> see below test_resize_copy_6

BOOST_AUTO_TEST_CASE(test_resize_copy_1) {
	std::vector<double> const source = {0.0, 1.0, 2.0, 3.0};  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)

	std::vector<double>     dest_v = {99.0, 99.0};  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)
	multi::array<double, 1> dest_a = {88.0, 88.0};

	BOOST_REQUIRE( dest_v.size() == 2 );
	BOOST_REQUIRE( dest_a.size() == 2 );

	resize_copy_1(source, dest_v);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3.0 );

	resize_copy_1(source, dest_a);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3.0 );
}

BOOST_AUTO_TEST_CASE(test_resize_copy_2) {
	std::vector<double> const source = {0.0, 1.0, 2.0, 3.0};  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)

	std::vector<double>     dest_v = {99.0, 99.0};  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)
	multi::array<double, 1> dest_a = {88.0, 88.0};

	BOOST_REQUIRE( dest_v.size() == 2 );
	BOOST_REQUIRE( dest_a.size() == 2 );

	resize_copy_2(source, dest_v);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3. );

	resize_copy_2(source, dest_a);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3. );
}

BOOST_AUTO_TEST_CASE(test_resize_copy_3) {
	std::vector<double> const source = {0.0, 1.0, 2.0, 3.0};  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)

	std::vector<double>     dest_v = {99.0, 99.0};  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)
	multi::array<double, 1> dest_a = {88.0, 88.0};

	BOOST_REQUIRE( dest_v.size() == 2 );
	BOOST_REQUIRE( dest_a.size() == 2 );

	resize_copy_3(source, dest_v);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3. );

	resize_copy_3(source, dest_a);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3. );
}

BOOST_AUTO_TEST_CASE(test_resize_copy_4) {
	std::vector<double> const source = {0.0, 1.0, 2.0, 3.0};  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)

	std::vector<double>     dest_v = {99.0, 99.0};  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)
	multi::array<double, 1> dest_a = {88.0, 88.0};

	BOOST_REQUIRE( dest_v.size() == 2 );
	BOOST_REQUIRE( dest_a.size() == 2 );

	resize_copy_4(source.begin(), source.end(), dest_v);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3.0 );

	resize_copy_4(source.begin(), source.end(), dest_a);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3.0 );
}

BOOST_AUTO_TEST_CASE(test_resize_copy_5) {
	std::vector<double> const source = {0.0, 1.0, 2.0, 3.0};  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)

	std::vector<double>     dest_v = {99.0, 99.0};  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)
	multi::array<double, 1> dest_a = {88.0, 88.0};

	BOOST_REQUIRE( dest_v.size() == 2 );
	BOOST_REQUIRE( dest_a.size() == 2 );

	resize_copy_5(source.begin(), source.end(), dest_v);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3.0 );

	resize_copy_5(source.begin(), source.end(), dest_a);

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3.0 );
}

BOOST_AUTO_TEST_CASE(test_resize_copy_6) {
	std::vector<double> const source = {0.0, 1.0, 2.0, 3.0};  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)

	std::vector<double>     dest_v = {99.0, 99.0};  // testing std::vector vs multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)
	multi::array<double, 1> dest_a = {88.0, 88.0};

	BOOST_REQUIRE( dest_v.size() == 2 );
	BOOST_REQUIRE( dest_a.size() == 2 );

	{  // look same code as below
		dest_v = decltype(dest_v)(source);
	}

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3.0 );

	{  // look same code as above
		dest_a = decltype(dest_a)(source);
	}

	BOOST_REQUIRE( dest_v.size() == 4 );
	BOOST_REQUIRE( dest_v[3] == 3.0 );
}

BOOST_AUTO_TEST_CASE(assign_equality) {
	{
		multi::array<double, 1> const AA = {1.0, 2.0, 3.0};
		std::vector<double> const aa = {1.0, 2.0, 3.0};  // NOLINT(fuchsia-default-arguments-calls)

		BOOST_REQUIRE( std::equal(AA.begin(), AA.end(), aa.begin() ) );
	}
	{
		multi::array<double, 1> const AA = {1.0, 2.0, 3.0};
		std::vector<double> const aa(AA.begin(), AA.end());  // NOLINT(fuchsia-default-arguments-calls)

		BOOST_REQUIRE( std::equal(AA.begin(), AA.end(), aa.begin() ) );
	}
	{
		multi::array<double, 1> const AA = {1.0, 2.0, 3.0};
		auto const aa(AA().operator std::vector<double>());

		BOOST_REQUIRE( std::equal(AA.begin(), AA.end(), aa.begin() ) );
	}
	// {
	//  multi::array<double, 1> const AA = {1.0, 2.0, 3.0};
	//  std::vector<double> const aa(AA);

	//  BOOST_REQUIRE( std::equal(AA.begin(), AA.end(), aa.begin() ) );
	// }
	{
		std::vector<double> const aa = {1.0, 2.0, 3.0};  // NOLINT(fuchsia-default-arguments-calls)
		multi::array<double, 1> const AA(aa.begin(), aa.end());

		BOOST_REQUIRE( std::equal(AA.begin(), AA.end(), aa.begin() ) );
	}
	{
		std::vector<double> const aa = {1.0, 2.0, 3.0};  // NOLINT(fuchsia-default-arguments-calls)
		multi::array<double, 1> const AA(aa);

		BOOST_REQUIRE( std::equal(AA.begin(), AA.end(), aa.begin() ) );
	}
}

BOOST_AUTO_TEST_CASE(construct_from_vector_2D) {
	{
		multi::array<double, 2> const               AA = {{1.0, 2.0}, {3.0, 4.0}};
		BOOST_REQUIRE( AA.num_elements() == 4 );

		std::vector<multi::array<double, 1>> const aa(AA.begin(), AA.end()); // NOLINT(fuchsia-default-arguments-calls)
	}
	{
		multi::array<double, 2> const               AA = {{1.0, 2.0}, {3.0, 4.0}};
		BOOST_REQUIRE( AA.num_elements() == 4 );

		auto const aa(AA().operator std::vector<std::vector<double>>());
	}
	#if not defined(__circle_build__)
	{
		multi::array<double, 2> const               AA = {{1.0, 2.0}, {3.0, 4.0}};
		BOOST_REQUIRE( AA.num_elements() == 4 );

		std::vector<std::vector<double>> const aa(AA);  // not working in circle https://gitlab.com/correaa/boost-multi/-/jobs/5016715421#L757
		// std::vector<std::vector<double>> const aa = AA;  // doesn't compile, needs implicit conversion
	}
	#endif
	{
		multi::array<double, 2> const               AA = {{1.0, 2.0}, {3.0, 4.0}};
		BOOST_REQUIRE( AA.num_elements() == 4 );

	}
	#if not defined(__circle_build__)
	{
		multi::array<double, 2> const              AA = {{1.0, 2.0}, {3.0, 4.0}};
		BOOST_REQUIRE( AA.num_elements() == 4 );

		std::vector<multi::array<double, 1>> const aa(AA);
	}
	#endif
	// {
	//  multi::array<double, 2>                    AA = {{1.0, 2.0}, {3.0, 4.0}};
	//  BOOST_REQUIRE( AA.num_elements() == 4 );

	//  multi::array<multi::array<double, 1>, 1> const aa(AA.begin(), AA.end());  // TODO(correaa)
	// }
	// {
	//  multi::array<double, 2>                    AA = {{1.0, 2.0}, {3.0, 4.0}};
	//  BOOST_REQUIRE( AA.num_elements() == 4 );

	//  multi::array<multi::array<double, 1>, 1> const aa(AA);  // TODO(correaa)
	// }
}
