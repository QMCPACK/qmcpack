// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi assignments"
#include<boost/test/unit_test.hpp>

#include<complex>

#include "multi/array.hpp"

namespace multi = boost::multi;

inline static constexpr auto make_ref(double* ptr) -> multi::array_ref<double, 2> {
	return multi::array_ref<double, 2>(ptr, {5, 7});
}

BOOST_AUTO_TEST_CASE(equality_1D) {
	multi::array<double, 1> arr  = {1., 2., 3.};
	multi::array<double, 1> arr2 = {1., 2., 3.};
	BOOST_REQUIRE(      arr == arr2   );
	BOOST_REQUIRE( not (arr != arr2) );

	BOOST_REQUIRE(      arr() == arr2() );
	BOOST_REQUIRE( not (arr() != arr2()) );
}

BOOST_AUTO_TEST_CASE(equality_2D) {
	multi::array<double, 2> arr = {
		{1., 2., 3.},
		{4., 5., 6.}
	};
	multi::array<double, 2> arr2 = {
		{1., 2., 3.},
		{4., 5., 6.}
	};
	BOOST_REQUIRE( arr == arr2 );
	BOOST_REQUIRE( not (arr != arr2) );

	BOOST_REQUIRE( arr() == arr2() );
	BOOST_REQUIRE( not (arr() != arr2()) );

	BOOST_REQUIRE( arr[0] == arr2[0] );
	BOOST_REQUIRE( not (arr[0] != arr2[0]) );
}

BOOST_AUTO_TEST_CASE(multi_copy_move) {
	multi::array<double, 2> arr({3, 3}, 0.);
	multi::array<double, 2> arr2 = arr;
	BOOST_REQUIRE( arr == arr2 );

	auto* arr_data = arr.data_elements();
	multi::array<double, 2> arr3 = std::move(arr);

	BOOST_REQUIRE( arr3.data_elements() == arr_data );

	multi::array<double, 2> arr4(std::move(arr2));
	BOOST_REQUIRE( size(arr4) == 3 );
}

#if 1
BOOST_AUTO_TEST_CASE(range_assignment) {
{
	auto ext = multi::make_extension_t(10L);
	multi::array<double, 1> vec(ext.begin(), ext.end());
	BOOST_REQUIRE( ext.size() == vec.size() );
	BOOST_REQUIRE( vec[1] = 10 );
}
{
	multi::array<double, 1> vec(multi::extensions_t<1>{multi::iextension{10}});
	auto ext = extension(vec);
	vec.assign(ext.begin(), ext.end());
	BOOST_REQUIRE( vec[1] == 1 );
}
}

BOOST_AUTO_TEST_CASE(rearranged_assignment) {
	multi::array<double, 4> tmp({14, 14, 7, 4});
	multi::array<double, 5> src({2, 14, 14, 7, 2}); src[0][1][2][3][1] = 99.;

	BOOST_REQUIRE( extensions(tmp.unrotated().partitioned(2).transposed().rotated()) == extensions(src) );
}

BOOST_AUTO_TEST_CASE(rvalue_assignments) {
	using complex = std::complex<double>;

	std::vector<double > const vec1(200, 99.);
	std::vector<complex>       vec2(200);
	auto linear1 = [&] {return multi::array_cptr<double, 1>(vec1.data(), 200);};
	auto linear2 = [&] {return multi::array_ptr<complex, 1>(vec2.data(), 200);};
	*linear2() = *linear1();
}

#if 0 // self-move-assigment is a standard warning in clang (-Wmove)
BOOST_AUTO_TEST_CASE(self_assigment) {
	multi::array<double, 1> A = {1., 2., 3.};
	A = std::move(A);
	std::cout << A[0] << std::endl;
	BOOST_REQUIRE( A.empty() );

	multi::array<double, 2> B = {{1., 2., 3.}, {2., 3., 4.}};
	B = std::move(B);
	BOOST_REQUIRE( B.empty() );
}
#endif

BOOST_AUTO_TEST_CASE(assignments) {
	{
		std::vector<double> vec( static_cast<std::size_t>(5*7), 99.);
		constexpr double val = 33.;
		multi::array<double, 2> arr({5, 7}, val);
		multi::array_ref<double, 2>(vec.data(), {5, 7}) = arr;
		BOOST_REQUIRE( vec[9] == val );
		BOOST_REQUIRE( not vec.empty() );
		BOOST_REQUIRE( not is_empty(arr) );
	}
	{
		std::vector<double> vec(5*7L, 99.);
		std::vector<double> wec(5*7L, 33.);

		multi::array_ptr<double, 2> Bp{wec.data(), {5, 7}};
		make_ref(vec.data()) = *Bp;
		make_ref(vec.data()) = Bp->sliced(0, 5);

		BOOST_REQUIRE( vec[9] == 33. );
	}
	{
		std::vector<double> vec(5*7L, 99.);
		std::vector<double> wec(5*7L, 33.);

		make_ref(vec.data()) = make_ref(wec.data());

		BOOST_REQUIRE( vec[9] == 33. );
	}
}

template<class T, class Allocator>
auto eye(multi::extensions_t<2> exts, Allocator alloc) {
	multi::array<T, 2, Allocator> ret(exts, 0., alloc);
	ret.diagonal().fill(1.);
	return ret;
}

template<class T>
auto eye(multi::extensions_t<2> exts) {return eye<T>(exts, std::allocator<T>{});}

BOOST_AUTO_TEST_CASE(assigment_temporary) {
	multi::array<double, 2> Id = eye<double>( multi::extensions_t<2>({3, 3}) );
	BOOST_REQUIRE( Id == eye<double>({3, 3}) );
	BOOST_REQUIRE( Id[1][1] == 1 );
	BOOST_REQUIRE( Id[1][0] == 0 );
}

#endif
