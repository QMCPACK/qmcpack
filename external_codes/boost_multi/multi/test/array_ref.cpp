// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi array reference"
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

#include<numeric>  // for std::iota

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(array_ref_from_carray) {
	// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays): test
	double arr[4][5] {
		{ 0.,  1.,  2.,  3.,  4.},
		{ 5.,  6.,  7.,  8.,  9.},
		{10., 11., 12., 13., 14.},
		{15., 16., 17., 18., 19.}
	};

	multi::array_ptr<double, 2> map{&arr};
	BOOST_REQUIRE( &map->operator[](1)[1] == &arr[1][1] );
	BOOST_REQUIRE( (*&arr)[1][1] == 6. );

	multi::array_ref<double, 2>&& mar = *map;

	BOOST_REQUIRE( &mar[1][1] == &arr[1][1] );

	mar[1][1] = 9.;
	BOOST_REQUIRE( &mar[1][1] == &arr[1][1] );

	auto const& a_const = arr;
//  double const(&a_const)[4][5] = a;
	BOOST_REQUIRE( &a_const[1][1] == &arr[1][1] );

	static_assert( decltype(mar(2, {1, 3}))::rank_v == 1 , "!");

	BOOST_REQUIRE( size(mar(2, {1, 3})) == 2 );
	BOOST_REQUIRE( &mar(2, {1, 3})[1] == &arr[2][2] );
}

BOOST_AUTO_TEST_CASE(array_ref_1D_reindexed) {
	std::array<std::string, 5> stdarr{ {"a", "b", "c", "d", "e"} };

	multi::array_ref<std::string, 1> mar = *multi::array_ptr<std::string, 1>(&stdarr);

	BOOST_REQUIRE( &mar[1] == &stdarr[1] );
	BOOST_REQUIRE( sizes(mar.reindexed(1)) == sizes(mar) );

	auto diff = &(mar.reindexed(1)[1]) - &mar[0];
	BOOST_REQUIRE( diff == 0 );

	BOOST_REQUIRE( &mar.blocked(2, 4)[2] == &mar[2] );
	for(auto idx : extension(mar.stenciled({2, 4}))) {
		BOOST_REQUIRE( &mar.stenciled({2, 4})[idx] == &mar[idx] );
	}

	multi::array<std::string, 1> arr({{2, 7}}, std::string{"xx"});
	BOOST_REQUIRE( size(arr) == 5 );
	BOOST_REQUIRE( extension(arr) == multi::iextension(2, 7) );
	arr[2] = "a";
	arr[3] = "b";
	arr[4] = "c";
	arr[5] = "d";
	arr[6] = "e";
	BOOST_REQUIRE( std::equal(arr.begin(), arr.end(), mar.begin(), mar.end()) );

	auto arrB = multi::array<std::string, 1>({"a", "b", "c", "d", "e"}).reindex(2);
	BOOST_REQUIRE( size(arrB) == 5 );
	BOOST_REQUIRE( arrB[2] == "a" );
	BOOST_REQUIRE( arrB[6] == "e" );
}

BOOST_AUTO_TEST_CASE(array_ref_of_nested_std_array_reindexed) {
	std::array<std::array<double, 5>, 4> arr = {{
		{{ 0.,  1.,  2.,  3.,  4.}},
		{{ 5.,  6.,  7.,  8.,  9.}},
		{{10., 11., 12., 13., 14.}},
		{{15., 16., 17., 18., 19.}}
	}};

	// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays): test type
	multi::array_ref<double, 2> mar = *multi::array_ptr<double, 2>(&arr);
	BOOST_REQUIRE( &mar[1][1] == &arr[1][1] );
}

BOOST_AUTO_TEST_CASE(array_ref_reindexed) {
	double (&&arr)[4][5] = {  // NOLINT(hicpp-avoid-c-arrays, modernize-avoid-c-arrays, cppcoreguidelines-avoid-c-arrays): test
		{ 0,  1,  2,  3,  4},
		{ 5,  6,  7,  8,  9},
		{10, 11, 12, 13, 14},
		{15, 16, 17, 18, 19}
	};

	// NOLINTNEXTLINE(hicpp-avoid-c-arrays, modernize-avoid-c-arrays, cppcoreguidelines-avoid-c-arrays): special type
	multi::array_ref<double, 2> mar = *multi::array_ptr<double, 2>(&arr);

	BOOST_REQUIRE( &mar[1][1] == &arr[1][1] );
	BOOST_REQUIRE( size(mar   .reindexed(1)) == size(mar) );
	BOOST_REQUIRE( size(mar[0].reindexed(1)) == size(mar[0]) );

	BOOST_REQUIRE( sizes(mar.reindexed(1)) == sizes(mar) );

	BOOST_REQUIRE( &mar.reindexed(1)[1][0] == &mar[0][0] );

	BOOST_REQUIRE( sizes(mar[0].reindexed(1)) == sizes(mar[0]) );
	BOOST_REQUIRE( mar[0].reindexed(1).extension().start () == mar[0].extension().start () + 1 );
	BOOST_REQUIRE( mar[0].reindexed(1).extension().finish() == mar[0].extension().finish() + 1 );

	auto diff = &mar[0].reindexed(1)[1] - &mar[0][0];
	BOOST_REQUIRE( diff == 0 );

//	BOOST_REQUIRE( &(((mar<<1).reindexed(2)>>1).reindexed(1))[1][2] == &mar[0][0] );
	BOOST_REQUIRE( &mar.reindexed(1, 2)[1][2] == &mar[0][0] );

	BOOST_REQUIRE( &mar.reindexed(1)({1, 5})[1][0] == &mar[0][0] );

	BOOST_REQUIRE(( sizes(mar.stenciled({2, 4})) == decltype(sizes(mar.stenciled({2, 4}))){2, 5} ));
	BOOST_REQUIRE( &mar.stenciled({2, 4})[2][0] == &mar[2][0] );
	BOOST_REQUIRE( &mar.stenciled({2, 4}, {1, 3})[2][1] == &mar[2][1] );

//	BOOST_REQUIRE( &mar[0][0] == mar.origin() ); // origin changed meaning in on 2020/Dec/16
//	BOOST_REQUIRE( mar.base() == mar.origin() );

//	BOOST_REQUIRE( mar.stenciled({2, 4}).origin() == mar.origin() );  // origin changed meaning in on 2020/Dec/16
	BOOST_REQUIRE( mar.stenciled({2, 4}).base()   != mar.base()   );

	BOOST_REQUIRE( &mar.stenciled({2, 4})[2][0] == mar.stenciled({2, 4}).base() );

	{
		multi::array<std::string, 2> arrB = {
			{"a", "b", "c", "d", "e"},
			{"f", "g", "h", "f", "g"},
			{"h", "i", "j", "k", "l"}
		};
		arrB.reindex(2);
		BOOST_REQUIRE( size(arrB) == 3 );
		BOOST_REQUIRE( arrB[2][0] == "a" );
	}
	{
		multi::array<std::string, 2> arrB = {
			{"a", "b", "c", "d", "e"},
			{"f", "g", "h", "f", "g"},
			{"h", "i", "j", "k", "l"}
		};
		arrB.reindex(2, 1);
		BOOST_REQUIRE( size(arrB) == 3 );
		BOOST_REQUIRE( arrB[2][1] == "a" );
	}
	{
		multi::array<std::string, 2> arrB = (multi::array<std::string, 2>
			{{"a", "b", "c", "d", "e"},
			 {"f", "g", "h", "f", "g"},
			 {"h", "i", "j", "k", "l"}})//.reindex(2, 1);
		;
		BOOST_REQUIRE( arrB.reindex(2).extension() == multi::iextension(2, 5) );
		auto exts = arrB.reindexed(2).extensions();

		multi::array<std::string, 2> arrC(exts);
		BOOST_REQUIRE( size(arrC) == 3 );
		BOOST_REQUIRE( size(arrC) == size(arrB) );

		BOOST_REQUIRE( arrC.extension().start()  == 2 );
		BOOST_REQUIRE( arrC.extension().finish() == 5 );
	}
}

BOOST_AUTO_TEST_CASE(array_ref_with_stencil) {
	std::array<std::array<double, 5>, 4> arr = {{
		{{ 0.,  1.,  2.,  3.,  4.}},
		{{ 5.,  6.,  7.,  8.,  9.}},
		{{10., 11., 12., 13., 14.}},
		{{15., 16., 17., 18., 19.}}
	}};
	// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test type
	auto const& mar = *multi::array_ptr<double, 2>(&arr);

	// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test type
	multi::array<double, 2> ss = {
		{ 0., +1.,  0.},
		{+1., -4., +1.},
		{ 0., +1.,  0.}
	};
	auto const& stencil = ss.reindexed(-1, -1);

	BOOST_REQUIRE( &stencil[-1][-1] == stencil.base() );

	multi::array<double, 2> gy(extensions(mar), 0.);

	 {
		auto xs = extensions(mar);
		for(auto eye = std::get<0>(xs).start() + 1; eye != std::get<0>(xs).finish()-1; ++eye) {
			for(auto jay = std::get<1>(xs).start() + 1; jay != std::get<1>(xs).finish() - 1; ++jay) {
				auto xt = extensions(stencil);
				for(auto kay : std::get<0>(xt)) {
					for(auto ell : std::get<1>(xt)) {
						gy[eye][jay] += stencil[kay][ell]*mar[eye + kay][jay + ell];
					}
				}
			}
		}
	}
}

BOOST_AUTO_TEST_CASE(array_ref_1D_from_vector) {
	std::vector<double> vec = {1, 2, 3};
	multi::array_ref<double, 1> aref({{1, 3}}, vec.data());
	BOOST_REQUIRE( aref.extension() == multi::iextension(1, 3) );
	BOOST_REQUIRE( &aref[1] == vec.data() );
}

BOOST_AUTO_TEST_CASE(array_ref_2D_from_vector) {
	std::vector<double> vec = {1, 2, 3, 4, 5, 6};
	multi::array_ref<double, 2> aref({2, 3}, vec.data());
	BOOST_REQUIRE( &aref[1][0] == &vec[3] );
}

BOOST_AUTO_TEST_CASE(array_ref_2D_from_vector_with_offset) {
	std::vector<double> vec = {
		1, 2, 3,
		4, 5, 6
	};
	multi::array_ref<double, 2> aref({multi::iextension(1, 3), multi::iextension(1, 4)}, vec.data());

	auto exts = aref.extensions();
	BOOST_REQUIRE( std::get<0>(exts) == multi::iextension(1, 3) );
	BOOST_REQUIRE( std::get<1>(exts).start()  == 1 );
	BOOST_REQUIRE( std::get<1>(exts).finish() == 4 );
	BOOST_REQUIRE( std::get<1>(exts) == multi::iextension(1, 4) );
	BOOST_REQUIRE( exts == decltype(exts)(multi::iextension(1, 3), multi::iextension(1, 4)) );

	BOOST_REQUIRE( &aref[1][1] == vec.data() );
}

BOOST_AUTO_TEST_CASE(array_2D_with_offset) {
	multi::array<double, 2> arr({multi::iextension(1, 3), multi::iextension(2, 5)}, 1.2);

	BOOST_REQUIRE( arr.extension().start()  == 1 );
	BOOST_REQUIRE( arr.extension().finish() == 3 );
}

BOOST_AUTO_TEST_CASE(array_ref_1D) {
	std::array<std::string, 5> arr = {{"a", "b", "c", "d", "e"}};

	multi::array_ref<std::string, 1>&& mar = *multi::array_ptr<std::string, 1>{&arr};
//  multi::Array<std::string(&)[1]> mar = *multi::Array<std::string(*)[1]>(&a);

	BOOST_REQUIRE(  extension(mar).first() == 0 );
	BOOST_REQUIRE(  extension(mar).last()  == 5 );

	auto&& mar1 = mar.reindexed(1);

	BOOST_REQUIRE( extension(mar1).size() == extension(mar).size() );

	BOOST_REQUIRE( mar1.extension() == extension(mar1) );
	BOOST_REQUIRE(  extension(mar1).first() == 1 );
	BOOST_REQUIRE(  mar1.extension().first() == 1 );
	BOOST_REQUIRE(  mar1.extension().last()  == 6 );
	BOOST_REQUIRE( *extension(mar1).begin() == 1 );

	BOOST_REQUIRE( size(mar1) == size(mar) );
	BOOST_REQUIRE( mar1.layout().extension().start() == 1 );
	BOOST_REQUIRE( extension(mar1).start() == 1 );
	BOOST_REQUIRE( &mar1[1]     == &arr[0] );  // NOLINT(readability-container-data-pointer) test access
	BOOST_REQUIRE(  mar1.base() == &arr[0] );  // NOLINT(readability-container-data-pointer) test access
	BOOST_REQUIRE(  mar1.base() ==  arr.data() );
}

BOOST_AUTO_TEST_CASE(array_ref_original_tests_carray) {
	double darr[4][5] = {{1., 2.}, {2., 3.}};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
	multi::array_ref<double, 2               >  ref(&darr[0][0], {4, 5});
	multi::array_ref<double, 2, double const*> cref(&darr[0][0], {4, 5});
	multi::array_ref<double const, 2> crefc(&darr[0][0], {4, 5});
	multi::array_cref<double, 2> ref2(&darr[0][0], {4, 5});

	BOOST_REQUIRE( &ref[1][2] == &cref[1][2] );
	BOOST_REQUIRE( &ref[1][2] == &crefc[1][2] );
	BOOST_REQUIRE( &ref[1][2] == &ref2[1][2] );

	ref[1][1] = 2.;

	double darr2[4][5] = {{1., 2.}, {2., 3.}};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type

	auto const& dd = static_cast<double const(&)[4][5]>(darr2);  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
	BOOST_REQUIRE( &(dd[1][2]) == &(darr2[1][2]) );
	BOOST_REQUIRE(( & ref[1].static_array_cast<double, double const*>()[1] == &ref[1][1] ));
	BOOST_REQUIRE(( &multi::static_array_cast<double, double const*>(ref[1])[1] == &ref[1][1] ));
}

BOOST_AUTO_TEST_CASE(array_ref_original_tests_const_carray) {
	double const d2D[4][5] = {{1., 2.}, {2., 3.}};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
	multi::array_ref<double, 2, const double*> d2Rce(&d2D[0][0], {4, 5});
	BOOST_REQUIRE( &d2Rce[2][3] == &d2D[2][3] );
	BOOST_REQUIRE( d2Rce.size() == 4 );
	BOOST_REQUIRE( num_elements(d2Rce) == 20 );
}

BOOST_AUTO_TEST_CASE(array_ref_original_tests_const_carray_string) {
	#if not defined(__circle_build__)  // circle 170 crashes https://github.com/seanbaxter/circle/issues/114
	std::string const dc3D[4][2][3] = {  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
		{{"A0a", "A0b", "A0c"}, {"A1a", "A1b", "A1c"}},
		{{"B0a", "B0b", "B0c"}, {"B1a", "B1b", "B1c"}},
		{{"C0a", "C0b", "C0c"}, {"C1a", "C1b", "C1c"}},
		{{"D0a", "D0b", "D0c"}, {"D1a", "D1b", "D1c"}},
	};
	multi::array_cref<std::string, 3> cref(&dc3D[0][0][0], {4, 2, 3});
	BOOST_REQUIRE( num_elements(cref) == 24 and cref[2][1][1] == "C1b" );
	auto const& A2 = cref.sliced(0, 3).rotated()[1].sliced(0, 2).unrotated();
	BOOST_REQUIRE( multi::rank<std::decay_t<decltype(A2)>>{} == 2 and num_elements(A2) == 6 );

	BOOST_REQUIRE( std::get<0>(sizes(A2)) == 3 and std::get<1>(sizes(A2)) == 2 );

	auto const& A3 = cref({0, 3}, 1, {0, 2});
	BOOST_REQUIRE( multi::rank<std::decay_t<decltype(A3)>>{} == 2 and num_elements(A3) == 6 );

	BOOST_REQUIRE( A2.layout()[2][1] == &A2[2][1] - A2.base() );
	BOOST_REQUIRE( A2.rotated().layout()[1][2] == &A2.rotated()[1][2] - A2.rotated().base() );
	#endif
}

BOOST_AUTO_TEST_CASE(array_ref_rebuild_2D) {
	double d2D[4][5] = {{1., 2.}, {2., 3.}};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
	multi::array_ref<double, 2> d2R(&d2D[0][0], {4, 5});
	auto&& d2B = d2R();
	auto&& d2B_ref = multi::ref(d2B.begin(), d2B.end());

	BOOST_REQUIRE( d2B.base()   == d2B_ref.base() );
	BOOST_REQUIRE( d2B.layout() == d2B_ref.layout() );
	BOOST_REQUIRE( &d2R() == &multi::ref(d2B.begin(), d2B.end()) );
}

BOOST_AUTO_TEST_CASE(array_ref_rebuild_1D) {
	double d1D[5] = {1., 2., 3., 4., 5.};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
	multi::array_ref<double, 1> d1R(&d1D[0], {5});
	auto&& d1B = d1R();
	auto&& d1B_ref = multi::ref(d1B.begin(), d1B.end());

	BOOST_REQUIRE( d1B.base()   == d1B_ref.base() );
	BOOST_REQUIRE( d1B.layout() == d1B_ref.layout() );
	BOOST_REQUIRE( &d1R() == &multi::ref(d1B.begin(), d1B.end()) );
}

BOOST_AUTO_TEST_CASE(array_ref_move_assigment_2D) {
	{
		multi::array<double, 2> arr ({5, 4}); std::iota(arr.elements().begin(), arr.elements().end(),  0.);
		multi::array<double, 2> arr2({5, 4}); std::iota(arr2.elements().begin(), arr2.elements().end(), 10.);

		multi::array_ref<double, 2>&& Aref{{5, 4}, arr.data_elements()};
		multi::array_ref<double, 2>&& Bref{{5, 4}, arr2.data_elements()};

		Bref = Aref;

		BOOST_REQUIRE( arr2 == arr );
	}
	{
		multi::array<double, 2> arr ({5, 4}); std::iota(arr.elements().begin(), arr.elements().end(),  0.);
		multi::array<double, 2> arr2({5, 4}); std::iota(arr2.elements().begin(), arr2.elements().end(), 10.);

		multi::array_ref<double, 2>&& ref2{{5, 4}, arr2.data_elements()};

		ref2 = multi::array_ref<double, 2>{{5, 4}, arr.data_elements()};

		BOOST_REQUIRE( arr2 == arr );
	}
	{
		multi::array<double, 2> arr ({5, 4}); std::iota(arr.elements().begin(), arr.elements().end(),  0.);
		multi::array<double, 2> arr2({5, 4}); std::iota(arr2.elements().begin(), arr2.elements().end(), 10.);

		multi::array_ref<double, 2>&& ref {{5, 4}, arr .data_elements()};
		multi::array_ref<double, 2>&& ref2{{5, 4}, arr2.data_elements()};

		ref2 = std::move(ref);

		BOOST_REQUIRE( arr2 == arr );
	}
}
