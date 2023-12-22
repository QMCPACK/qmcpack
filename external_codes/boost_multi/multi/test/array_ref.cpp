// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2023 Alfredo A. Correa

#include<boost/test/unit_test.hpp>

#include <multi/array.hpp>

#include <iostream>  // for std::cout
#include <numeric>  // for std::iota
#if defined(__cpp_lib_span) and (__cpp_lib_span >=  202002L)
#include <span>
#endif


namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(array_ref_from_carray) {
	// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays): test
	double arr[4][5] = {
		{ 0.0,  1.0,  2.0,  3.0,  4.0},
		{ 5.0,  6.0,  7.0,  8.0,  9.0},
		{10.0, 11.0, 12.0, 13.0, 14.0},
		{15.0, 16.0, 17.0, 18.0, 19.0},
	};

	multi::array_ptr<double, 2> const map{&arr};
	BOOST_REQUIRE( &map->operator[](1)[1] == &arr[1][1] );
	BOOST_REQUIRE( (*&arr)[1][1] == 6.0 );

	multi::array_ref<double, 2>&& mar = *map;

	BOOST_REQUIRE( &mar[1][1] == &arr[1][1] );

	mar[1][1] = 9.0;
	BOOST_REQUIRE( &mar[1][1] == &arr[1][1] );

	auto const& a_const = arr;
	//  double const(&a_const)[4][5] = a;
	BOOST_REQUIRE(&a_const[1][1] == &arr[1][1]);

	static_assert(decltype(mar(2, {1, 3}))::rank_v == 1);

	BOOST_REQUIRE( size(mar(2, {1, 3})) == 2 );
	BOOST_REQUIRE( &mar(2, {1, 3})[1] == &arr[2][2] );
}


BOOST_AUTO_TEST_CASE(array_ref_1D_reindexed) {
	using namespace std::string_literals;  // NOLINT(build/namespaces) for literal "string"s
	std::array<std::string, 5> stdarr = {{"a"s, "b"s, "c"s, "d"s, "e"s}};

	multi::array_ref<std::string, 1> mar = *multi::array_ptr<std::string, 1>(&stdarr);

	BOOST_REQUIRE( &mar[1] == &stdarr[1] );
	BOOST_REQUIRE( sizes(mar.reindexed(1)) == sizes(mar) );

	auto diff = &(mar.reindexed(1)[1]) - &mar[0];
	BOOST_REQUIRE( diff == 0 );

	BOOST_REQUIRE( &mar.blocked(2, 4)[2] == &mar[2] );
	for(auto idx : extension(mar.stenciled({2, 4}))) {
		BOOST_REQUIRE( &mar.stenciled({2, 4})[idx] == &mar[idx] );
	}

	multi::array<std::string, 1> arr({{2, 7}}, std::string{"xx"});  // std::string NOLINT(fuchsia-default-arguments-calls)
	BOOST_REQUIRE( size(arr) == 5 );
	BOOST_REQUIRE( extension(arr) == multi::iextension(2, 7) );
	arr[2] = "a";
	arr[3] = "b";
	arr[4] = "c";
	arr[5] = "d";
	arr[6] = "e";
	BOOST_REQUIRE( std::equal(arr.begin(), arr.end(), mar.begin(), mar.end()) );

	auto arrB = multi::array<std::string, 1>({"a", "b", "c", "d", "e"}).reindex(2);  // std::string NOLINT(fuchsia-default-arguments-calls)
	BOOST_REQUIRE( size(arrB) == 5 );
	BOOST_REQUIRE( arrB[2] == "a" );
	BOOST_REQUIRE( arrB[6] == "e" );
}

BOOST_AUTO_TEST_CASE(array_ref_of_nested_std_array_reindexed) {
	std::array<std::array<double, 5>, 4> arr = {{
		{{ 0.0,  1.0,  2.0,  3.0,  4.0}},
		{{ 5.0,  6.0,  7.0,  8.0,  9.0}},
		{{10.0, 11.0, 12.0, 13.0, 14.0}},
		{{15.0, 16.0, 17.0, 18.0, 19.0}}
	}};


	// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays): test type
	multi::array_ref<double, 2> mar = *multi::array_ptr<double, 2>(&arr);
	BOOST_REQUIRE( &mar[1][1] == &arr[1][1] );
}

BOOST_AUTO_TEST_CASE(array_ref_reindexed) {
	double (&&arr)[4][5] = {  // NOLINT(hicpp-avoid-c-arrays, modernize-avoid-c-arrays, cppcoreguidelines-avoid-c-arrays): test
		{ 0.0,  1.0,  2.0,  3.0,  4.0},
		{ 5.0,  6.0,  7.0,  8.0,  9.0},
		{10.0, 11.0, 12.0, 13.0, 14.0},
		{15.0, 16.0, 17.0, 18.0, 19.0}
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

	//  BOOST_REQUIRE( &(((mar<<1).reindexed(2)>>1).reindexed(1))[1][2] == &mar[0][0] );
	BOOST_REQUIRE( &mar.reindexed(1, 2)[1][2] == &mar[0][0] );

	BOOST_REQUIRE( &mar.reindexed(1)({1, 5})[1][0] == &mar[0][0] );

	BOOST_REQUIRE(( sizes(mar.stenciled({2, 4})) == decltype(sizes(mar.stenciled({2, 4}))){2, 5} ));
	BOOST_REQUIRE( &mar.stenciled({2, 4})[2][0] == &mar[2][0] );
	BOOST_REQUIRE( &mar.stenciled({2, 4}, {1, 3})[2][1] == &mar[2][1] );

	//  BOOST_REQUIRE( &mar[0][0] == mar.origin() ); // origin changed meaning in on 2020/Dec/16
	//  BOOST_REQUIRE( mar.base() == mar.origin() );

	//  BOOST_REQUIRE( mar.stenciled({2, 4}).origin() == mar.origin() );  // origin changed meaning in on 2020/Dec/16
	BOOST_REQUIRE( mar.stenciled({2, 4}).base()   != mar.base()   );

	BOOST_REQUIRE( &mar.stenciled({2, 4})[2][0] == mar.stenciled({2, 4}).base() );

	{
		multi::array<std::string, 2> arrB = {
			{"a", "b", "c", "d", "e"},  // std::string NOLINT(fuchsia-default-arguments-calls)
			{"f", "g", "h", "f", "g"},  // std::string NOLINT(fuchsia-default-arguments-calls)
			{"h", "i", "j", "k", "l"},  // std::string NOLINT(fuchsia-default-arguments-calls)
		};
		arrB.reindex(2);
		BOOST_REQUIRE( size(arrB) == 3 );
		BOOST_REQUIRE( arrB[2][0] == "a" );
	}
	{
		multi::array<std::string, 2> arrB = {
			{"a", "b", "c", "d", "e"},  // std::string NOLINT(fuchsia-default-arguments-calls)
			{"f", "g", "h", "f", "g"},  // std::string NOLINT(fuchsia-default-arguments-calls)
			{"h", "i", "j", "k", "l"},  // std::string NOLINT(fuchsia-default-arguments-calls)
		};
		arrB.reindex(2, 1);
		BOOST_REQUIRE( size(arrB) == 3 );
		BOOST_REQUIRE( arrB[2][1] == "a" );
	}
	{
		using namespace std::string_literals;  // NOLINT(build/namespaces) for literal "string"s
		multi::array<std::string, 2> arrB = (multi::array<std::string, 2>{
			{"a"s, "b"s, "c"s, "d"s, "e"s},
			{"f"s, "g"s, "h"s, "f"s, "g"s},
			{"h"s, "i"s, "j"s, "k"s, "l"s}
		});  // .reindex(2, 1);  // std::string NOLINT(fuchsia-default-arguments-calls)

		BOOST_REQUIRE( arrB.reindex(2).extension() == multi::iextension(2, 5) );
		auto exts = arrB.reindexed(2).extensions();

		multi::array<std::string, 2> const arrC(exts);
		BOOST_REQUIRE( size(arrC) == 3 );
		BOOST_REQUIRE( size(arrC) == size(arrB) );

		BOOST_REQUIRE( arrC.extension().start()  == 2 );
		BOOST_REQUIRE( arrC.extension().finish() == 5 );
	}
}

BOOST_AUTO_TEST_CASE(array_ref_with_stencil) {
	std::array<std::array<double, 5>, 4> arr = {{
		{{ 0.0,  1.0,  2.0,  3.0,  4.0}},
		{{ 5.0,  6.0,  7.0,  8.0,  9.0}},
		{{10.0, 11.0, 12.0, 13.0, 14.0}},
		{{15.0, 16.0, 17.0, 18.0, 19.0}}
	}};
	// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test type
	auto const& mar = *multi::array_ptr<double, 2>(&arr);
	BOOST_REQUIRE( mar.size() == 4 );

	// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test type
	multi::array<double, 2> ss = {
		{ 0.0, +1.0,  0.0},
		{+1.0, -4.0, +1.0},
		{ 0.0, +1.0,  0.0}
	};
	auto const& stencil = ss.reindexed(-1, -1);

	BOOST_REQUIRE( stencil.size() == 3 );
	BOOST_REQUIRE( &stencil[-1][-1] == stencil.base() );
}

BOOST_AUTO_TEST_CASE(array_ref_1D_from_vector) {
	std::vector<double> vec = {1.0, 2.0, 3.0};  // std::vector NOLINT(fuchsia-default-arguments-calls)
	multi::array_ref<double, 1> aref({{1, 3}}, vec.data());
	BOOST_REQUIRE( aref.extension() == multi::iextension(1, 3) );
	BOOST_REQUIRE( &aref[1] == vec.data() );
}

BOOST_AUTO_TEST_CASE(array_ref_2D_from_vector) {
	std::vector<double> vec = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};  // std::string NOLINT(fuchsia-default-arguments-calls)
	multi::array_ref<double, 2> aref({2, 3}, vec.data());
	BOOST_REQUIRE( &aref[1][0] == &vec[3] );
}

BOOST_AUTO_TEST_CASE(array_ref_2D_from_vector_with_offset) {
	std::vector<double> vec = {  // NOLINT(fuchsia-default-arguments-calls)
		1.0, 2.0, 3.0,
		4.0, 5.0, 6.0
	};
	multi::array_ref<double, 2> aref({multi::iextension(1, 3), multi::iextension(1, 4)}, vec.data());

	{
		auto exts = aref.extensions();
		BOOST_REQUIRE( std::get<0>(exts) == multi::iextension(1, 3) );
		BOOST_REQUIRE( std::get<1>(exts).start()  == 1 );
		BOOST_REQUIRE( std::get<1>(exts).finish() == 4 );
		BOOST_REQUIRE( std::get<1>(exts) == multi::iextension(1, 4) );
		BOOST_REQUIRE( exts == decltype(exts)(multi::iextension(1, 3), multi::iextension(1, 4)) );
	}
	{
		auto const exts = aref.extensions();
		BOOST_REQUIRE( std::get<0>(exts) == multi::iextension(1, 3) );
		BOOST_REQUIRE( std::get<1>(exts).start()  == 1 );
		BOOST_REQUIRE( std::get<1>(exts).finish() == 4 );
		BOOST_REQUIRE( std::get<1>(exts) == multi::iextension(1, 4) );
		BOOST_REQUIRE( exts == decltype(exts)(multi::iextension(1, 3), multi::iextension(1, 4)) );
	}
	{
		BOOST_REQUIRE( std::get<0>(aref.extensions()) == multi::iextension(1, 3) );
		BOOST_REQUIRE( std::get<1>(aref.extensions()).start()  == 1 );
		BOOST_REQUIRE( std::get<1>(aref.extensions()).finish() == 4 );
		BOOST_REQUIRE( std::get<1>(aref.extensions()) == multi::iextension(1, 4) );
		BOOST_REQUIRE( aref.extensions() == decltype(aref.extensions())(multi::iextension(1, 3), multi::iextension(1, 4)) );
	}
	{
		auto ss = aref.sizes();
		BOOST_REQUIRE( std::get<0>(ss) == 2 );
		BOOST_REQUIRE( std::get<1>(ss) == 3 );
		BOOST_REQUIRE( ss == decltype(ss)(2, 3) );
	}
	{
		auto [nn, mm] = aref.sizes();
		BOOST_REQUIRE( nn == 2 );
		BOOST_REQUIRE( mm == 3 );
	}
	{
		auto const ss = aref.sizes();
		BOOST_REQUIRE( std::get<0>(ss) == 2 );
		BOOST_REQUIRE( std::get<1>(ss) == 3 );
		BOOST_REQUIRE( ss == decltype(ss)(2, 3) );
	}
	{
		BOOST_REQUIRE( std::get<0>(aref.sizes()) == 2 );
		BOOST_REQUIRE( std::get<1>(aref.sizes()) == 3 );
		BOOST_REQUIRE( aref.sizes() == decltype(aref.sizes())(2, 3) );
	}
	{
		auto const ss = aref.sizes();
		using std::get;
		BOOST_REQUIRE( get<0>(ss) == 2 );
		BOOST_REQUIRE( get<1>(ss) == 3 );
		BOOST_REQUIRE( ss == decltype(ss)(2, 3) );
	}
	{
		using std::get;
		BOOST_REQUIRE( get<0>(aref.sizes()) == 2 );
		BOOST_REQUIRE( get<1>(aref.sizes()) == 3 );
		BOOST_REQUIRE( aref.sizes() == decltype(aref.sizes())(2, 3) );
	}
	#if __cplusplus >= 202002L  // GCC: use of function template name with no prior declaration in function call with explicit template arguments is a C++20 extension
	{
		auto const ss = aref.sizes();
		BOOST_REQUIRE( get<0>(ss) == 2 );
		BOOST_REQUIRE( get<1>(ss) == 3 );
		BOOST_REQUIRE( ss == decltype(ss)(2, 3) );
	}
	{
		BOOST_REQUIRE( get<0>(aref.sizes()) == 2 );
		BOOST_REQUIRE( get<1>(aref.sizes()) == 3 );
		BOOST_REQUIRE( aref.sizes() == decltype(aref.sizes())(2, 3) );
	}
	#endif

	BOOST_REQUIRE( &aref[1][1] == vec.data() );
}

BOOST_AUTO_TEST_CASE(array_2D_with_offset) {
	multi::array<double, 2> const arr({multi::iextension(1, 3), multi::iextension(2, 5)}, 1.2);

	BOOST_REQUIRE( arr.extension().start()  == 1 );
	BOOST_REQUIRE( arr.extension().finish() == 3 );
}

BOOST_AUTO_TEST_CASE(array_ref_1D) {
	std::array<std::string, 5> arr = {{"a", "b", "c", "d", "e"}};  // std::string NOLINT(fuchsia-default-arguments-calls)

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
	double darr[4][5] = {{1.0, 2.0}, {2.0, 3.0}};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
	multi::array_ref <double      , 2               >  ref (&darr[0][0], {4, 5});
	multi::array_ref <double      , 2, double const*> cref (&darr[0][0], {4, 5});
	multi::array_ref <double const, 2               > crefc(&darr[0][0], {4, 5});
	multi::array_cref<double      , 2               >  ref2(&darr[0][0], {4, 5});

	BOOST_REQUIRE( &ref[1][2] == &cref [1][2] );
	BOOST_REQUIRE( &ref[1][2] == &crefc[1][2] );
	BOOST_REQUIRE( &ref[1][2] == & ref2[1][2] );

	ref[1][1] = 2.0;

	double darr2[4][5] = {{1.0, 0.0}, {2.0, 3.0}};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
	darr2[1][0] = 2.0;

	// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
	auto const& dd = std::as_const(darr2);

	BOOST_REQUIRE( &(dd[1][2]) == &(darr2[1][2]) );
	BOOST_REQUIRE(( & ref[1].static_array_cast<double, double const*>()[1] == &ref[1][1] ));
	BOOST_REQUIRE(( &multi::static_array_cast<double, double const*>(ref[1])[1] == &ref[1][1] ));
}

BOOST_AUTO_TEST_CASE(array_ref_cast_carray) {
	double darr[2][2] = {{1.0, 2.0}, {2.0, 3.0}};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
	multi::array_ref <double, 2> ref(&darr[0][0], {2, 2});

	auto&& other_darr = static_cast<double(&)[2][2]>(ref);  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type

	double(& other_darr2)[2][2] = static_cast<double(&)[2][2]>(ref);  // NOLINT(hicpp-use-auto,modernize-use-auto,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
	double(& other_darr3)[2][2](ref);  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type

	BOOST_REQUIRE( &ref        [1][0] == &darr[1][0] );
	BOOST_REQUIRE( &other_darr [1][0] == &darr[1][0] );
	BOOST_REQUIRE( &other_darr2[1][0] == &darr[1][0] );
	BOOST_REQUIRE( &other_darr3[1][0] == &darr[1][0] );

	try {
		double(& other_darr4)[3][3](ref);  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type

		BOOST_REQUIRE( &other_darr4[1][0] == &darr[1][0] );
	} catch(...) {}
}

BOOST_AUTO_TEST_CASE(array_ref_original_tests_const_carray) {
	double const d2D[4][5] = {{1.0, 2.0}, {2.0, 3.0}};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
	multi::array_ref<double, 2, double const*> d2Rce(&d2D[0][0], {4, 5});

	BOOST_REQUIRE( &d2Rce[2][3] == &d2D[2][3] );
	BOOST_REQUIRE( d2Rce.size() == 4 );
	BOOST_REQUIRE( num_elements(d2Rce) == 20 );
}

BOOST_AUTO_TEST_CASE(array_ref_original_tests_const_carray_string) {
	std::string const dc3D[4][2][3] = {  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
		{{"A0a", "A0b", "A0c"}, {"A1a", "A1b", "A1c"}},  // std::string NOLINT(fuchsia-default-arguments-calls)
		{{"B0a", "B0b", "B0c"}, {"B1a", "B1b", "B1c"}},  // std::string NOLINT(fuchsia-default-arguments-calls)
		{{"C0a", "C0b", "C0c"}, {"C1a", "C1b", "C1c"}},  // std::string NOLINT(fuchsia-default-arguments-calls)
		{{"D0a", "D0b", "D0c"}, {"D1a", "D1b", "D1c"}},  // std::string NOLINT(fuchsia-default-arguments-calls)
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
}

BOOST_AUTO_TEST_CASE(array_ref_sizes_assingment) {
	multi::array_cref<int, 3> const cref(nullptr, {4, 2, 3});
	{
		auto [sizes1, sizes2, sizes3] = cref.sizes();
		BOOST_REQUIRE( sizes1 == 4 );
		BOOST_REQUIRE( sizes2 == 2 );
		BOOST_REQUIRE( sizes3 == 3 );
	}
	{
		auto sizes1 = std::get<0>(cref.sizes());
		auto sizes2 = std::get<1>(cref.sizes());
		auto sizes3 = std::get<2>(cref.sizes());
		BOOST_REQUIRE( sizes1 == 4 );
		BOOST_REQUIRE( sizes2 == 2 );
		BOOST_REQUIRE( sizes3 == 3 );
	}
	{
		multi::size_t sizes1, sizes2, sizes3;  // NOLINT(readability-isolate-declaration,cppcoreguidelines-init-variables) test a bad idiom
		multi::tie(sizes1, sizes2, sizes3) = cref.sizes();

		BOOST_REQUIRE( sizes1 == 4 );
		BOOST_REQUIRE( sizes2 == 2 );
		BOOST_REQUIRE( sizes3 == 3 );
	}
	{
		auto const [sizes1, sizes2, sizes3] = cref.sizes();

		BOOST_REQUIRE( sizes1 == 4 );
		BOOST_REQUIRE( sizes2 == 2 );
		BOOST_REQUIRE( sizes3 == 3 );
	}
	{
		// NOLINTNEXTLINE(runtime/int)
		long sizes1, sizes2, sizes3;  // NOLINT(google-runtime-int,readability-isolate-declaration,cppcoreguidelines-init-variables) test bad idiom
		multi::tie(sizes1, sizes2, sizes3) = cref.sizes();

		BOOST_REQUIRE( sizes1 == 4 );
		BOOST_REQUIRE( sizes2 == 2 );
		BOOST_REQUIRE( sizes3 == 3 );
	}
	{
		// NOLINTNEXTLINE(runtime/int)
		long long sizes1, sizes2, sizes3;  // NOLINT(google-runtime-int,readability-isolate-declaration,cppcoreguidelines-init-variables) test bad idiom
		multi::tie(sizes1, sizes2, sizes3) = cref.sizes();

		BOOST_REQUIRE( sizes1 == 4 );
		BOOST_REQUIRE( sizes2 == 2 );
		BOOST_REQUIRE( sizes3 == 3 );
	}
	{
		int64_t sizes1, sizes2, sizes3;  // NOLINT(readability-isolate-declaration,cppcoreguidelines-init-variables) test bad idiom
		multi::tie(sizes1, sizes2, sizes3) = cref.sizes();

		BOOST_REQUIRE( sizes1 == 4 );
		BOOST_REQUIRE( sizes2 == 2 );
		BOOST_REQUIRE( sizes3 == 3 );
	}
}

BOOST_AUTO_TEST_CASE(array_ref_rebuild_2D) {
	double d2D[4][5] = {{1.0, 2.0}, {2.0, 3.0}};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
	multi::array_ref<double, 2> d2R(&d2D[0][0], {4, 5});
	auto&& d2B = d2R();
	auto&& d2B_ref = multi::ref(d2B.begin(), d2B.end());

	BOOST_REQUIRE(  d2B[0][0]    ==  d2B_ref[0][0] );
	BOOST_REQUIRE( &d2B[0][0]    == &d2B_ref[0][0] );

	BOOST_REQUIRE(  d2B.base()   ==  d2B_ref.base() );
	BOOST_REQUIRE(  d2B.layout() ==  d2B_ref.layout() );

	BOOST_REQUIRE( &d2R() == &multi::ref(d2B.begin(), d2B.end()) );
}

BOOST_AUTO_TEST_CASE(array_ref_rebuild_1D) {
	double d1D[5] = {1.0, 2.0, 3.0, 4.0, 5.0};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
	multi::array_ref<double, 1> d1R(&d1D[0], {5});
	auto&& d1B     = d1R();
	auto&& d1B_ref = multi::ref(d1B.begin(), d1B.end());

	BOOST_REQUIRE( d1B.base()   == d1B_ref.base() );
	BOOST_REQUIRE( d1B.layout() == d1B_ref.layout() );
	BOOST_REQUIRE( &d1R() == &multi::ref(d1B.begin(), d1B.end()) );
}

BOOST_AUTO_TEST_CASE(array_ref_move_assigment_2D) {
	{
		multi::array<double, 2> arr({5, 4});
		std::iota(arr.elements().begin(), arr.elements().end(), 0.);

		multi::array<double, 2> arr2({5, 4});
		std::iota(arr2.elements().begin(), arr2.elements().end(), 10.);

		multi::array_ref<double, 2>&& Aref{{5, 4}, arr .data_elements()};
		multi::array_ref<double, 2>&& Bref{{5, 4}, arr2.data_elements()};

		Bref = Aref;

		BOOST_REQUIRE( arr2 == arr );
	}
	{
		multi::array<double, 2> arr({5, 4});
		std::iota(arr.elements().begin(), arr.elements().end(), 0.0);

		multi::array<double, 2> arr2({5, 4});
		std::iota(arr2.elements().begin(), arr2.elements().end(), 10.0);

		multi::array_ref<double, 2>&& ref2{{5, 4}, arr2.data_elements()};

		ref2 = multi::array_ref<double, 2>({5, 4}, arr.data_elements());

		BOOST_REQUIRE( arr2 == arr );
	}
	{
		multi::array<double, 2> arr({5, 4});
		std::iota(arr.elements().begin(), arr.elements().end(), 0.0);

		multi::array<double, 2> arr2({5, 4});
		std::iota(arr2.elements().begin(), arr2.elements().end(), 10.0);

		auto&& ref  = multi::array_ref<double, 2>({5, 4}, arr.data_elements());
		auto&& ref2 = multi::array_ref<double, 2>({5, 4}, arr2.data_elements());

		ref2 = std::move(ref);

		BOOST_REQUIRE( arr2 == arr );
	}
}

auto f1d5(double(&carr)[5]) -> double;  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
auto f1d5(double(&carr)[5]) -> double {  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
	return carr[1];
}

void f2d54(double(&carr)[5][4]);  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
void f2d54(double(&carr)[5][4]) {  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
	BOOST_REQUIRE(carr[0][1] == 1.0);
}

BOOST_AUTO_TEST_CASE(array_ref_conversion_1D) {
	multi::array<double, 1> arr({5}, double{});
	BOOST_REQUIRE( arr.size() == 5 );
	std::iota(arr.elements().begin(), arr.elements().end(), 0.0);

	{
		auto& carr = static_cast<double(&)[5]>(arr);  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
		BOOST_REQUIRE( &carr[3] == &arr[3] );

		BOOST_REQUIRE(f1d5(static_cast<double(&)[5]>(arr)) == 1.0);  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
	}
	{
		double(&carr)[5](arr);  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
		BOOST_REQUIRE( &carr[3] == &arr[3] );
	//  f1d5((double(&)[5])(arr));  // this will warn with -Wold-style-cast  NOLINT(google-readability-casting,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
	}
}

BOOST_AUTO_TEST_CASE(array_ref_conversion_2D) {
	multi::array<double, 2> arr({5, 4});
	std::iota(arr.elements().begin(), arr.elements().end(), 0.0);

	{
		auto& carr = static_cast<double(&)[5][4]>(arr);  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
		BOOST_REQUIRE( &carr[3][2] == &arr[3][2] );

		f2d54(static_cast<double(&)[5][4]>(arr));  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
	}
	{
		double(&carr)[5][4](arr);  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
		BOOST_REQUIRE( &carr[3][2] == &arr[3][2] );
	//  f2d54((double(&)[5][4])(arr));  // this will warn with -Wold-style-cast NOLINT(google-readability-casting,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
	}
}

BOOST_AUTO_TEST_CASE(as_span) {
	#if defined(__cpp_lib_span) and (__cpp_lib_span >=  202002L)
	auto printMe0 = [](std::span<int> rng) {
        std::cout << "rng.size(): " << rng.size() << '\n';  // (4)
		std::for_each(rng.begin(), rng.end(), [](auto const& elem) {std::cout << elem << ' ';});
        std::cout << "\n\n";
	};
	#endif

	auto printMe1 = [](multi::array_ref<int, 1> rng) {
        std::cout << "rng.size(): " << rng.size() << '\n';  // (4)
		std::for_each(rng.begin(), rng.end(), [](auto const& elem) {std::cout << elem << ' ';});
        std::cout << "\n\n";
	};

	auto printMe2 = [](multi::array_ptr<int, 1> ptr) {
        std::cout << "ptr->size(): " << ptr->size() << '\n';  // (4)
		std::for_each(ptr->begin(), ptr->end(), [](auto const& elem) {std::cout << elem << ' ';});
        std::cout << "\n\n";
	};

	#if defined(__cpp_lib_span) and (__cpp_lib_span >=  202002L)
    {
        int arr[] = {1, 2, 3, 4};
        printMe0(arr);

        std::vector vec = {1, 2, 3, 4, 5};
        printMe0(vec);

        std::array<int, 6> arr2 = {{1, 2, 3, 4, 5, 6}};
        printMe0(arr2);
    }
	#endif
    {
        int arr[] = {1, 2, 3, 4};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test c-arrays

        printMe1(multi::array_ref<int, 1>{arr});
        printMe1(                         arr );

        std::vector<int> vec = {1, 2, 3, 4, 5};  // NOLINT(fuchsia-default-arguments-calls)

        printMe1( *multi::array_ptr<int, 1>{vec.data(), 5} );

        std::array<int, 6> arr2 = {{1, 2, 3, 4, 5, 6}};

		printMe1(arr2);
		printMe1(*multi::array_ptr<int, 1>{ arr2.data(), {6} });


        multi::static_array<int, 1> marr({10}, 99);
        printMe1(*multi::array_ptr<int, 1>{marr.data_elements(), 10});

	    // printMe2(&marr);
    }
    {
        int arr[] = {1, 2, 3, 4};  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test c-arrays
        printMe2(multi::array_ptr<int, 1>{&arr});
        printMe2(                         &arr );

        std::vector<int> vec = {1, 2, 3, 4, 5};  // NOLINT(fuchsia-default-arguments-calls)
        printMe2( {vec.data(), 5} );

        std::array<int, 6> arr2 = {{1, 2, 3, 4, 5, 6}};
	//  printMe2(&arr2);  // this crashes clang-tidy
		printMe2({arr2.data(), {6} });

    //  multi::static_array<int, 1> marr({10}, 99);
    //  printMe2(&marr);  // TODO(correaa) make this work
    }
}

BOOST_AUTO_TEST_CASE(diagonal) {
	double (&&arr)[4][3] = {  // NOLINT(hicpp-avoid-c-arrays, modernize-avoid-c-arrays, cppcoreguidelines-avoid-c-arrays): test
		{ 0.0,  1.0,  2.0},
		{ 5.0,  1.0,  7.0},
		{10.0, 11.0,  2.0},
		{99.0, 99.0, 99.9}
	};

	// NOLINTNEXTLINE(hicpp-avoid-c-arrays, modernize-avoid-c-arrays, cppcoreguidelines-avoid-c-arrays): special type
	multi::array_ref<double, 2> mar = *multi::array_ptr<double, 2>(&arr);

	BOOST_REQUIRE( &mar({0, 3}, {0, 3}).diagonal()[0] == &arr[0][0] );
	BOOST_REQUIRE( &mar({0, 3}, {0, 3}).diagonal()[1] == &arr[1][1] );
	BOOST_REQUIRE( &mar({0, 3}, {0, 3}).diagonal()[2] == &arr[2][2] );

	auto sum = 0.0;
	for(auto const& aii : mar.diagonal() ) {sum += aii;}  // NOLINT(altera-unroll-loops) test for-range loop
	BOOST_REQUIRE( sum == mar[0][0] + mar[1][1] + mar[2][2]);
}
