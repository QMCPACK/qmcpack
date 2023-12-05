// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2021 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi layout"
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"
#include "multi/utility.hpp"

#include "multi/detail/tuple_zip.hpp"

//#include <boost/archive/xml_iarchive.hpp>
//#include <boost/archive/xml_oarchive.hpp>

#include<tuple>

namespace multi = boost::multi;

static auto second_finish(multi::extensions_t<3> exts) {
	return std::get<1>(exts).finish();
}

BOOST_AUTO_TEST_CASE(extensions_3D) {
	BOOST_REQUIRE( 20 == second_finish( multi::extensions_t<3>  { {0, 10}, {0, 20}, {0, 30} }  ) );
	BOOST_REQUIRE( 20 == second_finish( multi::extensions_t<3>( { {0, 10}, {0, 20}, {0, 30} } )) );
	BOOST_REQUIRE( 20 == second_finish(                         { {0, 10}, {0, 20}, {0, 30} }  ) );

	multi::extensions_t<3> exts{ {0, 10}, {0, 20}, {0, 30} };
	BOOST_REQUIRE( 20 == second_finish(exts                                                     ) );
}

//BOOST_AUTO_TEST_CASE(extensions_1D) {
//	BOOST_REQUIRE( multi::extensions_t<1>( { {0, 10} } ) == multi::extensions_t<1>( { {0, 10} } ) );
//}

#if 0
BOOST_AUTO_TEST_CASE(serialize_extensions) {
	multi::extensions_t<3> x{51, 52, 53};
	std::stringstream ss;
	{
		boost::archive::xml_oarchive xoa{ss};
		xoa<< BOOST_SERIALIZATION_NVP(x);
	}
	{
		std::cerr<< ss.str() << std::endl;
		boost::archive::xml_iarchive xia{ss};
		multi::extensions_t<3> x2{51, 52, 53};
		xia>> BOOST_SERIALIZATION_NVP(x2);
		BOOST_REQUIRE(x == x2);
	}
}
#endif

BOOST_AUTO_TEST_CASE(extensions_to_linear) {
	multi::extensions_t<3> exts{4, 5, 3};
	BOOST_REQUIRE( exts.to_linear(0, 0, 0) ==  0 );
	BOOST_REQUIRE( exts.to_linear(0, 0, 1) ==  1 );
	BOOST_REQUIRE( exts.to_linear(0, 0, 2) ==  2 );
	BOOST_REQUIRE( exts.to_linear(0, 1, 0) ==  3 );
	BOOST_REQUIRE( exts.to_linear(0, 1, 1) ==  4 );
	BOOST_REQUIRE( exts.to_linear(0, 1, 2) ==  5 );
	BOOST_REQUIRE( exts.to_linear(1, 0, 0) == 15 );

	for(int eye = 0; eye != 4; ++eye) {
		for(int jay = 0; jay != 5; ++jay) {
			for(int kay = 0; kay != 3; ++kay) {
				BOOST_REQUIRE(( exts.from_linear(exts.to_linear(eye, jay, kay)) == decltype(exts.from_linear(exts.to_linear(eye, jay, kay))){eye, jay, kay} ));
			}
		}
	}

	BOOST_REQUIRE( exts.to_linear(4, 0, 0) == exts.num_elements() );

	for(int idx = 0; idx != exts.num_elements(); ++idx) {
		BOOST_REQUIRE( std::apply([&](auto... indices) {return exts.to_linear(indices...);}, exts.from_linear(idx)) == idx );
	}
}

BOOST_AUTO_TEST_CASE(extensions_layout_to_linear) {
	multi::array<double, 3> arr({40, 50, 80});
	auto&& sub = arr({10, 30}, {20, 32}, {60, 75});

	for(int i = 0; i != 10; ++i) {
		for(int j = 0; j != 12; ++j) {
			for(int k = 0; k != 15; ++k) {
				BOOST_REQUIRE( &  sub.base()  [sub.layout()(i, j, k)] == &sub(i, j, k) );
				BOOST_REQUIRE( &*(sub.base() + sub.layout()(i, j, k)) == &sub(i, j, k) );
			}
		}
	}
}

BOOST_AUTO_TEST_CASE(extensions_layout_to_linear_2) {
	multi::array<double, 3> arr({40, 50, 80});
	auto&& sub = arr({10, 30}, {20, 32}, {60, 75});

	auto const& rot = sub.rotated();
	auto exts = rot.extensions();
	for(auto eye : std::get<0>(exts)) {
		for(auto jay : std::get<1>(exts)) {
			for(auto kay : std::get<2>(exts)) {
				BOOST_REQUIRE( &  rot.base()  [rot.layout()(eye, jay, kay)] == &rot(eye, jay, kay) );
				BOOST_REQUIRE( &*(rot.base() + rot.layout()(eye, jay, kay)) == &rot(eye, jay, kay) );
			}
		}
	}
}

BOOST_AUTO_TEST_CASE(linearize) {
	multi::array<double, 3> arr({10, 20, 30});

	BOOST_REQUIRE((  25 % extensions(arr) == decltype(  25 % extensions(arr)){0, 0, 25} ));
	BOOST_REQUIRE((  55 % extensions(arr) == decltype(  55 % extensions(arr))(0, 1, 25) ));
	BOOST_REQUIRE(( 655 % extensions(arr) == decltype( 655 % extensions(arr))(1, 1, 25) ));
	BOOST_REQUIRE((1255 % extensions(arr) == decltype(1255 % extensions(arr))(2, 1, 25) ));

	auto const point = arr.extensions().from_linear(655);
//  BOOST_REQUIRE( p == std::make_tuple(1, 1, 25) );
	using multi::detail::get;
	BOOST_REQUIRE( get<0>(point) ==  1 );
	BOOST_REQUIRE( get<1>(point) ==  1 );
	BOOST_REQUIRE( get<2>(point) == 25 );
}

BOOST_AUTO_TEST_CASE(layout_0) {
	multi::array<double, 3> arr(
#if defined(__INTEL_COMPILER) or (defined(__GNUC__) and (__GNUC__ < 6))
		multi::extensions_t<3>
#endif
		{51, 52, 53}
	);
	BOOST_REQUIRE( size(arr) == 51      ); BOOST_REQUIRE( arr.size() == 51       );
	BOOST_REQUIRE( size(arr[0]) == 52   ); BOOST_REQUIRE( arr[0].size() == 52    );
	BOOST_REQUIRE( size(arr[0][0]) == 53); BOOST_REQUIRE( arr[0][0].size() == 53 );
}

BOOST_AUTO_TEST_CASE(layout_1) {
	//NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays): testing feature
	double arr[50][50][50];
	using multi::size;
	BOOST_REQUIRE( size(arr) == 50 );

	using multi::extension;
	BOOST_REQUIRE(( extension(arr) == multi::index_extension{0, 50} ));
	BOOST_REQUIRE(( extension(arr) == multi::iextension{0, 50}      ));
	BOOST_REQUIRE(( extension(arr) == multi::irange{0, 50} ));
}

BOOST_AUTO_TEST_CASE(layout_2) {
	std::array<std::array<std::array<double, 50>, 50>, 50> arr{};
	using multi::size;
	BOOST_REQUIRE( size(arr) == 50 );

	using multi::extension;
	BOOST_REQUIRE(( extension(arr) == multi::index_extension{0, 50} ));
	BOOST_REQUIRE(( extension(arr) == multi::iextension{0, 50}      ));
	BOOST_REQUIRE(( extension(arr) == multi::irange{0, 50}          ));
}

BOOST_AUTO_TEST_CASE(layout_3) {
	multi::array<double, 2> arr(
//#if defined(__INTEL_COMPILER) or (defined(__GNUC__) and (__GNUC__ < 6))
//		multi::extensions_t<2>
//#endif
		{50, 50}
	);
	BOOST_REQUIRE( size(arr) == 50 ); BOOST_REQUIRE( arr.size() == 50 );
	BOOST_REQUIRE( arr[0].sliced(10, 20).size() == 10 );
	BOOST_REQUIRE( size(arr[0].sliced(10, 20))  == 10 );

	static_assert( decltype(arr(0, {10, 20}))::rank_v  == 1 , "!");

	BOOST_REQUIRE( size(arr(0, {10, 20})) == 10 );

	BOOST_REQUIRE(      arr.layout() == arr.layout()  );
	BOOST_REQUIRE( not (arr.layout() <  arr.layout()) );
}

BOOST_AUTO_TEST_CASE(layout) {
{
	multi::array<double, 2> A2 = {
		{1., 2., 3.},
		{4., 5., 6.},
		{7., 8., 9.}
	};

	BOOST_REQUIRE( size(A2) == 3 );

	multi::array<int, 2> B2(
#if defined(__INTEL_COMPILER) or (defined(__GNUC__) and (__GNUC__ < 6))
		multi::extensions_t<2>
#endif
		{4, 4}
	);
	BOOST_REQUIRE( size(B2) == 4 );
	B2[3][3] = 99.;

	auto B2copy =+ B2({0, 2}, {0, 2});

	BOOST_REQUIRE( &B2copy[1][1] != &B2({0, 2}, {0, 2})[1][1] );

	std::array<std::array<decltype(B2({0, 2}, {0, 2})), 2>, 2> B2blk = {{
		{{ B2({0, 2}, {0, 2}), B2({0, 2}, {2, 4}) }},
		{{ B2({2, 4}, {0, 2}), B2({2, 4}, {2, 4}) }}
	}};

	BOOST_REQUIRE( &B2blk[1][1][1][1] == &B2[3][3] );
}
{
	// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays, hicpp-avoid-c-arrays,-warnings-as-errors, modernize-avoid-c-arrays,-warnings-as-errors): test
	double arr[3][4][5] = {};
	using multi::dimensionality;
	static_assert(dimensionality(arr)==3, "!");
	using multi::extensions;
	auto xA = extensions(arr);

	BOOST_REQUIRE( size(std::get<0>(xA)) == 3 );
	BOOST_REQUIRE( size(std::get<1>(xA)) == 4 );
	BOOST_REQUIRE( size(std::get<2>(xA)) == 5 );

	static_assert( multi::stride(arr)       == 20 , "!" );

//  static_assert( multi::stride(arr)       ==  5 , "!" );
	static_assert( multi::stride(arr[1])    ==  5 , "!" );
	static_assert( multi::stride(arr[0][0]) ==  1 , "!" );

	multi::array<double, 3> AA({3, 4, 5});
	using multi::layout;
	BOOST_REQUIRE( layout(AA) == layout(arr) );

	BOOST_REQUIRE( AA     .stride() == 20 );
}
{
	std::array<std::array<std::array<double, 5>, 4>, 3> arr = {};
#if defined(__circle_build__)  // circle doesn't see dimensionality as a constexpr "cannot access value of A at compile time;"
	       assert( multi::dimensionality(arr) == 3 );
#else  // other compilers ok
	static_assert( multi::dimensionality(arr) == 3 );
#endif

	using multi::extensions;
	auto xA = extensions(arr);
	using std::get;
	BOOST_REQUIRE( size(std::get<0>(xA)) == 3 );
	BOOST_REQUIRE( size(std::get<1>(xA)) == 4 );
	BOOST_REQUIRE( size(std::get<2>(xA)) == 5 );

	multi::array<double, 3> AA({3, 4, 5});
	using multi::layout;
	BOOST_REQUIRE( layout(AA) == layout(arr) );

	BOOST_REQUIRE( AA.stride() == 20 );

#if defined(__circle_build__)  // circle doesn't recognize this as a constexpr "cannot access value of A at compile time;"
	       assert( multi::stride(arr) == 20);
#else  // other compilers ok
	static_assert( multi::stride(arr) == 20);
#endif

	BOOST_REQUIRE( multi::stride(arr[0])    == 5 );
	BOOST_REQUIRE( multi::stride(arr[1])    == 5 );
	BOOST_REQUIRE( multi::stride(arr[0][0]) == 1 );
//		assert( stride(arr) == 20 );
//		assert( stride(arr[0]) == 20 );
}
{
	multi::array<double, 2> B2 = {
		{1.},
		{2.},
		{3.}
	};
	BOOST_REQUIRE( size(B2) == 3 );
	BOOST_REQUIRE( size(rotated(B2)) == 1 ); BOOST_REQUIRE( size(B2[0]) == 1);
	BOOST_REQUIRE( B2   .stride() == 1 );
	BOOST_REQUIRE( B2[0].stride() == 1 );
}
}

BOOST_AUTO_TEST_CASE(multi_layout_with_offset) {
	{
		multi::layout_t<1> l1(multi::iextension(2, 5));
		BOOST_REQUIRE( l1.extension().start()  == 2 );
		BOOST_REQUIRE( l1.extension().finish() == 5 );
	}
	{
		boost::multi::layout_t<2>::extensions_type exts{
			multi::iextension(2, 5),
			multi::iextension(0, 5)
		};
		multi::layout_t<2> l2(exts);
		BOOST_REQUIRE( l2.extension().start()  == std::get<0>(exts).start()  );
		BOOST_REQUIRE( l2.extension().finish() == std::get<0>(exts).finish() );
	}
	 {
		multi::layout_t<2> l2({multi::iextension(0, 3), multi::iextension(2, 7)});
		BOOST_REQUIRE( std::get<1>(l2.extensions()).start()  == 2 );
		BOOST_REQUIRE( std::get<1>(l2.extensions()).finish() == 7 );
	}
}

BOOST_AUTO_TEST_CASE(multi_layout_part1) {
{
	multi::layout_t<0> lyt;
	static_assert( decltype(lyt)::rank_v==0 , "!");
	BOOST_REQUIRE( num_elements(lyt) == 1 );
}{
	multi::iextensions<0> exts{};
	multi::layout_t<0> lyt(exts);
	BOOST_REQUIRE(lyt.num_elements() == 1);
}{  multi::layout_t<1> lyt{};
	static_assert( decltype(lyt)::rank_v == 1 , "!");
	BOOST_REQUIRE( num_elements(lyt) == 0 );
	BOOST_REQUIRE( size(lyt) == 0 );
	BOOST_REQUIRE( size(extension(lyt))==0 );
	BOOST_REQUIRE( stride(lyt)!=0 );
	BOOST_REQUIRE( is_empty(lyt) );
}{
	multi::layout_t<2> lyt({2, 10});
	static_assert( decltype(lyt)::rank_v == 2 , "!");
	BOOST_REQUIRE( num_elements(lyt) == 20 );
	BOOST_REQUIRE( size(lyt) == 2 );
	BOOST_REQUIRE( size(extension(lyt))==2 );
	BOOST_REQUIRE( stride(lyt)==10 );
	BOOST_REQUIRE( not is_empty(lyt) );
} {
	multi::layout_t<1> lyt(multi::iextensions<1>{20});
	static_assert( decltype(lyt)::rank_v == 1 , "!");
	BOOST_REQUIRE( num_elements(lyt) == 20 );
	BOOST_REQUIRE( size(lyt) == 20 );
	BOOST_REQUIRE( stride(lyt) == 1 );
}
}

BOOST_AUTO_TEST_CASE(multi_layout_part2) {
{
	multi::layout_t<1> lyt(multi::iextensions<1>{1});
	static_assert( decltype(lyt)::rank_v ==1 , "!");
	BOOST_REQUIRE( num_elements(lyt) == 1 );
	BOOST_REQUIRE( size(lyt) == 1 );
	BOOST_REQUIRE( stride(lyt) == 1 );
} {
	multi::layout_t<2> lyt({1, 10});
	static_assert( decltype(lyt)::rank_v ==2 , "!");
	BOOST_REQUIRE( num_elements(lyt) == 10 );
	BOOST_REQUIRE( size(lyt) == 1);
	BOOST_REQUIRE( not is_empty(lyt) );
	BOOST_REQUIRE( size(extension(lyt))==1 );
	BOOST_REQUIRE( stride(lyt)== 10 );//std::numeric_limits<std::ptrdiff_t>::max() );
	using std::get;
	BOOST_REQUIRE( get<0>(strides(lyt)) == 10);
	BOOST_REQUIRE( get<1>(strides(lyt)) == 1 );
}
}

BOOST_AUTO_TEST_CASE(multi_layout_part3) {
{
	multi::layout_t<2> lyt({10, 1});
	static_assert( decltype(lyt)::rank_v ==2 , "!");
	BOOST_REQUIRE( num_elements(lyt) == 10 );
	BOOST_REQUIRE( size(lyt) == 10 );
	using std::get;
	BOOST_REQUIRE( get<0>(strides(lyt)) == 1 );
	BOOST_REQUIRE( get<1>(strides(lyt)) == 1 );
}{  multi::layout_t<2> lyt{};
	BOOST_REQUIRE( dimensionality(lyt)==2 );
	BOOST_REQUIRE( num_elements(lyt) == 0 );
	BOOST_REQUIRE( size(lyt) == 0 );
	BOOST_REQUIRE( size(extension(lyt))==0 );
	BOOST_REQUIRE( stride(lyt)!=0 );
	BOOST_REQUIRE( is_empty(lyt) );
}{  multi::layout_t<3> lyt{}; BOOST_REQUIRE( num_elements(lyt) == 0 );
}{	multi::layout_t<3> lyt({{0, 10}, {0, 10}, {0, 10}}); BOOST_REQUIRE( num_elements(lyt) == 1000 );
}{	multi::layout_t<3> lyt({{10}, {10}, {10}}); BOOST_REQUIRE( num_elements(lyt) == 1000 );
}{	multi::layout_t<3> lyt({10, 10, 10}); BOOST_REQUIRE( num_elements(lyt) == 1000 );
}{	multi::layout_t<3> lyt({multi::index_extension{0, 10}, {0, 10}, {0, 10}}); BOOST_REQUIRE( num_elements(lyt) == 1000 );
}{	multi::layout_t<3> lyt(multi::layout_t<3>::extensions_type{{0, 10}, {0, 10}, {0, 10}}); BOOST_REQUIRE( num_elements(lyt) == 1000 );
}
}

BOOST_AUTO_TEST_CASE(layout_to_offset) {
	multi::layout_t<3> lyt({10, 20, 30});
	multi::array<double, 3> arr({10, 20, 30});
	BOOST_REQUIRE( lyt[0][0][0] == &arr[0][0][0] - arr.data_elements() );
	BOOST_REQUIRE( lyt[0][0][1] == &arr[0][0][1] - arr.data_elements() );
	BOOST_REQUIRE( lyt[0][0][2] == &arr[0][0][2] - arr.data_elements() );

	BOOST_TEST_REQUIRE( lyt[0][1][2] == &arr[0][1][2] - arr.data_elements() );
	BOOST_TEST_REQUIRE( lyt[3][1][2] == &arr[3][1][2] - arr.data_elements() );
}

BOOST_AUTO_TEST_CASE(layout_to_offset_sub) {
	multi::array<double, 3> arr({10, 20, 30});
	auto&& sub = arr({2, 6}, {4, 8}, {10, 20});
	auto const lyt = sub.layout();
	BOOST_REQUIRE( lyt[0][0][0] == &sub[0][0][0] - base(sub) );
	BOOST_REQUIRE( lyt[0][0][1] == &sub[0][0][1] - base(sub) );
	BOOST_REQUIRE( lyt[0][0][2] == &sub[0][0][2] - base(sub) );
	BOOST_REQUIRE( lyt[0][1][2] == &sub[0][1][2] - base(sub) );
	BOOST_REQUIRE( lyt[3][1][2] == &sub[3][1][2] - base(sub) );
}

BOOST_AUTO_TEST_CASE(continued_part1) {
{
	multi::layout_t<3> lyt(multi::layout_t<3>::extensions_type{{0, 10}, {0, 10}, {0, 10}});
	BOOST_REQUIRE( num_elements(lyt) == 1000);
}
{	multi::layout_t<3> lyt({multi::iextension{0, 10}, multi::iextension{0, 10}, multi::iextension{0, 10}}); BOOST_REQUIRE(lyt.num_elements() == 1000);
}{	multi::layout_t<3> lyt({multi::iextension{10}, multi::iextension{10}, multi::iextension{10}}); BOOST_REQUIRE( num_elements(lyt) == 1000);
}{	multi::layout_t<3> lyt({10, 10, multi::iextension{10}}); BOOST_REQUIRE( num_elements(lyt) == 1000 );
}{
	multi::layout_t<1> lyt;
	BOOST_REQUIRE( size(lyt) == 0 );
}{
	multi::layout_t<1> lyt({{0, 10}});
	BOOST_REQUIRE( size(lyt) == 10 );
	BOOST_REQUIRE( extension(lyt).start () ==  0 );
	BOOST_REQUIRE( extension(lyt).finish() == 10 );

	lyt.reindex(1);
	BOOST_REQUIRE( size(lyt) == 10 );
	BOOST_REQUIRE( extension(lyt).start () ==  1 );
	BOOST_REQUIRE( extension(lyt).finish() == 11 );
}{
	multi::layout_t<2> lyt;
	BOOST_REQUIRE( size(lyt) == 0 );
}
{
	multi::layout_t<2> lyt({{0, 10}, {0, 20}});
	BOOST_REQUIRE( size(lyt) == 10 );
	BOOST_REQUIRE( extension(lyt).start () ==  0 );
	BOOST_REQUIRE( extension(lyt).finish() == 10 );

	lyt.reindex(1);
	BOOST_REQUIRE( extension(lyt).start () ==  1 );
	BOOST_REQUIRE( extension(lyt).finish() == 11 );

	lyt.rotate().reindex(3).unrotate();
	BOOST_TEST_REQUIRE( extension(lyt).start () ==  1 );
	BOOST_TEST_REQUIRE( extension(lyt).finish() == 11 );

	BOOST_TEST_REQUIRE( std::get<0>(extensions(lyt)).start () == 1 );
	BOOST_TEST_REQUIRE( std::get<1>(extensions(lyt)).start () == 3 );
	BOOST_TEST_REQUIRE( std::get<1>(extensions(lyt)).finish() == 23 );
}
}

BOOST_AUTO_TEST_CASE(continued_part2) {
	multi::layout_t<3> lyt({{0, 10}, {0, 20}, {0, 30}});

	BOOST_REQUIRE( not lyt.empty() );

	BOOST_REQUIRE( stride(lyt) == lyt.stride() );
	BOOST_REQUIRE( offset(lyt) == lyt.offset() );
	BOOST_REQUIRE( nelems(lyt) == lyt.nelems() );

	BOOST_REQUIRE( stride(lyt) == 20*30L );
	BOOST_REQUIRE( offset(lyt) == 0 );
	BOOST_REQUIRE( nelems(lyt) == 10*20L*30L );

	BOOST_REQUIRE( lyt.stride() == stride(lyt) );
	BOOST_REQUIRE( lyt.offset() == offset(lyt) );
	BOOST_REQUIRE( lyt.nelems() == nelems(lyt) );

	using boost::multi::detail::get;
	BOOST_REQUIRE( get<1>(lyt.strides()) == 30     );
	BOOST_REQUIRE( get<1>(lyt.offsets()) ==  0     );
	BOOST_REQUIRE( get<1>(lyt.nelemss()) == 20*30L );

	BOOST_REQUIRE( get<2>(lyt.strides()) ==  1 );
	BOOST_REQUIRE( get<2>(lyt.offsets()) ==  0 );
	BOOST_REQUIRE( get<2>(lyt.nelemss()) == 30 );
}

BOOST_AUTO_TEST_CASE(continued_part3) {
	multi::layout_t<3> lyt({{0, 10}, {0, 20}, {0, 30}});

	BOOST_REQUIRE( lyt.num_elements() == num_elements(lyt) );
	BOOST_REQUIRE( lyt.size() == size(lyt) );
	BOOST_REQUIRE( lyt.extension() == extension(lyt) );

	BOOST_REQUIRE( num_elements(lyt) == 10*20L*30L );
	BOOST_REQUIRE( size(lyt) == 10 );
	BOOST_REQUIRE( extension(lyt).first() == 0 );
	BOOST_REQUIRE( extension(lyt).last() == 10 );

	BOOST_REQUIRE( std::get<0>(lyt.extensions()) == lyt.extension() );

	boost::multi::extensions_t<2> exts2;

	using boost::multi::detail::get;
	using std::get;

	BOOST_REQUIRE( get<0>(exts2).is_empty() );

//	BOOST_REQUIRE( std::get<0>(L.sizes()) == L.size(0) );
//	BOOST_REQUIRE( std::get<0>(L.extensions()) == L.extension(0) );

	BOOST_REQUIRE(( get<0>(lyt.extensions()) == multi::index_extension{0, 10} ));

	BOOST_REQUIRE( get<0>(lyt.extensions()).first() ==  0 );
	BOOST_REQUIRE( get<0>(lyt.extensions()).last()  == 10 );

//  BOOST_REQUIRE( L.size(1) == 20 );
	BOOST_REQUIRE( get<1>(lyt.extensions()).first() ==  0 );
	BOOST_REQUIRE( get<1>(lyt.extensions()).last()  == 20 );

//  BOOST_REQUIRE( L.size(2) == 30 );
	BOOST_REQUIRE( get<2>(lyt.extensions()).first() ==  0 );
	BOOST_REQUIRE( get<2>(lyt.extensions()).last()  == 30 );

	using std::get;
	BOOST_REQUIRE( get<0>(strides(lyt)) == lyt.stride() );

	auto const& strides = lyt.strides();
	BOOST_REQUIRE( get<0>(strides) == lyt.stride() );
}

BOOST_AUTO_TEST_CASE(continued) {
{
	multi::layout_t<3> lyt;
	BOOST_REQUIRE( size(lyt) == 0 );
}
{
	multi::layout_t<3> lyt( {{0, 10}, {0, 20}, {0, 30}} );
	BOOST_REQUIRE( stride(lyt) == 20*30L );
}
{
	multi::layout_t<1> lyt({{0, 10}});
	BOOST_REQUIRE( extension(lyt).first() == 0 );
	BOOST_REQUIRE( extension(lyt).last() == 10 );
}
{
	multi::layout_t<1> lyt({{8, 18}});
	BOOST_REQUIRE( extension(lyt).first() == 8 );
	BOOST_REQUIRE( extension(lyt).last() == 18 );
}
{
	multi::layout_t<2> lyt({{0, 10}, {0, 20}});
	BOOST_REQUIRE( extension(lyt).first() == 0 );
	BOOST_REQUIRE( extension(lyt).last() == 10 );
}
{
	multi::layout_t<2> lyt( {{0, 10}, {11, 31}} );
	BOOST_REQUIRE( size(lyt) == 10   );
	BOOST_REQUIRE( stride(lyt) == 20 );
	BOOST_REQUIRE( offset(lyt) == 0 );
}
{
	multi::layout_t<2> lyt( {{8, 18}, {0, 20}} );
	BOOST_REQUIRE( size(lyt) == 10 );
	BOOST_REQUIRE( stride(lyt) == 20 );
}
{
	multi::layout_t<3> lyt({{0, 3}, {0, 5}, {10, 17}});
	BOOST_REQUIRE( stride(lyt) == 5*7L );
	BOOST_REQUIRE( stride(lyt.sub().sub()) == 1 );
}
{
	multi::layout_t<3> lyt({{0, 10}, {0, 20}, {0, 30}});
	BOOST_REQUIRE( stride(lyt) == 20*30L );
	BOOST_REQUIRE( offset(lyt) == 0 );
	BOOST_REQUIRE( nelems(lyt) == 10*20L*30L );
}
{
	multi::layout_t<3> lyt({{10, 20}, {10, 30}, {10, 40}});
	BOOST_REQUIRE( stride(lyt) == 20*30L );
}
{
	auto const ttt = boost::multi::tuple<int, int, int>{1, 2, 3};
	auto const arr = std::apply([](auto... elems) {return std::array<int, 3>{{elems...}};}, ttt);
	BOOST_REQUIRE(arr[1] == 2);
}
}

//BOOST_AUTO_TEST_CASE(tuple_zip_test) {  // TODO(correaa) make it work
//	auto t1 = std::make_tuple( 1,  2,  3);
//	auto t2 = std::make_tuple(10, 20, 30);
//	auto t3 = std::make_tuple(std::string{"10"}, std::string{"20"}, std::string{"30"});
//	auto t123 = boost::multi::detail::tuple_zip(t1, t2, t3);
//	BOOST_REQUIRE( std::get<2>(std::get<0>(t123)) == std::string{"10"} );
//}

BOOST_AUTO_TEST_CASE(extensions_from_linear_1d) {
	multi::extensions_t<1> exts{11};

	auto ijk = exts.from_linear(9);

	using multi::detail::get;
	BOOST_TEST_REQUIRE( get<0>(ijk) == 9 );

	multi::layout_t<1> lyt{exts};
	BOOST_TEST_REQUIRE( lyt[get<0>(ijk)] == 9 );
	BOOST_TEST_REQUIRE( lyt(get<0>(ijk)) == 9 );

	BOOST_TEST_REQUIRE( lyt(std::get<0>(lyt.extensions().from_linear(9))) == 9 );

	BOOST_TEST_REQUIRE( std::apply(lyt, lyt.extensions().from_linear(9)) == 9 );
}

BOOST_AUTO_TEST_CASE(extensions_from_linear_2d_structured_binding) {
	multi::extensions_t<2> exts{3, 5};
	auto [eye, jay] = exts.from_linear(7);

	BOOST_TEST_REQUIRE( eye == 1 );
	BOOST_TEST_REQUIRE( jay == 2 );
//	BOOST_TEST_REQUIRE( std::apply(l, l.extensions().from_linear(9)) == 9 );
}

BOOST_AUTO_TEST_CASE(extensions_from_linear_2d_std_get) {
	multi::extensions_t<2> exts{3, 5};
	auto eye = std::get<0>(exts.from_linear(7));
	auto jay = std::get<1>(exts.from_linear(7));
	BOOST_TEST_REQUIRE( eye == 1 );
	BOOST_TEST_REQUIRE( jay == 2 );
}

BOOST_AUTO_TEST_CASE(extensions_from_linear_2d_std_get_using) {
	multi::extensions_t<2> exts{3, 5};
	using std::get;
	auto const eye = get<0>(exts.from_linear(7));
	auto const jay = get<1>(exts.from_linear(7));
	BOOST_TEST_REQUIRE( eye == 1 );
	BOOST_TEST_REQUIRE( jay == 2 );
}

BOOST_AUTO_TEST_CASE(extensions_from_linear_2d_get_using) {
	multi::extensions_t<2> exts{3, 5};

	using multi::detail::get;

	auto eye = get<0>(exts.from_linear(7));
	auto jay = get<1>(exts.from_linear(7));
	BOOST_TEST_REQUIRE( eye == 1 );
	BOOST_TEST_REQUIRE( jay == 2 );
}

BOOST_AUTO_TEST_CASE(extensions_from_linear_2d) {
	multi::extensions_t<2> exts{3, 5};

	auto ij = exts.from_linear(7);

	using multi::detail::get;

	BOOST_TEST_REQUIRE( get<0>(ij) == 1 );
	BOOST_TEST_REQUIRE( get<1>(ij) == 2 );

	multi::layout_t<2> lyt{exts};
	BOOST_TEST_REQUIRE( lyt[get<0>(ij)][get<1>(ij)] == 7 );
}

BOOST_AUTO_TEST_CASE(extensions_from_linear_3d_std_get) {
	multi::extensions_t<3> exts{11, 13, 17};

	BOOST_TEST_REQUIRE( std::get<0>(exts.from_linear( 0)) ==  0 );
	BOOST_TEST_REQUIRE( std::get<1>(exts.from_linear( 0)) ==  0 );
	BOOST_TEST_REQUIRE( std::get<2>(exts.from_linear( 0)) ==  0 );

	BOOST_TEST_REQUIRE( std::get<0>(exts.from_linear( 1)) ==  0 );
	BOOST_TEST_REQUIRE( std::get<1>(exts.from_linear( 1)) ==  0 );
	BOOST_TEST_REQUIRE( std::get<2>(exts.from_linear( 1)) ==  1 );

	BOOST_TEST_REQUIRE( std::get<0>(exts.from_linear(16)) ==  0 );
	BOOST_TEST_REQUIRE( std::get<1>(exts.from_linear(16)) ==  0 );
	BOOST_TEST_REQUIRE( std::get<2>(exts.from_linear(16)) == 16 );

	BOOST_TEST_REQUIRE( std::get<0>(exts.from_linear(17)) ==  0 );
	BOOST_TEST_REQUIRE( std::get<1>(exts.from_linear(17)) ==  1 );
	BOOST_TEST_REQUIRE( std::get<2>(exts.from_linear(17)) ==  0 );

	BOOST_TEST_REQUIRE( std::get<0>(exts.from_linear(18)) ==  0 );
	BOOST_TEST_REQUIRE( std::get<1>(exts.from_linear(18)) ==  1 );
	BOOST_TEST_REQUIRE( std::get<2>(exts.from_linear(18)) ==  1 );

	multi::layout_t<3> lyt{exts};
	using std::get;
	BOOST_TEST_REQUIRE( lyt[get<0>(exts.from_linear(19))][get<1>(exts.from_linear(19))][get<2>(exts.from_linear(19))] == 19 );
	BOOST_TEST_REQUIRE( lyt(get<0>(exts.from_linear(19)), get<1>(exts.from_linear(19)), get<2>(exts.from_linear(19))) == 19 );
}

BOOST_AUTO_TEST_CASE(extensions_from_linear_3d_std_get_using) {
	multi::extensions_t<3> exts{11, 13, 17};

	using std::get;

	BOOST_TEST_REQUIRE( get<0>(exts.from_linear( 0)) ==  0 );
	BOOST_TEST_REQUIRE( get<1>(exts.from_linear( 0)) ==  0 );
	BOOST_TEST_REQUIRE( get<2>(exts.from_linear( 0)) ==  0 );

	BOOST_TEST_REQUIRE( get<0>(exts.from_linear( 1)) ==  0 );
	BOOST_TEST_REQUIRE( get<1>(exts.from_linear( 1)) ==  0 );
	BOOST_TEST_REQUIRE( get<2>(exts.from_linear( 1)) ==  1 );

	BOOST_TEST_REQUIRE( get<0>(exts.from_linear(16)) ==  0 );
	BOOST_TEST_REQUIRE( get<1>(exts.from_linear(16)) ==  0 );
	BOOST_TEST_REQUIRE( get<2>(exts.from_linear(16)) == 16 );

	BOOST_TEST_REQUIRE( get<0>(exts.from_linear(17)) ==  0 );
	BOOST_TEST_REQUIRE( get<1>(exts.from_linear(17)) ==  1 );
	BOOST_TEST_REQUIRE( get<2>(exts.from_linear(17)) ==  0 );

	BOOST_TEST_REQUIRE( get<0>(exts.from_linear(18)) ==  0 );
	BOOST_TEST_REQUIRE( get<1>(exts.from_linear(18)) ==  1 );
	BOOST_TEST_REQUIRE( get<2>(exts.from_linear(18)) ==  1 );

	BOOST_TEST_REQUIRE( get<0>(exts.from_linear(19)) ==  0 );
	BOOST_TEST_REQUIRE( get<1>(exts.from_linear(19)) ==  1 );
	BOOST_TEST_REQUIRE( get<2>(exts.from_linear(19)) ==  2 );

	multi::layout_t<3> lyt{exts};
	BOOST_TEST_REQUIRE( lyt[get<0>(exts.from_linear(19))][get<1>(exts.from_linear(19))][get<2>(exts.from_linear(19))] == 19 );
	BOOST_TEST_REQUIRE( lyt(get<0>(exts.from_linear(19)), get<1>(exts.from_linear(19)), get<2>(exts.from_linear(19))) == 19 );
}

BOOST_AUTO_TEST_CASE(extensions_from_linear_3d_struct_bind) {
	multi::extensions_t<3> exts{11, 13, 17};

	using std::get;
	{
	auto [eye, jay, kay] = exts.from_linear( 0);
	BOOST_TEST_REQUIRE( eye ==  0 );
	BOOST_TEST_REQUIRE( jay ==  0 );
	BOOST_TEST_REQUIRE( kay ==  0 );
	}
	{
	auto [eye, jay, kay] = exts.from_linear( 1);
	BOOST_TEST_REQUIRE( eye ==  0 );
	BOOST_TEST_REQUIRE( jay ==  0 );
	BOOST_TEST_REQUIRE( kay ==  1 );
	}
	{
	auto [eye, jay, kay] = exts.from_linear(16);
	BOOST_TEST_REQUIRE( eye ==  0 );
	BOOST_TEST_REQUIRE( jay ==  0 );
	BOOST_TEST_REQUIRE( kay == 16 );
	}
	{
	auto [eye, jay, kay] = exts.from_linear(17);
	BOOST_TEST_REQUIRE( eye ==  0 );
	BOOST_TEST_REQUIRE( jay ==  1 );
	BOOST_TEST_REQUIRE( kay ==  0 );
	}
	{
	auto [eye, jay, kay] = exts.from_linear(18);
	BOOST_TEST_REQUIRE( eye ==  0 );
	BOOST_TEST_REQUIRE( jay ==  1 );
	BOOST_TEST_REQUIRE( kay ==  1 );

	multi::layout_t<3> lyt{exts};
	BOOST_TEST_REQUIRE( lyt[eye][jay][kay] == 18 );
	BOOST_TEST_REQUIRE( lyt(eye, jay, kay) == 18 );
	}
}

BOOST_AUTO_TEST_CASE(extensions_from_linear_3d) {
	multi::extensions_t<3> exts{11, 13, 17};

	auto ijk = exts.from_linear(19);

	{
	using std::get;
	BOOST_TEST_REQUIRE( get<0>(exts.from_linear(19)) == 0 );
	BOOST_TEST_REQUIRE( get<1>(exts.from_linear(19)) == 1 );
	BOOST_TEST_REQUIRE( get<2>(exts.from_linear(19)) == 2 );
	}
	{
	using std::get;
//  using multi::detail::get;
	BOOST_TEST_REQUIRE( get<0>(ijk) == 0 );
	BOOST_TEST_REQUIRE( get<1>(ijk) == 1 );
	BOOST_TEST_REQUIRE( get<2>(ijk) == 2 );

	multi::layout_t<3> lyt{exts};

	BOOST_TEST_REQUIRE( lyt[get<0>(ijk)][get<1>(ijk)][get<2>(ijk)] == 19 );
	BOOST_TEST_REQUIRE( lyt(get<0>(ijk), get<1>(ijk), get<2>(ijk)) == 19 );
	}
}

