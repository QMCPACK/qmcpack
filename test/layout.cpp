// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2021 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi layout"
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"
#include "multi/utility.hpp"

#include "multi/detail/tuple_zip.hpp"

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

#include<tuple>

namespace multi = boost::multi;

static auto second_finish(multi::extensions_t<3> x) {
	return std::get<1>(x).finish();
}

BOOST_AUTO_TEST_CASE(extensions_3D) {
	BOOST_REQUIRE( 20 == second_finish( multi::extensions_t<3>  { {0, 10}, {0, 20}, {0, 30} }  ) );
	BOOST_REQUIRE( 20 == second_finish( multi::extensions_t<3>( { {0, 10}, {0, 20}, {0, 30} } )) );
	BOOST_REQUIRE( 20 == second_finish(                         { {0, 10}, {0, 20}, {0, 30} }  ) );

	multi::extensions_t<3> x3{ {0, 10}, {0, 20}, {0, 30} };
	BOOST_REQUIRE( 20 == second_finish(x3                                                     ) );
}

//BOOST_AUTO_TEST_CASE(extensions_1D) {
//	BOOST_REQUIRE( multi::extensions_t<1>( { {0, 10} } ) == multi::extensions_t<1>( { {0, 10} } ) );
//}

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

BOOST_AUTO_TEST_CASE(extensions_to_linear) {
	multi::extensions_t<3> x{4, 5, 3};
	BOOST_REQUIRE( x.to_linear(0, 0, 0) ==  0 );
	BOOST_REQUIRE( x.to_linear(0, 0, 1) ==  1 );
	BOOST_REQUIRE( x.to_linear(0, 0, 2) ==  2 );
	BOOST_REQUIRE( x.to_linear(0, 1, 0) ==  3 );
	BOOST_REQUIRE( x.to_linear(0, 1, 1) ==  4 );
	BOOST_REQUIRE( x.to_linear(0, 1, 2) ==  5 );
	BOOST_REQUIRE( x.to_linear(1, 0, 0) == 15 );

	for(int i = 0; i != 4; ++i) {
		for(int j = 0; j != 5; ++j) {
			for(int k = 0; k != 3; ++k) {
				BOOST_REQUIRE(( x.from_linear(x.to_linear(i, j, k)) ==  decltype(x.from_linear(x.to_linear(i, j, k))){i, j, k} ));
			}
		}
	}

	BOOST_REQUIRE( x.to_linear(4, 0, 0) == x.num_elements() );

	for(int n = 0; n != x.num_elements(); ++n) {
		BOOST_REQUIRE( std::apply([&](auto... e){return x.to_linear(e...);}, x.from_linear(n)) == n );
	}
}

BOOST_AUTO_TEST_CASE(extensions_layout_to_linear) {
	multi::array<double, 3> A({40, 50, 80});
	auto&& B = A({10, 30}, {20, 32}, {60, 75});

	for(int i = 0; i != 10; ++i) {
		for(int j = 0; j != 12; ++j) {
			for(int k = 0; k != 15; ++k) {
				BOOST_REQUIRE( &  B.base()  [B.layout()(i, j, k)] == &B(i, j, k) );
				BOOST_REQUIRE( &*(B.base() + B.layout()(i, j, k)) == &B(i, j, k) );
			}
		}
	}
}

BOOST_AUTO_TEST_CASE(extensions_layout_to_linear_2) {
	multi::array<double, 3> A({40, 50, 80});
	auto&& B = A({10, 30}, {20, 32}, {60, 75});

	auto const& C = B.rotated();
	auto Cx = C.extensions();
	for(auto i : std::get<0>(Cx)) {
		for(auto j : std::get<1>(Cx)) {
			for(auto k : std::get<2>(Cx)) {
				BOOST_REQUIRE( &  C.base()  [C.layout()(i, j, k)] == &C(i, j, k) );
				BOOST_REQUIRE( &*(C.base() + C.layout()(i, j, k)) == &C(i, j, k) );
			}
		}
	}
}

BOOST_AUTO_TEST_CASE(linearize) {
	multi::array<double, 3> A({10, 20, 30});

	BOOST_REQUIRE((  25 % extensions(A) == decltype(  25 % extensions(A)){0, 0, 25} ));
	BOOST_REQUIRE((  55 % extensions(A) == decltype(  55 % extensions(A))(0, 1, 25) ));
	BOOST_REQUIRE(( 655 % extensions(A) == decltype( 655 % extensions(A))(1, 1, 25) ));
	BOOST_REQUIRE((1255 % extensions(A) == decltype(1255 % extensions(A))(2, 1, 25) ));

	auto const p = A.extensions().from_linear(655);
//  BOOST_REQUIRE( p == std::make_tuple(1, 1, 25) );
	using multi::detail::get;
	BOOST_REQUIRE( get<0>(p) ==  1 );
	BOOST_REQUIRE( get<1>(p) ==  1 );
	BOOST_REQUIRE( get<2>(p) == 25 );
}

BOOST_AUTO_TEST_CASE(layout_0) {
	multi::array<double, 3> A3(
#if defined(__INTEL_COMPILER) or (defined(__GNUC__) and (__GNUC__ < 6))
		multi::extensions_t<3>
#endif
		{51, 52, 53}
	);
	BOOST_REQUIRE( size(A3) == 51      ); BOOST_REQUIRE( A3.size() == 51       );
	BOOST_REQUIRE( size(A3[0]) == 52   ); BOOST_REQUIRE( A3[0].size() == 52    );
	BOOST_REQUIRE( size(A3[0][0]) == 53); BOOST_REQUIRE( A3[0][0].size() == 53 );
}

BOOST_AUTO_TEST_CASE(layout_1) {
	//NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays): testing feature
	double DA[50][50][50];
	using multi::size;
	BOOST_REQUIRE( size(DA) == 50 );

	using multi::extension;
	BOOST_REQUIRE(( extension(DA) == multi::index_extension{0, 50} ));
	BOOST_REQUIRE(( extension(DA) == multi::iextension{0, 50}      ));
	BOOST_REQUIRE(( extension(DA) == multi::irange{0, 50} ));
}

BOOST_AUTO_TEST_CASE(layout_2) {
	std::array<std::array<std::array<double, 50>, 50>, 50> DA{};
	using multi::size;
	BOOST_REQUIRE( size(DA) == 50 );

	using multi::extension;
	BOOST_REQUIRE(( extension(DA) == multi::index_extension{0, 50} ));
	BOOST_REQUIRE(( extension(DA) == multi::iextension{0, 50}      ));
	BOOST_REQUIRE(( extension(DA) == multi::irange{0, 50}          ));
}

BOOST_AUTO_TEST_CASE(layout_3) {
	multi::array<double, 2> B2(
//#if defined(__INTEL_COMPILER) or (defined(__GNUC__) and (__GNUC__ < 6))
//		multi::extensions_t<2>
//#endif
		{50, 50}
	);
	BOOST_REQUIRE( size(B2) == 50 ); BOOST_REQUIRE( B2.size() == 50 );
	BOOST_REQUIRE( B2[0].sliced(10, 20).size() == 10 );
	BOOST_REQUIRE( size(B2[0].sliced(10, 20))  == 10 );

	static_assert( decltype(B2(0, {10, 20}))::rank_v  == 1 , "!");

	BOOST_REQUIRE( size(B2(0, {10, 20})) == 10 );

	BOOST_REQUIRE(      B2.layout() == B2.layout()  );
	BOOST_REQUIRE( not (B2.layout() <  B2.layout()) );
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
	double A[3][4][5] = {};
	using multi::dimensionality;
	static_assert(dimensionality(A)==3, "!");
	using multi::extensions;
	auto xA = extensions(A);

	BOOST_REQUIRE( size(std::get<0>(xA)) == 3 );
	BOOST_REQUIRE( size(std::get<1>(xA)) == 4 );
	BOOST_REQUIRE( size(std::get<2>(xA)) == 5 );

	static_assert( multi::stride(A)       == 20 , "!" );

//  static_assert( multi::stride(A)       ==  5 , "!" );
	static_assert( multi::stride(A[1])    ==  5 , "!" );
	static_assert( multi::stride(A[0][0]) ==  1 , "!" );

	multi::array<double, 3> AA({3, 4, 5});
	using multi::layout;
	BOOST_REQUIRE( layout(AA) == layout(A) );

	BOOST_REQUIRE( AA     .stride() == 20 );
}
{
	std::array<std::array<std::array<double, 5>, 4>, 3> A = {};
#if defined(__circle_build__)  // circle doesn't see dimensionality as a constexpr "cannot access value of A at compile time;"
	       assert( multi::dimensionality(A) == 3 );
#else  // other compilers ok
	static_assert( multi::dimensionality(A) == 3 );
#endif

	using multi::extensions;
	auto xA = extensions(A);
	using std::get;
	BOOST_REQUIRE( size(std::get<0>(xA)) == 3 );
	BOOST_REQUIRE( size(std::get<1>(xA)) == 4 );
	BOOST_REQUIRE( size(std::get<2>(xA)) == 5 );

	multi::array<double, 3> AA({3, 4, 5});
	using multi::layout;
	BOOST_REQUIRE( layout(AA) == layout(A) );

	BOOST_REQUIRE( AA.stride() == 20 );

#if defined(__circle_build__)  // circle doesn't recognize this as a constexpr "cannot access value of A at compile time;"
	       assert( multi::stride(A) == 20);
#else  // other compilers ok
	static_assert( multi::stride(A) == 20);
#endif

	BOOST_REQUIRE( multi::stride(A[0])    == 5 );
	BOOST_REQUIRE( multi::stride(A[1])    == 5 );
	BOOST_REQUIRE( multi::stride(A[0][0]) == 1 );
//		assert( stride(A) == 20 );
//		assert( stride(A[0]) == 20 );
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
		boost::multi::layout_t<2>::extensions_type x{
			multi::iextension(2, 5),
			multi::iextension(0, 5)
		};
		multi::layout_t<2> l2(x);
		BOOST_REQUIRE( l2.extension().start()  == std::get<0>(x).start()  );
		BOOST_REQUIRE( l2.extension().finish() == std::get<0>(x).finish() );
	}
	 {
		multi::layout_t<2> l2({multi::iextension(0, 3), multi::iextension(2, 7)});
		BOOST_REQUIRE( std::get<1>(l2.extensions()).start()  == 2 );
		BOOST_REQUIRE( std::get<1>(l2.extensions()).finish() == 7 );
	}
}

BOOST_AUTO_TEST_CASE(multi_layout_part1) {
{
	multi::layout_t<0> L;
	static_assert( decltype(L)::rank_v==0 , "!");
	BOOST_REQUIRE( num_elements(L) == 1 );
}{
	multi::iextensions<0> x{};
	multi::layout_t<0> L(x);
	BOOST_REQUIRE(L.num_elements() == 1);
}{  multi::layout_t<1> L{};
	static_assert( decltype(L)::rank_v == 1 , "!");
	BOOST_REQUIRE( num_elements(L) == 0 );
	BOOST_REQUIRE( size(L) == 0 );
	BOOST_REQUIRE( size(extension(L))==0 );
	BOOST_REQUIRE( stride(L)!=0 );
	BOOST_REQUIRE( is_empty(L) );
}{
	multi::layout_t<2> L({2, 10});
	static_assert( decltype(L)::rank_v == 2 , "!");
	BOOST_REQUIRE( num_elements(L) == 20 );
	BOOST_REQUIRE( size(L) == 2 );
	BOOST_REQUIRE( size(extension(L))==2 );
	BOOST_REQUIRE( stride(L)==10 );
	BOOST_REQUIRE( not is_empty(L) );
} {
	multi::layout_t<1> L(multi::iextensions<1>{20});
	static_assert( decltype(L)::rank_v == 1 , "!");
	BOOST_REQUIRE( num_elements(L) == 20 );
	BOOST_REQUIRE( size(L) == 20 );
	BOOST_REQUIRE( stride(L) == 1 );
}
}

BOOST_AUTO_TEST_CASE(multi_layout_part2) {
{
	multi::layout_t<1> L(multi::iextensions<1>{1});
	static_assert( decltype(L)::rank_v ==1 , "!");
	BOOST_REQUIRE( num_elements(L) == 1 );
	BOOST_REQUIRE( size(L) == 1 );
	BOOST_REQUIRE( stride(L) == 1 );
} {
	multi::layout_t<2> L({1, 10});
	static_assert( decltype(L)::rank_v ==2 , "!");
	BOOST_REQUIRE( num_elements(L) == 10 );
	BOOST_REQUIRE( size(L) == 1);
	BOOST_REQUIRE( not is_empty(L) );
	BOOST_REQUIRE( size(extension(L))==1 );
	BOOST_REQUIRE( stride(L)== 10 );//std::numeric_limits<std::ptrdiff_t>::max() );
	using std::get;
	BOOST_REQUIRE( get<0>(strides(L)) == 10);
	BOOST_REQUIRE( get<1>(strides(L)) == 1 );
}
}

BOOST_AUTO_TEST_CASE(multi_layout_part3) {
{
	multi::layout_t<2> L({10, 1});
	static_assert( decltype(L)::rank_v ==2 , "!");
	BOOST_REQUIRE( num_elements(L) == 10 );
	BOOST_REQUIRE( size(L) == 10 );
	using std::get;
	BOOST_REQUIRE( get<0>(strides(L)) == 1 );
	BOOST_REQUIRE( get<1>(strides(L)) == 1 );
}{  multi::layout_t<2> L{};
	BOOST_REQUIRE( dimensionality(L)==2 );
	BOOST_REQUIRE( num_elements(L) == 0 );
	BOOST_REQUIRE( size(L) == 0 );
	BOOST_REQUIRE( size(extension(L))==0 );
	BOOST_REQUIRE( stride(L)!=0 );
	BOOST_REQUIRE( is_empty(L) );
}{  multi::layout_t<3> L{}; BOOST_REQUIRE( num_elements(L) == 0 );
}{	multi::layout_t<3> L({{0, 10}, {0, 10}, {0, 10}}); BOOST_REQUIRE( num_elements(L) == 1000 );
}{	multi::layout_t<3> L({{10}, {10}, {10}}); BOOST_REQUIRE( num_elements(L) == 1000 );
}{	multi::layout_t<3> L({10, 10, 10}); BOOST_REQUIRE( num_elements(L) == 1000 );
}{	multi::layout_t<3> L({multi::index_extension{0, 10}, {0, 10}, {0, 10}}); BOOST_REQUIRE( num_elements(L) == 1000 );
}{	multi::layout_t<3> L(multi::layout_t<3>::extensions_type{{0, 10}, {0, 10}, {0, 10}}); BOOST_REQUIRE( num_elements(L) == 1000 );
}
}

BOOST_AUTO_TEST_CASE(layout_to_offset) {
	multi::layout_t<3> L({10, 20, 30});
	multi::array<double, 3> A({10, 20, 30});
	BOOST_REQUIRE( L[0][0][0] == &A[0][0][0] - A.data_elements() );
	BOOST_REQUIRE( L[0][0][1] == &A[0][0][1] - A.data_elements() );
	BOOST_REQUIRE( L[0][0][2] == &A[0][0][2] - A.data_elements() );

	BOOST_TEST_REQUIRE( L[0][1][2] == &A[0][1][2] - A.data_elements() );
	BOOST_TEST_REQUIRE( L[3][1][2] == &A[3][1][2] - A.data_elements() );
}

BOOST_AUTO_TEST_CASE(layout_to_offset_sub) {
	multi::array<double, 3> A({10, 20, 30});
	auto&& s = A({2, 6}, {4, 8}, {10, 20});
	auto l = s.layout();
	BOOST_REQUIRE( l[0][0][0] == &s[0][0][0] - base(s) );
	BOOST_REQUIRE( l[0][0][1] == &s[0][0][1] - base(s) );
	BOOST_REQUIRE( l[0][0][2] == &s[0][0][2] - base(s) );
	BOOST_REQUIRE( l[0][1][2] == &s[0][1][2] - base(s) );
	BOOST_REQUIRE( l[3][1][2] == &s[3][1][2] - base(s) );
}

BOOST_AUTO_TEST_CASE(continued_part1) {
{
	multi::layout_t<3> L(multi::layout_t<3>::extensions_type{{0, 10}, {0, 10}, {0, 10}});
	BOOST_REQUIRE( num_elements(L) == 1000);
}
{	multi::layout_t<3> L({multi::iextension{0, 10}, multi::iextension{0, 10}, multi::iextension{0, 10}}); BOOST_REQUIRE(L.num_elements() == 1000);
}{	multi::layout_t<3> L({multi::iextension{10}, multi::iextension{10}, multi::iextension{10}}); BOOST_REQUIRE( num_elements(L) == 1000);
}{	multi::layout_t<3> L({10, 10, multi::iextension{10}}); BOOST_REQUIRE( num_elements(L) == 1000 );
}{
	multi::layout_t<1> L;
	BOOST_REQUIRE( size(L) == 0 );
}{
	multi::layout_t<1> L({{0, 10}});
	BOOST_REQUIRE( size(L) == 10 );
	BOOST_REQUIRE( extension(L).start () ==  0 );
	BOOST_REQUIRE( extension(L).finish() == 10 );

	L.reindex(1);
	BOOST_REQUIRE( size(L) == 10 );
	BOOST_REQUIRE( extension(L).start () ==  1 );
	BOOST_REQUIRE( extension(L).finish() == 11 );
}{
	multi::layout_t<2> L;
	BOOST_REQUIRE( size(L) == 0 );
}
{
	multi::layout_t<2> L({{0, 10}, {0, 20}});
	BOOST_REQUIRE( size(L) == 10 );
	BOOST_REQUIRE( extension(L).start () ==  0 );
	BOOST_REQUIRE( extension(L).finish() == 10 );

	L.reindex(1);
	BOOST_REQUIRE( extension(L).start () ==  1 );
	BOOST_REQUIRE( extension(L).finish() == 11 );

	L.rotate().reindex(3).unrotate();
	BOOST_TEST_REQUIRE( extension(L).start () ==  1 );
	BOOST_TEST_REQUIRE( extension(L).finish() == 11 );

	BOOST_TEST_REQUIRE( std::get<0>(extensions(L)).start () == 1 );
	BOOST_TEST_REQUIRE( std::get<1>(extensions(L)).start () == 3 );
	BOOST_TEST_REQUIRE( std::get<1>(extensions(L)).finish() == 23 );
}
}

BOOST_AUTO_TEST_CASE(continued_part2) {
	multi::layout_t<3> L({{0, 10}, {0, 20}, {0, 30}});

	BOOST_REQUIRE( not L.empty() );

	BOOST_REQUIRE( stride(L) == L.stride() );
	BOOST_REQUIRE( offset(L) == L.offset() );
	BOOST_REQUIRE( nelems(L) == L.nelems() );

	BOOST_REQUIRE( stride(L) == 20*30L );
	BOOST_REQUIRE( offset(L) == 0 );
	BOOST_REQUIRE( nelems(L) == 10*20L*30L );

	BOOST_REQUIRE( L.stride() == stride(L) );
	BOOST_REQUIRE( L.offset() == offset(L) );
	BOOST_REQUIRE( L.nelems() == nelems(L) );

	using boost::multi::detail::get;
	BOOST_REQUIRE( get<1>(L.strides()) == 30     );
	BOOST_REQUIRE( get<1>(L.offsets()) ==  0     );
	BOOST_REQUIRE( get<1>(L.nelemss()) == 20*30L );

	BOOST_REQUIRE( get<2>(L.strides()) ==  1 );
	BOOST_REQUIRE( get<2>(L.offsets()) ==  0 );
	BOOST_REQUIRE( get<2>(L.nelemss()) == 30 );
}

BOOST_AUTO_TEST_CASE(continued_part3) {
	multi::layout_t<3> L({{0, 10}, {0, 20}, {0, 30}});

	BOOST_REQUIRE( L.num_elements() == num_elements(L) );
	BOOST_REQUIRE( L.size() == size(L) );
	BOOST_REQUIRE( L.extension() == extension(L) );

	BOOST_REQUIRE( num_elements(L) == 10*20L*30L );
	BOOST_REQUIRE( size(L) == 10 );
	BOOST_REQUIRE( extension(L).first() == 0 );
	BOOST_REQUIRE( extension(L).last() == 10 );

	BOOST_REQUIRE( std::get<0>(L.extensions()) == L.extension() );


	boost::multi::extensions_t<2> x2;

	using boost::multi::detail::get;
	using std::get;

	BOOST_REQUIRE( get<0>(x2).is_empty() );

//	BOOST_REQUIRE( std::get<0>(L.sizes()) == L.size(0) );
//	BOOST_REQUIRE( std::get<0>(L.extensions()) == L.extension(0) );

	BOOST_REQUIRE(( get<0>(L.extensions()) == multi::index_extension{0, 10} ));

	BOOST_REQUIRE( get<0>(L.extensions()).first() ==  0 );
	BOOST_REQUIRE( get<0>(L.extensions()).last()  == 10 );

//  BOOST_REQUIRE( L.size(1) == 20 );
	BOOST_REQUIRE( get<1>(L.extensions()).first() ==  0 );
	BOOST_REQUIRE( get<1>(L.extensions()).last()  == 20 );

//  BOOST_REQUIRE( L.size(2) == 30 );
	BOOST_REQUIRE( get<2>(L.extensions()).first() ==  0 );
	BOOST_REQUIRE( get<2>(L.extensions()).last()  == 30 );

	using std::get;
	BOOST_REQUIRE( get<0>(strides(L)) == L.stride() );

	auto const& strides = L.strides();
	BOOST_REQUIRE( get<0>(strides) == L.stride() );
}

BOOST_AUTO_TEST_CASE(continued) {
{
	multi::layout_t<3> L;
	BOOST_REQUIRE( size(L) == 0 );
}
{
	multi::layout_t<3> L( {{0, 10}, {0, 20}, {0, 30}} );
	BOOST_REQUIRE( stride(L) == 20*30L );
}
{
	multi::layout_t<1> L({{0, 10}});
	BOOST_REQUIRE( extension(L).first() == 0 );
	BOOST_REQUIRE( extension(L).last() == 10 );
}
{
	multi::layout_t<1> L({{8, 18}});
	BOOST_REQUIRE( extension(L).first() == 8 );
	BOOST_REQUIRE( extension(L).last() == 18 );
}
{
	multi::layout_t<2> L({{0, 10}, {0, 20}});
	BOOST_REQUIRE( extension(L).first() == 0 );
	BOOST_REQUIRE( extension(L).last() == 10 );
}
{
	multi::layout_t<2> L( {{0, 10}, {11, 31}} );
	BOOST_REQUIRE( size(L) == 10   );
	BOOST_REQUIRE( stride(L) == 20 );
	BOOST_REQUIRE( offset(L) == 0 );
}
{
	multi::layout_t<2> L( {{8, 18}, {0, 20}} );
	BOOST_REQUIRE( size(L) == 10 );
	BOOST_REQUIRE( stride(L) == 20 );
}
{
	multi::layout_t<3> L({{0, 3}, {0, 5}, {10, 17}});
	BOOST_REQUIRE( stride(L) == 5*7L );
	BOOST_REQUIRE( stride(L.sub().sub()) == 1 );
}
{
	multi::layout_t<3> L({{0, 10}, {0, 20}, {0, 30}});
	BOOST_REQUIRE( stride(L) == 20*30L );
	BOOST_REQUIRE( offset(L) == 0 );
	BOOST_REQUIRE( nelems(L) == 10*20L*30L );
}
{
	multi::layout_t<3> L({{10, 20}, {10, 30}, {10, 40}});
	BOOST_REQUIRE( stride(L) == 20*30L );
}
{
	auto const ttt = boost::multi::tuple<int, int, int>{1, 2, 3};
	auto const arr = std::apply([](auto... es) {return std::array<int, 3>{{es...}};}, ttt);
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
	multi::extensions_t<1> x{11};

	auto ijk = x.from_linear(9);

	using multi::detail::get;
	BOOST_TEST_REQUIRE( get<0>(ijk) == 9 );

	multi::layout_t<1> l{x};
	BOOST_TEST_REQUIRE( l[get<0>(ijk)] == 9 );
	BOOST_TEST_REQUIRE( l(get<0>(ijk)) == 9 );

//	BOOST_TEST_REQUIRE( l(std::get<0>(l.extensions().from_linear(9))) == 9 );

#if(__cplusplus >= 201703)
//	BOOST_TEST_REQUIRE( std::apply(l, l.extensions().from_linear(9)) == 9 );
#endif
}

BOOST_AUTO_TEST_CASE(extensions_from_linear_2d) {
	multi::extensions_t<2> x{3, 5};

	auto ij = x.from_linear(7);

	using multi::detail::get;

	BOOST_TEST_REQUIRE( get<0>(ij) == 1 );
	BOOST_TEST_REQUIRE( get<1>(ij) == 2 );

	multi::layout_t<2> l{x};
	BOOST_TEST_REQUIRE( l[get<0>(ij)][get<1>(ij)] == 7 );
//	BOOST_TEST_REQUIRE( l(std::get<0>(ij), std::get<1>(ij)) == l[std::get<0>(ij)](std::get<1>(ij)) );
//	BOOST_TEST_REQUIRE( l[std::get<0>(ij)](std::get<1>(ij)) == l[std::get<0>(ij)][std::get<1>(ij)] );

//	BOOST_TEST_REQUIRE( l(std::get<0>(ij), std::get<1>(ij)) == 7 );

//	BOOST_TEST_REQUIRE( l(std::get<0>(l.extensions().from_linear(7)), std::get<1>(l.extensions().from_linear(7))) == 7 );

#if(__cplusplus >= 201703)
	auto [i, j] = x.from_linear(7);

	BOOST_TEST_REQUIRE( i == 1 );
	BOOST_TEST_REQUIRE( j == 2 );
//	BOOST_TEST_REQUIRE( std::apply(l, l.extensions().from_linear(9)) == 9 );
#endif
}

BOOST_AUTO_TEST_CASE(extensions_from_linear_3d) {
	multi::extensions_t<3> x{11, 13, 17};

	auto ijk = x.from_linear(19);

	using multi::detail::get;
	BOOST_TEST_REQUIRE( get<0>(ijk) == 0 );
	BOOST_TEST_REQUIRE( get<1>(ijk) == 1 );
	BOOST_TEST_REQUIRE( get<2>(ijk) == 2 );

	multi::layout_t<3> l{x};

	BOOST_TEST_REQUIRE( l[get<0>(ijk)][get<1>(ijk)][get<2>(ijk)] == 19 );
//	BOOST_TEST_REQUIRE( l(std::get<0>(ijk), std::get<1>(ijk), std::get<2>(ijk)) == 19 );
}

