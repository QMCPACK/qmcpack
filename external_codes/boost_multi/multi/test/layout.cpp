// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// © Alfredo A. Correa 2018-2021

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi layout"
#define BOOST_TEST_DYN_LINK
#include<boost/test/unit_test.hpp>

#include "../array.hpp"
#include "../utility.hpp"

#include "../detail/tuple_zip.hpp"

#include<tuple>

namespace multi = boost::multi;

//BOOST_AUTO_TEST_CASE(tuple_to_extensions){
//	std::tuple<int, int> t{3, 5};
//	auto x = std::apply([](auto... e){return multi::extensions_t<2>{e...};}, t);
//	BOOST_REQUIRE( x.num_elements() == 15 );
//}

BOOST_AUTO_TEST_CASE(linearize) {
	multi::array<double, 3> A({10, 20, 30});
	BOOST_REQUIRE(  25 % extensions(A) == std::make_tuple(0, 0, 25) );
	BOOST_REQUIRE(  55 % extensions(A) == std::make_tuple(0, 1, 25) );
	BOOST_REQUIRE( 655 % extensions(A) == std::make_tuple(1, 1, 25) );
	BOOST_REQUIRE(1255 % extensions(A) == std::make_tuple(2, 1, 25) );

	std::tuple<multi::index, multi::index, multi::index> p = A.extensions().from_linear(655);
	BOOST_REQUIRE( p == std::make_tuple(1, 1, 25) );
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
#if defined(__INTEL_COMPILER) or (defined(__GNUC__) and (__GNUC__ < 6))
		multi::extensions_t<2>
#endif
		{50, 50}
	);
	BOOST_REQUIRE( size(B2) == 50 ); BOOST_REQUIRE( B2.size() == 50 );
	BOOST_REQUIRE( B2[0].sliced(10, 20).size() == 10 );
	BOOST_REQUIRE( size(B2[0].sliced(10, 20))  == 10 );

	static_assert( decltype(B2(0, {10, 20}))::rank_v  == 1 , "!");

	BOOST_REQUIRE( size(B2(0, {10, 20})) == 10 );
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

	std::array<std::array<decltype(B2({0, 2}, {0, 2})), 2>, 2> B2blk = {
		{
			{B2({0, 2}, {0, 2}), B2({0, 2}, {2, 4})},
			{B2({2, 4}, {0, 2}), B2({2, 4}, {2, 4})}
		}
	};

	BOOST_REQUIRE( &B2blk[1][1][1][1] == &B2[3][3] );
}
{
	// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays, hicpp-avoid-c-arrays,-warnings-as-errors, modernize-avoid-c-arrays,-warnings-as-errors): test
	double A[3][4][5] = {};
	using multi::dimensionality;
	static_assert(dimensionality(A)==3, "!");
	using multi::extensions;
	auto xA = extensions(A);
	using std::get;
	BOOST_REQUIRE( size(get<0>(xA)) == 3 );
	BOOST_REQUIRE( size(get<1>(xA)) == 4 );
	BOOST_REQUIRE( size(get<2>(xA)) == 5 );

	multi::array<double, 3> AA({3, 4, 5});
	using multi::layout;
	BOOST_REQUIRE( layout(AA) == layout(A) );

	using multi::stride;
	BOOST_REQUIRE( stride(AA) == 20 );
	static_assert( stride(A) == 20 , "!" );
	static_assert( stride(A[0]) == 5 , "!" );
	static_assert( stride(A[1]) == 5 , "!" );
	static_assert( stride(A[0][0]) == 1 , "!" );
//		assert( stride(A) == 20 );
//		assert( stride(A[0]) == 20 );
}
{
	std::array<std::array<std::array<double, 5>, 4>, 3> A = {};
	using multi::dimensionality;
	static_assert(dimensionality(A)==3, "!");
	using multi::extensions;
	auto xA = extensions(A);
	using std::get;
	BOOST_REQUIRE( size(get<0>(xA)) == 3 );
	BOOST_REQUIRE( size(get<1>(xA)) == 4 );
	BOOST_REQUIRE( size(get<2>(xA)) == 5 );

	multi::array<double, 3> AA({3, 4, 5});
	using multi::layout;
	BOOST_REQUIRE( layout(AA) == layout(A) );

	using multi::stride;
	BOOST_REQUIRE( stride(AA) == 20 );
	static_assert( stride(A) == 20 , "!" );
	BOOST_REQUIRE( stride(A[0]) == 5 );
	BOOST_REQUIRE( stride(A[1]) == 5 );
	BOOST_REQUIRE( stride(A[0][0]) == 1 );
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
	BOOST_REQUIRE( stride(B2) == 1 );
	BOOST_REQUIRE( stride(B2[0]) == 1 );
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
	BOOST_REQUIRE( L[0][0][0] == &A[0][0][0] - data_elements(A) );
	BOOST_REQUIRE( L[0][0][1] == &A[0][0][1] - data_elements(A) );
	BOOST_REQUIRE( L[0][0][2] == &A[0][0][2] - data_elements(A) );
	BOOST_REQUIRE( L[0][1][2] == &A[0][1][2] - data_elements(A) );
	BOOST_REQUIRE( L[3][1][2] == &A[3][1][2] - data_elements(A) );
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

	BOOST_REQUIRE( stride(L) == 20*30 );
	BOOST_REQUIRE( offset(L) == 0 );
	BOOST_REQUIRE( nelems(L) == 10*20*30 );

	BOOST_REQUIRE( L.stride(0) == stride(L) );
	BOOST_REQUIRE( L.offset(0) == offset(L) );
	BOOST_REQUIRE( L.nelems(0) == nelems(L) );

	BOOST_REQUIRE( L.stride(1) == 30 );
	BOOST_REQUIRE( L.offset(1) == 0 );
	BOOST_REQUIRE( L.nelems(1) == 20*30 );

	BOOST_REQUIRE( L.stride(2) == 1 );
	BOOST_REQUIRE( L.offset(2) == 0 );
	BOOST_REQUIRE( L.nelems(2) == 30 );
}

BOOST_AUTO_TEST_CASE(continued_part3) {
	multi::layout_t<3> L({{0, 10}, {0, 20}, {0, 30}});

	BOOST_REQUIRE( L.num_elements() == num_elements(L) );
	BOOST_REQUIRE( L.size() == size(L) );
	BOOST_REQUIRE( L.extension() == extension(L) );

	BOOST_REQUIRE( num_elements(L) == 10*20*30 );
	BOOST_REQUIRE( size(L) == 10 );
	BOOST_REQUIRE( extension(L).first() == 0 );
	BOOST_REQUIRE( extension(L).last() == 10 );

	BOOST_REQUIRE( L.size(1) == 20 );
	BOOST_REQUIRE( L.extension(1).first() == 0 );
	BOOST_REQUIRE( L.extension(1).last() == 20 );

	BOOST_REQUIRE( L.size(2) == 30 );
	BOOST_REQUIRE( L.extension(2).first() == 0 );
	BOOST_REQUIRE( L.extension(2).last() == 30 );

	using std::get;
	BOOST_REQUIRE( get<0>(strides(L)) == L.stride(0) );
	BOOST_REQUIRE( get<1>(strides(L)) == L.stride(1) );
	BOOST_REQUIRE( get<2>(strides(L)) == L.stride(2) );

	auto const& strides = L.strides();
	BOOST_REQUIRE( get<0>(strides) == L.stride(0) );
}

BOOST_AUTO_TEST_CASE(continued) {
{
	multi::layout_t<3> L;
	BOOST_REQUIRE( size(L) == 0 );
}{
	multi::layout_t<3> L( {{0, 10}, {0, 20}, {0, 30}} );
	BOOST_REQUIRE( stride(L) == 20*30 );
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
	BOOST_REQUIRE( stride(L) == 5*7 );
	BOOST_REQUIRE( stride(L.sub().sub()) == 1 );
}
{
	multi::layout_t<3> L({{0, 10}, {0, 20}, {0, 30}});
	BOOST_REQUIRE( stride(L) == 20*30 );
	BOOST_REQUIRE( offset(L) == 0 );
	BOOST_REQUIRE( nelems(L) == 10*20*30 );
}
{
	multi::layout_t<3> L({{10, 20}, {10, 30}, {10, 40}});
	BOOST_REQUIRE( stride(L) == 20*30 );
}
{
	std::tuple<int, int, int> ttt = {1, 2, 3};
	auto arrr = boost::multi::detail::to_array(ttt);
	BOOST_REQUIRE(arrr[1] == 2);
}
}

BOOST_AUTO_TEST_CASE(tuple_zip_test) {
	auto t1 = std::make_tuple( 1,  2,  3);
	auto t2 = std::make_tuple(10, 20, 30);
	auto t3 = std::make_tuple(std::string{"10"}, std::string{"20"}, std::string{"30"});
	auto t123 = boost::multi::detail::tuple_zip(t1, t2, t3);
	BOOST_REQUIRE( std::get<2>(std::get<0>(t123)) == std::string{"10"} );
}

