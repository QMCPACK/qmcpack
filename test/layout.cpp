// Copyright 2018-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for range, layout_t, get, extensions_t

#include <array>     // for array, array<>::value_type
#include <iterator>  // for size
#if __cplusplus > 201703L
#if __has_include(<ranges>)
#include <ranges>  // NOLINT(misc-include-cleaner) IWYU pragma: keep
#endif
#endif
#include <tuple>     // for make_tuple, tuple_element<>::type
// IWYU pragma: no_include <type_traits>
// IWYU pragma: no_#include <version>

namespace multi = boost::multi;

namespace {
auto second_finish(multi::extensions_t<3> exts) {
	using std::get;  // workaround: function call with explicit template arguments is a C++20 extension [-Wc++20-extensions]
	return get<1>(exts).last();
}
}  // namespace

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
BOOST_AUTO_TEST_CASE(extensions_3D) {
	BOOST_TEST( 20 == second_finish( multi::extensions_t<3>  { {0, 10}, {0, 20}, {0, 30} }  ) );
	BOOST_TEST( 20 == second_finish( multi::extensions_t<3>( { {0, 10}, {0, 20}, {0, 30} } )) );
	BOOST_TEST( 20 == second_finish(                         { {0, 10}, {0, 20}, {0, 30} }  ) );

	multi::extensions_t<3> const exts({0, 10}, {0, 20}, {0, 30});
	BOOST_TEST( 20 == second_finish(exts) );
}


BOOST_AUTO_TEST_CASE(extensions_to_linear) {
	multi::extensions_t<3> exts{4, 5, 3};
	BOOST_TEST( exts.to_linear(0, 0, 0) ==  0 );
	BOOST_TEST( exts.to_linear(0, 0, 1) ==  1 );
	BOOST_TEST( exts.to_linear(0, 0, 2) ==  2 );
	BOOST_TEST( exts.to_linear(0, 1, 0) ==  3 );
	BOOST_TEST( exts.to_linear(0, 1, 1) ==  4 );
	BOOST_TEST( exts.to_linear(0, 1, 2) ==  5 );
	BOOST_TEST( exts.to_linear(1, 0, 0) == 15 );

	for(int eye = 0; eye != 4; ++eye) {
		for(int jay = 0; jay != 5; ++jay) {
			for(int kay = 0; kay != 3; ++kay) {  // NOLINT(altera-unroll-loops)
				BOOST_TEST(( exts.from_linear(exts.to_linear(eye, jay, kay)) == decltype(exts.from_linear(exts.to_linear(eye, jay, kay))){eye, jay, kay} ));
			}
		}
	}

	BOOST_TEST( exts.to_linear(4, 0, 0) == exts.num_elements() );

	for(int idx = 0; idx != exts.num_elements(); ++idx) {  // NOLINT(altera-unroll-loops)
		BOOST_TEST( std::apply([&](auto... indices) { return exts.to_linear(indices...);}, exts.from_linear(idx)) == idx );
	}
}

BOOST_AUTO_TEST_CASE(extensions_layout_to_linear) {
	multi::array<double, 3> arr(
	#ifdef _MSC_VER  // problem with MSVC 14.3 c++17
		multi::extensions_t<3>
	#endif
		{40, 50, 80}
	);

	auto&& sub = arr({10, 30}, {20, 32}, {60, 75});

	#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunknown-warning-option"
	#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
	#endif

	for(int i = 0; i != 10; ++i) {
		for(int j = 0; j != 12; ++j) {
			for(int k = 0; k != 15; ++k) {  // NOLINT(altera-unroll-loops)
				BOOST_TEST( &  sub.base()  [sub.layout()(i, j, k)] == &sub(i, j, k) );
				BOOST_TEST( &*(sub.base() + sub.layout()(i, j, k)) == &sub(i, j, k) );
			}
		}
	}

	#if defined(__clang__)
	#pragma clang diagnostic pop
	#endif
}

BOOST_AUTO_TEST_CASE(extensions_layout_to_linear_2) {
	multi::array<double, 3> arr(
	#ifdef _MSC_VER  // problem with MSVC 14.3 c++17
		multi::extensions_t<3>
	#endif
		{40, 50, 80}
	);

	auto&& sub = arr({10, 30}, {20, 32}, {60, 75});

	auto const& rot = sub.rotated();

	auto const [is, js, ks] = rot.extensions();

	#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunknown-warning-option"
	#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
	#endif

	for(auto const i : is) {
		for(auto const j : js) {
			for(auto const k : ks) {  // NOLINT(altera-unroll-loops)
				BOOST_TEST( &  rot.base()  [rot.layout()(i, j, k)] == &rot(i, j, k) );
				BOOST_TEST( &*(rot.base() + rot.layout()(i, j, k)) == &rot(i, j, k) );
			}
		}
	}

	#if defined(__clang__)
	#pragma clang diagnostic pop
	#endif
}

BOOST_AUTO_TEST_CASE(linearize) {
	multi::array<double, 3> const arr(
	#ifdef _MSC_VER  // problem with MSVC 14.3 c++17
		multi::extensions_t<3>
	#endif
		{10, 20, 30}
	);

	BOOST_TEST((  25 % extensions(arr) == decltype(  25 % extensions(arr)){0, 0, 25} ));
	BOOST_TEST((  55 % extensions(arr) == decltype(  55 % extensions(arr))(0, 1, 25) ));
	BOOST_TEST(( 655 % extensions(arr) == decltype( 655 % extensions(arr))(1, 1, 25) ));
	BOOST_TEST((1255 % extensions(arr) == decltype(1255 % extensions(arr))(2, 1, 25) ));

	auto const point = arr.extensions().from_linear(655);
	//  BOOST_TEST( p == std::make_tuple(1, 1, 25) );
	using multi::detail::get;
	BOOST_TEST( get<0>(point) ==  1 );
	BOOST_TEST( get<1>(point) ==  1 );
	BOOST_TEST( get<2>(point) == 25 );
}

BOOST_AUTO_TEST_CASE(layout_tuple_2d) {
	multi::extensions_t<2> const x1({51, 52});
	multi::extensions_t<2> const x2({multi::iextension(0, 51), multi::iextension(0, 52)});

	BOOST_TEST( x1 == x2 );

	multi::extensions_t<2> const x3(std::make_tuple(multi::iextension{0, 51}, multi::iextension{0, 52}));
	BOOST_TEST( x1 == x3 );

	multi::extensions_t<2> const x4 = std::make_tuple(multi::iextension{0, 51}, multi::iextension{0, 52});
	BOOST_TEST( x1 == x4 );

	multi::extensions_t<2> const x5 = std::tuple{
		multi::iextension{0, 51},
		multi::iextension{0, 52}
	};
	BOOST_TEST( x1 == x5 );

	multi::extensions_t<2> const x6 = std::tuple{51, 52};
	BOOST_TEST( x1 == x6 );

	multi::extensions_t<2> const x7{51, 52};
	BOOST_TEST( x1 == x7 );

	multi::extensions_t<2> const x8 = {51, 52};
	BOOST_TEST( x1 == x8 );

	auto const x9 = multi::extensions_t<2>{51, 52};
	BOOST_TEST( x1 == x9 );

	// multi::extensions_t x10{51, 52, 53};  // TODO(correaa) should it work?
	// BOOST_TEST( x1 == x10 );
}

BOOST_AUTO_TEST_CASE(layout_tuple_3d) {
	multi::extensions_t<3> const x1({51, 52, 53});
	multi::extensions_t<3> const x2({
		multi::iextension{0, 51},
		multi::iextension{0, 52},
		multi::iextension{0, 53}
	});
	BOOST_TEST( x1 == x2 );

	multi::extensions_t<3> const x3(std::make_tuple(multi::iextension{0, 51}, multi::iextension{0, 52}, multi::iextension{0, 53}));
	BOOST_TEST( x1 == x3 );

	multi::extensions_t<3> const x4 = std::make_tuple(multi::iextension{0, 51}, multi::iextension{0, 52}, multi::iextension{0, 53});
	BOOST_TEST( x1 == x4 );

	multi::extensions_t<3> const x5 = std::tuple{
		multi::iextension{0, 51},
		multi::iextension{0, 52},
		multi::iextension{0, 53}
	};
	BOOST_TEST( x1 == x5 );

	multi::extensions_t<3> const x6 = std::tuple{51, 52, 53};
	BOOST_TEST( x1 == x6 );

	multi::extensions_t<3> const x7{51, 52, 53};
	BOOST_TEST( x1 == x7 );

	// multi::extensions_t x8{51, 52, 53};  // TODO(correaa) should it work?
	// BOOST_TEST( x1 == x8 );
}

BOOST_AUTO_TEST_CASE(layout_0) {
	multi::array<double, 3> arr(
	#ifdef _MSC_VER  // problem with MSVC 14.3 c++17
		multi::extensions_t<3>
	#endif
		{51, 52, 53}
	);

	BOOST_TEST( size(arr)  == 51       );
	BOOST_TEST( arr.size() == 51       );

	BOOST_TEST( size(arr[0])  == 52    );
	BOOST_TEST( arr[0].size() == 52    );

	BOOST_TEST( size(arr[0][0])  == 53 );
	BOOST_TEST( arr[0][0].size() == 53 );
}

BOOST_AUTO_TEST_CASE(layout_1) {
	// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays): testing feature
	double arr[25][25][25];  // this can overflow the stack: double arr[50][50][50];  50*50*50*8bytes = 1MB

	using multi::size;
	BOOST_TEST( size(arr) == 25 );

	using multi::extension;

	BOOST_TEST(( extension(arr) == multi::index_extension{0, 25} ));
	BOOST_TEST(( extension(arr) == multi::iextension     {0, 25} ));
	// BOOST_TEST(( extension(arr) == multi::irange{0, 25} ));
}

BOOST_AUTO_TEST_CASE(layout_2) {
	std::array<std::array<std::array<double, 25>, 25>, 25> const arr{};
	using multi::size;
	BOOST_TEST( size(arr) == 25 );

	using multi::extension;
	BOOST_TEST(( extension(arr) == multi::index_extension{0, 25} ));
	BOOST_TEST(( extension(arr) == multi::iextension     {0, 25} ));
}

BOOST_AUTO_TEST_CASE(layout_3) {
	multi::array<double, 2> arr(
	#ifdef _MSC_VER  // problem with MSVC 14.3 c++17
		multi::extensions_t<2>
	#endif
		{50, 50}
	);
	BOOST_TEST( size(arr)  == 50 );
	BOOST_TEST( arr.size() == 50 );

	BOOST_TEST( arr[0].sliced(10, 20).size() == 10 );
	BOOST_TEST( size(arr[0].sliced(10, 20))  == 10 );

	static_assert(decltype(arr(0, {10, 20}))::rank_v == 1);

	BOOST_TEST( size(arr(0, {10, 20})) == 10 );

	BOOST_TEST(      arr.layout() == arr.layout()  );
	BOOST_TEST( !(arr.layout() <  arr.layout()) );
}

BOOST_AUTO_TEST_CASE(layout_AA) {
	multi::array<int, 2> const A2 = {
		{1, 2, 3},
		{4, 5, 6},
		{7, 8, 9},
	};

#if defined(__cpp_lib_ranges) && (__cpp_lib_ranges >= 201911L)
#if !defined(__clang_major__) || ((__clang_major__ < 14) && (__clang_major__ != 10))
#if !defined(__NVCC__)
#if !defined(_MSC_VER)
	static_assert( std::ranges::random_access_range<decltype(A2.extension())> );

	auto tiA2 = std::views::transform(
		// A2.extension(),
		std::views::iota(0L, A2.size()),
		[](auto idx) noexcept {return idx;}
	);
	BOOST_TEST( *tiA2.begin() == 0 );
	BOOST_TEST( tiA2[0] == 0 );
#endif
#endif
#endif
#endif

	BOOST_TEST( size(A2) == 3 );

	multi::array<int, 2> B2(
	#ifdef _MSC_VER  // problem with MSVC 14.3 c++17
		multi::extensions_t<2>
	#endif
		{4, 4}
	);

	BOOST_TEST( size(B2) == 4 );
	B2[3][3] = 99;

	multi::array<int, 2> B2copy{B2({0, 2}, {0, 2})};

	auto B2copy2 = B2({0, 2}, {0, 2}).decay();

	BOOST_TEST( &B2copy[1][1] != &B2({0, 2}, {0, 2})[1][1] );

	// clang-format off
	std::array<std::array<decltype(B2({0, 2}, {0, 2})), 2>, 2> B2blk = {{
		{{B2({0, 2}, {0, 2}), B2({0, 2}, {2, 4})}},
		{{B2({2, 4}, {0, 2}), B2({2, 4}, {2, 4})}},
	}};
	// clang-format on

	BOOST_TEST( &B2blk[1][1][1][1] == &B2[3][3] );
}

// BOOST_AUTO_TEST_CASE(layout_BB) {
//  {
//      // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
//      double arr[3][4][5] = {};
//      using multi::dimensionality;
//      static_assert(dimensionality(arr) == 3);
//      using multi::extensions;
//      auto xA = extensions(arr);

//      BOOST_TEST( size(std::get<0>(xA)) == 3 );
//      BOOST_TEST( size(std::get<1>(xA)) == 4 );
//      BOOST_TEST( size(std::get<2>(xA)) == 5 );

//      static_assert(multi::stride(arr) == 20);

//      static_assert(multi::stride(arr[1]) == 5);
//      static_assert(multi::stride(arr[0][0]) == 1);

//      multi::array<double, 3> AA({3, 4, 5});
//      using multi::layout;
//      BOOST_TEST( layout(AA) == layout(arr) );

//      BOOST_TEST( AA.stride() == 20 );
//  }
//  {
//      std::array<std::array<std::array<double, 5>, 4>, 3> arr = {};
//      static_assert(multi::dimensionality(arr) == 3);

//      using multi::extensions;
//      auto xA = extensions(arr);
//      using std::get;
//      BOOST_TEST( size(std::get<0>(xA)) == 3 );
//      BOOST_TEST( size(std::get<1>(xA)) == 4 );
//      BOOST_TEST( size(std::get<2>(xA)) == 5 );

//      multi::array<double, 3> AA({3, 4, 5});
//      using multi::layout;
//      BOOST_TEST( layout(AA) == layout(arr) );

//      BOOST_TEST( AA.stride() == 20 );

//      static_assert(multi::stride(arr) == 20);

//      BOOST_TEST( multi::stride(arr[0])    == 5 );
//      BOOST_TEST( multi::stride(arr[1])    == 5 );
//      BOOST_TEST( multi::stride(arr[0][0]) == 1 );
//  }
//  {
//      multi::array<double, 2> const B2 = {
//          {1.0},
//          {2.0},
//          {3.0},
//      };
//      BOOST_TEST( size(B2) == 3 );
//      BOOST_TEST( B2.rotated().size() == 1 );
//      BOOST_TEST( size(B2[0]) == 1);
//      BOOST_TEST( B2   .stride() == 1 );
//      BOOST_TEST( B2[0].stride() == 1 );
//  }
// }

// BOOST_AUTO_TEST_CASE(multi_layout_with_offset) {
//  static_assert( std::is_trivially_default_constructible_v< multi::layout_t<0> > );
//  static_assert( std::is_trivially_default_constructible_v< multi::layout_t<1> > );
//  static_assert( std::is_trivially_default_constructible_v< multi::layout_t<2> > );

//  static_assert( std::is_trivially_copyable_v< multi::layout_t<2> > );

//  {
//      multi::layout_t<1> const l1(multi::iextension(2, 5));
//      BOOST_TEST( l1.extension().first()  == 2 );
//      BOOST_TEST( l1.extension().last() == 5 );
//  }
//  {
//      boost::multi::layout_t<2>::extensions_type const exts{
//          multi::iextension(2, 5),
//          multi::iextension(0, 5)
//      };
//      multi::layout_t<2> const l2(exts);
//      BOOST_TEST( l2.extension().first()  == std::get<0>(exts).first()  );
//      BOOST_TEST( l2.extension().last () == std::get<0>(exts).last() );
//  }
//  {
//      multi::layout_t<2> const l2({multi::iextension(0, 3), multi::iextension(2, 7)});
//      BOOST_TEST( std::get<1>(l2.extensions()).first()  == 2 );
//      BOOST_TEST( std::get<1>(l2.extensions()).last() == 7 );
//  }
// }

// BOOST_AUTO_TEST_CASE(multi_layout_part1) {
//  {
//      multi::layout_t<0> const lyt;
//      static_assert(decltype(lyt)::rank_v == 0);
//      BOOST_TEST( num_elements(lyt) == 1 );
//  }
//  {
//      multi::iextensions<0> const exts{};
//      multi::layout_t<0> const    lyt(exts);
//      BOOST_TEST(lyt.num_elements() == 1);
//  }
//  {
//      multi::layout_t<1> const lyt{};
//      static_assert(decltype(lyt)::rank_v == 1);
//      BOOST_TEST( num_elements(lyt) == 0 );
//      BOOST_TEST( size(lyt) == 0 );
//      BOOST_TEST( size(extension(lyt)) == 0 );
//      BOOST_TEST( stride(lyt) != 0 );
//      BOOST_TEST( is_empty(lyt) );
//  }
//  {
//      multi::layout_t<2> const lyt({2, 10});
//      static_assert(decltype(lyt)::rank_v == 2);
//      BOOST_TEST( num_elements(lyt) == 20 );
//      BOOST_TEST( size(lyt) == 2 );
//      BOOST_TEST( size(extension(lyt)) == 2 );
//      BOOST_TEST( stride(lyt) == 10 );
//      BOOST_TEST( !is_empty(lyt) );
//  }
//  {
//      multi::layout_t<1> const lyt(multi::iextensions<1>{20});
//      static_assert(decltype(lyt)::rank_v == 1);
//      BOOST_TEST( num_elements(lyt) == 20 );
//      BOOST_TEST( size(lyt) == 20 );
//      BOOST_TEST( stride(lyt) == 1 );
//  }
// }

// BOOST_AUTO_TEST_CASE(multi_layout_part2) {
//  {
//      multi::layout_t<1> const lyt(multi::iextensions<1>{1});
//      static_assert(decltype(lyt)::rank_v == 1);
//      BOOST_TEST( num_elements(lyt) == 1 );
//      BOOST_TEST( size(lyt) == 1 );
//      BOOST_TEST( stride(lyt) == 1 );
//  }
//  {
//      multi::layout_t<2> const lyt({1, 10});
//      static_assert(decltype(lyt)::rank_v == 2);
//      BOOST_TEST( num_elements(lyt) == 10 );
//      BOOST_TEST( size(lyt) == 1);
//      BOOST_TEST( !is_empty(lyt) );
//      BOOST_TEST( size(extension(lyt)) == 1 );
//      BOOST_TEST( stride(lyt) == 10 );  // std::numeric_limits<std::ptrdiff_t>::max() );

//      using std::get;
//      BOOST_TEST( get<0>(strides(lyt)) == 10);
//      BOOST_TEST( get<1>(strides(lyt)) == 1 );
//  }
// }

// BOOST_AUTO_TEST_CASE(multi_layout_part3) {
//  {
//      multi::layout_t<2> const lyt({10, 1});
//      static_assert(decltype(lyt)::rank_v == 2);
//      BOOST_TEST( num_elements(lyt) == 10 );
//      BOOST_TEST( size(lyt) == 10 );
//      using std::get;
//      BOOST_TEST( get<0>(strides(lyt)) == 1 );
//      BOOST_TEST( get<1>(strides(lyt)) == 1 );
//  }
//  {
//      multi::layout_t<2> const lyt{};
//      BOOST_TEST( dimensionality(lyt) == 2 );
//      BOOST_TEST( num_elements(lyt) == 0 );
//      BOOST_TEST( size(lyt) == 0 );
//      BOOST_TEST( size(extension(lyt)) == 0 );
//      BOOST_TEST( stride(lyt) != 0 );
//      BOOST_TEST( is_empty(lyt) );
//  }
//  {
//      multi::layout_t<3> const lyt{};
//      BOOST_TEST( num_elements(lyt) == 0 );
//  }
//  {
//      multi::layout_t<3> const lyt({
//          {0, 10},
//          {0, 10},
//          {0, 10},
//      });
//      BOOST_TEST( num_elements(lyt) == 1000 );
//  }
//  {
//      multi::layout_t<3> const lyt({{10}, {10}, {10}});
//      BOOST_TEST( num_elements(lyt) == 1000 );
//  }
//  {
//      multi::layout_t<3> const lyt({10, 10, 10});
//      BOOST_TEST( num_elements(lyt) == 1000 );
//  }
//  {
//      multi::layout_t<3> const lyt({
//          multi::index_extension{0, 10},
//          {0, 10},
//          {0, 10},
//      });
//      BOOST_TEST( num_elements(lyt) == 1000 );
//  }
//  {
//      multi::layout_t<3> const lyt(multi::layout_t<3>::extensions_type{
//          {0, 10},
//          {0, 10},
//          {0, 10},
//      });
//      BOOST_TEST( num_elements(lyt) == 1000 );
//  }
// }

// BOOST_AUTO_TEST_CASE(layout_to_offset) {
//  multi::layout_t<3> const      lyt({10, 20, 30});
//  multi::array<double, 3> const arr({10, 20, 30});
//  BOOST_TEST( lyt[0][0][0] == &arr[0][0][0] - arr.data_elements() );
//  BOOST_TEST( lyt[0][0][1] == &arr[0][0][1] - arr.data_elements() );
//  BOOST_TEST( lyt[0][0][2] == &arr[0][0][2] - arr.data_elements() );

//  BOOST_TEST_REQUIRE(lyt[0][1][2] == &arr[0][1][2] - arr.data_elements());
//  BOOST_TEST_REQUIRE(lyt[3][1][2] == &arr[3][1][2] - arr.data_elements());
// }

// BOOST_AUTO_TEST_CASE(layout_to_offset_sub) {
//  multi::array<double, 3> arr({10, 20, 30});

//  auto&& sub = arr({2, 6}, {4, 8}, {10, 20});

//  auto const lyt = sub.layout();

//  BOOST_TEST( lyt[0][0][0] == &sub[0][0][0] - base(sub) );
//  BOOST_TEST( lyt[0][0][1] == &sub[0][0][1] - base(sub) );
//  BOOST_TEST( lyt[0][0][2] == &sub[0][0][2] - base(sub) );
//  BOOST_TEST( lyt[0][1][2] == &sub[0][1][2] - base(sub) );
//  BOOST_TEST( lyt[3][1][2] == &sub[3][1][2] - base(sub) );
// }

// BOOST_AUTO_TEST_CASE(continued_part1) {
//  {
//      multi::layout_t<3> const lyt(multi::layout_t<3>::extensions_type{
//          {0, 10},
//          {0, 10},
//          {0, 10},
//      });
//      BOOST_TEST( num_elements(lyt) == 1000);
//  }
//  {
//      multi::layout_t<3> const lyt({
//          multi::iextension{0, 10},
//          multi::iextension{0, 10},
//          multi::iextension{0, 10},
//      });
//      BOOST_TEST(lyt.num_elements() == 1000);
//  }
//  {
//      multi::layout_t<3> const lyt({multi::iextension{10}, multi::iextension{10}, multi::iextension{10}});
//      BOOST_TEST( num_elements(lyt) == 1000);
//  }
//  {
//      multi::layout_t<3> const lyt({10, 10, multi::iextension{10}});
//      BOOST_TEST( num_elements(lyt) == 1000 );
//  }
//  {
//      multi::layout_t<1> const lyt;
//      BOOST_TEST( size(lyt) == 0 );
//  }
//  {
//      multi::layout_t<1> lyt({
//          {0, 10},
//      });
//      BOOST_TEST( size(lyt) == 10 );
//      BOOST_TEST( extension(lyt).first() ==  0 );
//      BOOST_TEST( extension(lyt).last () == 10 );

//      lyt.reindex(1);
//      BOOST_TEST( size(lyt) == 10 );
//      BOOST_TEST( extension(lyt).first() ==  1 );
//      BOOST_TEST( extension(lyt).last () == 11 );
//  }
//  {
//      multi::layout_t<2> const lyt;
//      BOOST_TEST( size(lyt) == 0 );
//  }
//  {
//      multi::layout_t<2> lyt(multi::extensions_t<2>({
//          {0, 10},
//          {0, 20},
//      }));
//      BOOST_TEST( size(lyt) == 10 );
//      BOOST_TEST( extension(lyt).first() ==  0 );
//      BOOST_TEST( extension(lyt).last () == 10 );

//      lyt.reindex(1);
//      BOOST_TEST( extension(lyt).first() ==  1 );
//      BOOST_TEST( extension(lyt).last () == 11 );

//      lyt.rotate().reindex(3).unrotate();
//      BOOST_TEST_REQUIRE( extension(lyt).first() ==  1 );
//      BOOST_TEST_REQUIRE( extension(lyt).last () == 11 );

//      BOOST_TEST_REQUIRE( std::get<0>(extensions(lyt)).first() == 1 );
//      BOOST_TEST_REQUIRE( std::get<1>(extensions(lyt)).first() == 3 );
//      BOOST_TEST_REQUIRE( std::get<1>(extensions(lyt)).last () == 23 );
//  }
// }

// BOOST_AUTO_TEST_CASE(continued_part2) {
//  multi::layout_t<3> const lyt({
//      {0, 10},
//      {0, 20},
//      {0, 30},
//  });

//  BOOST_TEST( !lyt.empty() );

//  BOOST_TEST( stride(lyt) == lyt.stride() );
//  BOOST_TEST( offset(lyt) == lyt.offset() );
//  BOOST_TEST( nelems(lyt) == lyt.nelems() );

//  BOOST_TEST( stride(lyt) == 20*30L );
//  BOOST_TEST( offset(lyt) == 0 );
//  BOOST_TEST( nelems(lyt) == 10*20L*30L );

//  BOOST_TEST( lyt.stride() == stride(lyt) );
//  BOOST_TEST( lyt.offset() == offset(lyt) );
//  BOOST_TEST( lyt.nelems() == nelems(lyt) );

//  using boost::multi::detail::get;
//  BOOST_TEST( get<1>(lyt.strides()) == 30     );
//  BOOST_TEST( get<1>(lyt.offsets()) ==  0     );
//  BOOST_TEST( get<1>(lyt.nelemss()) == 20*30L );

//  BOOST_TEST( get<2>(lyt.strides()) ==  1 );
//  BOOST_TEST( get<2>(lyt.offsets()) ==  0 );
//  BOOST_TEST( get<2>(lyt.nelemss()) == 30 );
// }

// BOOST_AUTO_TEST_CASE(continued_part3) {
//  multi::layout_t<3> const lyt({
//      {0, 10},
//      {0, 20},
//      {0, 30},
//  });

//  BOOST_TEST( lyt.num_elements() == num_elements(lyt) );
//  BOOST_TEST( lyt.size() == size(lyt) );
//  BOOST_TEST( lyt.extension() == extension(lyt) );

//  BOOST_TEST( num_elements(lyt) == 10*20L*30L );
//  BOOST_TEST( size(lyt) == 10 );
//  BOOST_TEST( extension(lyt).first() == 0 );
//  BOOST_TEST( extension(lyt).last() == 10 );

//  BOOST_TEST( std::get<0>(lyt.extensions()) == lyt.extension() );

//  boost::multi::extensions_t<2> const exts2;

//  using boost::multi::detail::get;
//  using std::get;

//  BOOST_TEST( get<0>(exts2).is_empty() );

//  //  BOOST_TEST( std::get<0>(L.sizes()) == L.size(0) );
//  //  BOOST_TEST( std::get<0>(L.extensions()) == L.extension(0) );

//  BOOST_TEST(( get<0>(lyt.extensions()) == multi::index_extension{0, 10} ));

//  BOOST_TEST( get<0>(lyt.extensions()).first() ==  0 );
//  BOOST_TEST( get<0>(lyt.extensions()).last()  == 10 );

//  //  BOOST_TEST( L.size(1) == 20 );
//  BOOST_TEST( get<1>(lyt.extensions()).first() ==  0 );
//  BOOST_TEST( get<1>(lyt.extensions()).last()  == 20 );

//  //  BOOST_TEST( L.size(2) == 30 );
//  BOOST_TEST( get<2>(lyt.extensions()).first() ==  0 );
//  BOOST_TEST( get<2>(lyt.extensions()).last()  == 30 );

//  using std::get;
//  BOOST_TEST( get<0>(strides(lyt)) == lyt.stride() );

//  auto const& strides = lyt.strides();
//  BOOST_TEST( get<0>(strides) == lyt.stride() );
// }

// BOOST_AUTO_TEST_CASE(continued) {
//  {
//      multi::layout_t<3> const lyt;
//      BOOST_TEST( size(lyt) == 0 );
//  }
//  {
//      multi::layout_t<3> const lyt({
//          {0, 10},
//          {0, 20},
//          {0, 30},
//      });
//      BOOST_TEST( stride(lyt) == 20*30L );
//  }
//  {
//      multi::layout_t<1> const lyt({
//          {0, 10},
//      });
//      BOOST_TEST( extension(lyt).first() == 0 );
//      BOOST_TEST( extension(lyt).last() == 10 );
//  }
//  {
//      multi::layout_t<1> const lyt({
//          {8, 18},
//      });
//      BOOST_TEST( extension(lyt).first() == 8 );
//      BOOST_TEST( extension(lyt).last() == 18 );
//  }
//  {
//      multi::layout_t<2> const lyt(multi::extensions_t<2>({
//          {0, 10},
//          {0, 20},
//      }));
//      BOOST_TEST( extension(lyt).first() == 0 );
//      BOOST_TEST( extension(lyt).last() == 10 );
//  }
//  // {  // this is ambiguous in nvcc
//  //  multi::layout_t<2> const lyt({
//  //      {0, 10},
//  //      {0, 20},
//  //  });
//  //  BOOST_TEST( extension(lyt).first() == 0 );
//  //  BOOST_TEST( extension(lyt).last() == 10 );
//  // }
//  {
//      multi::layout_t<2> const lyt(multi::extensions_t<2>({
//          { 0, 10},
//          {11, 31},
//      }));
//      BOOST_TEST( size(lyt) == 10   );
//      BOOST_TEST( stride(lyt) == 20 );
//      BOOST_TEST( offset(lyt) == 0 );
//  }
//  {  // this is ambiguous in nvcc
//      multi::layout_t<2> const lyt(multi::extensions_t<2>({
//          { 0, 10},
//          {11, 31},
//      }));
//      BOOST_TEST( size(lyt) == 10   );
//      BOOST_TEST( stride(lyt) == 20 );
//      BOOST_TEST( offset(lyt) == 0 );
//  }
//  {
//      multi::layout_t<2> const lyt(multi::extensions_t<2>({
//          {8, 18},
//          {0, 20},
//      }));
//      BOOST_TEST( size(lyt) == 10 );
//      BOOST_TEST( stride(lyt) == 20 );
//  }
//  // {
//  //  multi::layout_t<3> const lyt(multi::extensions_t<3>({
//  //      { 0,  3},
//  //      { 0,  5},
//  //      {10, 17},
//  //  }));
//  //  BOOST_TEST( stride(lyt) == 5*7L );
//  //  BOOST_TEST( stride(lyt.sub().sub()) == 1 );
//  // }
//  {
//      multi::layout_t<3> const lyt({
//          {0, 10},
//          {0, 20},
//          {0, 30},
//      });
//      BOOST_TEST( size(lyt) == 10 );
//      BOOST_TEST( stride(lyt) == 20*30L );
//      BOOST_TEST( offset(lyt) == 0 );
//      BOOST_TEST( nelems(lyt) == 10*20L*30L );
//  }
//  {
//      multi::layout_t<3> const lyt({
//          {10, 20},
//          {10, 30},
//          {10, 40},
//      });
//      BOOST_TEST( stride(lyt) == 20*30L );
//  }
//  {
//      auto const ttt = boost::multi::tuple<int, int, int>{1, 2, 3};
//      auto const arr = std::apply([](auto... elems) { return std::array<int, 3>{{elems...}}; }, ttt);
//      BOOST_TEST(arr[1] == 2);
//  }
// }

// BOOST_AUTO_TEST_CASE(extensions_from_linear_1d) {
//  multi::extensions_t<1> const exts{11};

//  auto ijk = exts.from_linear(9);

//  using multi::detail::get;
//  BOOST_TEST_REQUIRE( get<0>(ijk) == 9 );

//  multi::layout_t<1> const lyt{exts};
//  BOOST_TEST_REQUIRE( lyt[get<0>(ijk)] == 9 );
//  BOOST_TEST_REQUIRE( lyt(get<0>(ijk)) == 9 );

//  BOOST_TEST_REQUIRE( lyt(std::get<0>(lyt.extensions().from_linear(9))) == 9 );

//  BOOST_TEST_REQUIRE( std::apply(lyt, lyt.extensions().from_linear(9)) == 9 );
// }

// BOOST_AUTO_TEST_CASE(extensions_from_linear_2d_structured_binding) {
//  multi::extensions_t<2> const exts{3, 5};
//  auto [eye, jay] = exts.from_linear(7);

//  BOOST_TEST_REQUIRE( eye == 1 );
//  BOOST_TEST_REQUIRE( jay == 2 );
//  //  BOOST_TEST_REQUIRE( std::apply(l, l.extensions().from_linear(9)) == 9 );
// }

// BOOST_AUTO_TEST_CASE(extensions_from_linear_2d_std_get) {
//  multi::extensions_t<2> const exts{3, 5};
//  auto                         eye = std::get<0>(exts.from_linear(7));
//  auto                         jay = std::get<1>(exts.from_linear(7));
//  BOOST_TEST_REQUIRE( eye == 1 );
//  BOOST_TEST_REQUIRE( jay == 2 );
// }

// BOOST_AUTO_TEST_CASE(extensions_from_linear_2d_std_get_using) {
//  multi::extensions_t<2> const exts{3, 5};
//  using std::get;
//  auto       fl  = exts.from_linear(7L);
//  auto const eye = get<0>(fl);
//  auto const jay = get<1>(fl);
//  BOOST_TEST_REQUIRE( eye == 1 );
//  BOOST_TEST_REQUIRE( jay == 2 );
// }

// BOOST_AUTO_TEST_CASE(extensions_from_linear_2d_get_using) {
//  multi::extensions_t<2> const exts{3, 5};

//  using multi::detail::get;

//  auto eye = get<0>(exts.from_linear(7));
//  auto jay = get<1>(exts.from_linear(7));
//  BOOST_TEST_REQUIRE( eye == 1 );
//  BOOST_TEST_REQUIRE( jay == 2 );
// }

// BOOST_AUTO_TEST_CASE(extensions_from_linear_2d) {
//  multi::extensions_t<2> const exts{3, 5};

//  auto ij = exts.from_linear(7);

//  using multi::detail::get;

//  BOOST_TEST_REQUIRE( get<0>(ij) == 1 );
//  BOOST_TEST_REQUIRE( get<1>(ij) == 2 );

//  multi::layout_t<2> const lyt{exts};
//  BOOST_TEST_REQUIRE( lyt[get<0>(ij)][get<1>(ij)] == 7 );
// }

// BOOST_AUTO_TEST_CASE(extensions_from_linear_3d_std_get) {
//  multi::extensions_t<3> const exts{11, 13, 17};

//  BOOST_TEST_REQUIRE( std::get<0>(exts.from_linear( 0)) ==  0 );
//  BOOST_TEST_REQUIRE( std::get<1>(exts.from_linear( 0)) ==  0 );
//  BOOST_TEST_REQUIRE( std::get<2>(exts.from_linear( 0)) ==  0 );

//  BOOST_TEST_REQUIRE( std::get<0>(exts.from_linear( 1)) ==  0 );
//  BOOST_TEST_REQUIRE( std::get<1>(exts.from_linear( 1)) ==  0 );
//  BOOST_TEST_REQUIRE( std::get<2>(exts.from_linear( 1)) ==  1 );

//  BOOST_TEST_REQUIRE( std::get<0>(exts.from_linear(16)) ==  0 );
//  BOOST_TEST_REQUIRE( std::get<1>(exts.from_linear(16)) ==  0 );
//  BOOST_TEST_REQUIRE( std::get<2>(exts.from_linear(16)) == 16 );

//  BOOST_TEST_REQUIRE( std::get<0>(exts.from_linear(17)) ==  0 );
//  BOOST_TEST_REQUIRE( std::get<1>(exts.from_linear(17)) ==  1 );
//  BOOST_TEST_REQUIRE( std::get<2>(exts.from_linear(17)) ==  0 );

//  BOOST_TEST_REQUIRE( std::get<0>(exts.from_linear(18)) ==  0 );
//  BOOST_TEST_REQUIRE( std::get<1>(exts.from_linear(18)) ==  1 );
//  BOOST_TEST_REQUIRE( std::get<2>(exts.from_linear(18)) ==  1 );

//  multi::layout_t<3> const lyt{exts};

//  using std::get;
//  BOOST_TEST_REQUIRE( lyt[get<0>(exts.from_linear(19))][get<1>(exts.from_linear(19))][get<2>(exts.from_linear(19))] == 19 );
//  BOOST_TEST_REQUIRE( lyt(get<0>(exts.from_linear(19)), get<1>(exts.from_linear(19)), get<2>(exts.from_linear(19))) == 19 );
// }

// BOOST_AUTO_TEST_CASE(extensions_from_linear_3d_std_get_using) {
//  multi::extensions_t<3> const exts{11, 13, 17};

//  using std::get;

//  BOOST_TEST_REQUIRE( get<0>(exts.from_linear( 0)) ==  0 );
//  BOOST_TEST_REQUIRE( get<1>(exts.from_linear( 0)) ==  0 );
//  BOOST_TEST_REQUIRE( get<2>(exts.from_linear( 0)) ==  0 );

//  BOOST_TEST_REQUIRE( get<0>(exts.from_linear( 1)) ==  0 );
//  BOOST_TEST_REQUIRE( get<1>(exts.from_linear( 1)) ==  0 );
//  BOOST_TEST_REQUIRE( get<2>(exts.from_linear( 1)) ==  1 );

//  BOOST_TEST_REQUIRE( get<0>(exts.from_linear(16)) ==  0 );
//  BOOST_TEST_REQUIRE( get<1>(exts.from_linear(16)) ==  0 );
//  BOOST_TEST_REQUIRE( get<2>(exts.from_linear(16)) == 16 );

//  BOOST_TEST_REQUIRE( get<0>(exts.from_linear(17)) ==  0 );
//  BOOST_TEST_REQUIRE( get<1>(exts.from_linear(17)) ==  1 );
//  BOOST_TEST_REQUIRE( get<2>(exts.from_linear(17)) ==  0 );

//  BOOST_TEST_REQUIRE( get<0>(exts.from_linear(18)) ==  0 );
//  BOOST_TEST_REQUIRE( get<1>(exts.from_linear(18)) ==  1 );
//  BOOST_TEST_REQUIRE( get<2>(exts.from_linear(18)) ==  1 );

//  BOOST_TEST_REQUIRE( get<0>(exts.from_linear(19)) ==  0 );
//  BOOST_TEST_REQUIRE( get<1>(exts.from_linear(19)) ==  1 );
//  BOOST_TEST_REQUIRE( get<2>(exts.from_linear(19)) ==  2 );

//  multi::layout_t<3> const lyt{exts};
//  BOOST_TEST_REQUIRE( lyt[get<0>(exts.from_linear(19))][get<1>(exts.from_linear(19))][get<2>(exts.from_linear(19))] == 19 );
//  BOOST_TEST_REQUIRE( lyt(get<0>(exts.from_linear(19)), get<1>(exts.from_linear(19)), get<2>(exts.from_linear(19))) == 19 );
// }

// BOOST_AUTO_TEST_CASE(extensions_from_linear_3d_struct_bind) {
//  multi::extensions_t<3> const exts{11, 13, 17};

//  using std::get;
//  {
//      auto [eye, jay, kay] = exts.from_linear(0);
//      BOOST_TEST_REQUIRE(eye == 0);
//      BOOST_TEST_REQUIRE(jay == 0);
//      BOOST_TEST_REQUIRE(kay == 0);
//  }
//  {
//      auto [eye, jay, kay] = exts.from_linear(1);
//      BOOST_TEST_REQUIRE(eye == 0);
//      BOOST_TEST_REQUIRE(jay == 0);
//      BOOST_TEST_REQUIRE(kay == 1);
//  }
//  {
//      auto [eye, jay, kay] = exts.from_linear(16);
//      BOOST_TEST_REQUIRE(eye == 0);
//      BOOST_TEST_REQUIRE(jay == 0);
//      BOOST_TEST_REQUIRE(kay == 16);
//  }
//  {
//      auto [eye, jay, kay] = exts.from_linear(17);
//      BOOST_TEST_REQUIRE(eye == 0);
//      BOOST_TEST_REQUIRE(jay == 1);
//      BOOST_TEST_REQUIRE(kay == 0);
//  }
//  {
//      auto [eye, jay, kay] = exts.from_linear(18);
//      BOOST_TEST_REQUIRE(eye == 0);
//      BOOST_TEST_REQUIRE(jay == 1);
//      BOOST_TEST_REQUIRE(kay == 1);

//      multi::layout_t<3> const lyt{exts};
//      BOOST_TEST_REQUIRE(lyt[eye][jay][kay] == 18);
//      BOOST_TEST_REQUIRE(lyt(eye, jay, kay) == 18);
//  }
// }

// BOOST_AUTO_TEST_CASE(extensions_from_linear_3d) {
//  multi::extensions_t<3> const exts{11, 13, 17};

//  auto ijk = exts.from_linear(19);

//  {
//      using std::get;
//      BOOST_TEST_REQUIRE(get<0>(exts.from_linear(19)) == 0);
//      BOOST_TEST_REQUIRE(get<1>(exts.from_linear(19)) == 1);
//      BOOST_TEST_REQUIRE(get<2>(exts.from_linear(19)) == 2);
//  }
//  {
//      using std::get;
//      //  using multi::detail::get;
//      BOOST_TEST_REQUIRE(get<0>(ijk) == 0);
//      BOOST_TEST_REQUIRE(get<1>(ijk) == 1);
//      BOOST_TEST_REQUIRE(get<2>(ijk) == 2);

//      multi::layout_t<3> const lyt{exts};

//      BOOST_TEST_REQUIRE(lyt[get<0>(ijk)][get<1>(ijk)][get<2>(ijk)] == 19);
//      BOOST_TEST_REQUIRE(lyt(get<0>(ijk), get<1>(ijk), get<2>(ijk)) == 19);
//  }
// }

// BOOST_AUTO_TEST_CASE(extension_1D_iteration) {
//  multi::extension_t const ext(10);
//  BOOST_TEST_REQUIRE(ext[0] == 0);
//  BOOST_TEST_REQUIRE(ext[1] == 1);
// }

// BOOST_AUTO_TEST_CASE(extensionS_1D_iteration) {
//  {
//      multi::extensions_t<1> const exts(10);
//      BOOST_TEST_REQUIRE(std::get<0>(exts[0]) == 0);
//      BOOST_TEST_REQUIRE(std::get<0>(exts[1]) == 1);
//  }
//  {
//      multi::extensions_t<1> const exts(multi::iextension{0, 10});
//      BOOST_TEST_REQUIRE(std::get<0>(exts[0]) == 0);
//      BOOST_TEST_REQUIRE(std::get<0>(exts[1]) == 1);
//  }
// }

// BOOST_AUTO_TEST_CASE(extensionS_2D_iteration) {
//  {
//      multi::extensions_t<2> exts({3, 5});
//      BOOST_TEST_REQUIRE(std::get<0>(exts[0]) == 0);
//      BOOST_TEST_REQUIRE(std::get<0>(exts[1]) == 1);
//  }
//  {
//      multi::extensions_t<2> exts({multi::iextension{0, 3}, multi::iextension{0, 5}});
//      BOOST_TEST_REQUIRE(std::get<0>(exts[0]) == 0);
//      BOOST_TEST_REQUIRE(std::get<0>(exts[1]) == 1);
//  }
// }

// BOOST_AUTO_TEST_CASE(layout_1D_iteration) {
//  multi::layout_t<1> const lyt{multi::extensions_t<1>(10)};
//  BOOST_TEST( lyt[0] == 0 );
//  BOOST_TEST( lyt[1] == 1 );
//  BOOST_TEST( lyt[2] == 2 );

//  //  BOOST_TEST_REQUIRE(std::get<0>(exts[0]) == 0);
//  //  BOOST_TEST_REQUIRE(std::get<0>(exts[1]) == 1);
// }

// BOOST_AUTO_TEST_CASE(layout_2D_iteration) {
//  multi::layout_t<2> const lyt{multi::extensions_t<2>({5, 3})};
//  BOOST_TEST( lyt[0][0] == 0 );
//  BOOST_TEST( lyt[0][1] == 1 );
//  BOOST_TEST( lyt[0][2] == 2 );

//  BOOST_TEST( lyt[1][0] == 3 );
//  BOOST_TEST( lyt[1][1] == 4 );
//  BOOST_TEST( lyt[1][2] == 5 );

//  //  BOOST_TEST_REQUIRE(std::get<0>(exts[0]) == 0);
//  //  BOOST_TEST_REQUIRE(std::get<0>(exts[1]) == 1);
// }

return boost::report_errors();
}
