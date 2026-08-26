// Copyright 2018-2025 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for range, layout_t, get, extents_t

#include <boost/core/lightweight_test.hpp>

#include <algorithm>  // for copy
#include <array>      // for array, array<>::value_type
#include <cstddef>    // for ptrdiff_t, size_t  // IWYU pragma: keep
#include <iterator>   // for size

#if (__cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L)) && __has_include(<ranges>)
#if !defined(__clang_major__) || (__clang_major__ != 16)
#include <ranges>  // IWYU pragma: keep
#endif
#endif

#include <tuple>   // for make_tuple, tuple_element<>::type
#include <vector>  // for vector
// IWYU pragma: no_include <version>
#include <type_traits>

namespace multi = boost::multi;

// all layout types must stay trivial (trivially copyable + trivial default ctor) so default
// construction is free and the types remain cheap/memcpy-able/device-friendly.
// (std::is_trivial is deprecated in C++26, hence the two-trait spelling it recommends)
template<class T> constexpr bool is_trivial_v = std::is_trivially_copyable_v<T> && std::is_trivially_default_constructible_v<T>;
static_assert(is_trivial_v<multi::layout_t<0>>);
static_assert(is_trivial_v<multi::layout_t<1>>);
static_assert(is_trivial_v<multi::layout_t<2>>);
static_assert(is_trivial_v<multi::layout_t<3>>);
static_assert(is_trivial_v<multi::layout_t<4>>);

namespace {
auto second_finish(multi::extents_t<3> exts) {
	using std::get;  // workaround: function call with explicit template arguments is a C++20 extension [-Wc++20-extensions]
	return get<1>(exts).last();
}
}  // namespace

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	// BOOST_AUTO_TEST_CASE(extensions_3D)
	{
		BOOST_TEST( 20 == second_finish( multi::extents_t<3>  { {0, 10}, {0, 20}, {0, 30} }  ) );
		BOOST_TEST( 20 == second_finish( multi::extents_t<3>( { {0, 10}, {0, 20}, {0, 30} } )) );
		BOOST_TEST( 20 == second_finish(                         { {0, 10}, {0, 20}, {0, 30} }  ) );

		multi::extents_t<3> const exts({0, 10}, {0, 20}, {0, 30});
		BOOST_TEST( 20 == second_finish(exts) );
	}

	// BOOST_AUTO_TEST_CASE(extensions_to_linear)
	{
		multi::extents_t<3> exts{4, 5, 3};
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

	// BOOST_AUTO_TEST_CASE(contiguous_layout)
	{
		std::vector<int> vec(10, 99);  // NOLINT(fuchsia-default-arguments-calls)
		using ArrayRef = multi::array_ref<int, 1, int*, multi::contiguous_layout<>>;
		auto arr       = ArrayRef({static_cast<multi::ssize_t>(vec.size())}, vec.data());

		BOOST_TEST( &arr[1] == &vec[1] );

		static_assert(
			std::is_base_of_v<
				std::random_access_iterator_tag, ArrayRef::const_iterator::iterator_category>
		);

		// #if (__cplusplus >= 202002L)
		//      static_assert(
		//          std::is_base_of_v<
		//              std::contiguous_iterator_tag, ArrayRef::const_iterator::iterator_category>
		//      );
		// #endif
	}

	{
		std::vector<int> vec = {1, 2, 3, 4, 5};  // NOLINT(fuchsia-default-arguments-calls)
		std::vector<int> vec2(5);                // NOLINT(fuchsia-default-arguments-calls)

		multi::array_ref<int, 1, int*, multi::contiguous_layout<>> const arr(static_cast<std::ptrdiff_t>(vec.size()), vec.data());

		auto&& arr_d = arr.dropped(1);
		BOOST_TEST(  arr_d.size() == arr.size() - 1 );
		BOOST_TEST( &arr_d[0] == &arr[1] );

		auto&& arr_s = arr.sliced(1, 4);
		BOOST_TEST( &arr_s[0] == &arr[1] );

		static_assert(
			std::is_base_of_v<
				std::random_access_iterator_tag,
				decltype(arr.cbegin())::iterator_category>
		);

#if defined(__cplusplus) && (__cplusplus >= 202002L) && defined(__cpp_lib_ranges) && (!defined(__clang__) || __clang_major__ != 10)
		static_assert(
			std::is_base_of_v<
				std::contiguous_iterator_tag,
				decltype(arr.cbegin())::iterator_category>
		);

		static_assert(std::contiguous_iterator<decltype(arr.begin())>);
		static_assert(std::contiguous_iterator<decltype(arr.end())>);

		static_assert(std::contiguous_iterator<decltype(arr.cbegin())>);
		static_assert(std::contiguous_iterator<decltype(arr.cend())>);

		int const* beg(arr.begin());
		BOOST_TEST( beg == &arr.front() );
#endif

		// std::copy(arr.begin(), arr.end(), vec2.begin());
		BOOST_TEST( arr.cbegin().stride() == 1 );
		BOOST_TEST( arr.cend().stride() == 1 );

		auto size = arr.cend() - arr.cbegin();
		BOOST_TEST( size == 5 );
		BOOST_TEST( arr.size() == 5 );

		BOOST_TEST( arr.cbegin().stride() == 1 );

		BOOST_TEST( arr.cend().base() - arr.cbegin().base() == 5 );

		BOOST_TEST( (arr.cend().base() - arr.cbegin().base()) % arr.cbegin().stride() == 0 );

		// std::copy(arr.data_elements(), std::next(arr.data_elements(), arr.num_elements()), vec2.begin());
		BOOST_TEST( &*arr.cbegin() == arr.data_elements() );

		// vvv this is UB, never dereference an end iterator
		// BOOST_TEST( &*arr.cend() == std::next(arr.data_elements(), arr.num_elements()) );

		BOOST_TEST( *arr.cbegin() == 1 );
		BOOST_TEST( *std::next(arr.cbegin(), 1) == 2 );

		// multi::what(arr.cbegin());
		std::copy(arr.cbegin(), arr.cbegin() + arr.size(), vec2.begin());
		// std::copy_n(arr.begin(), arr.size(), vec2.begin());
		// std::copy(arr.cbegin(), arr.cend(), vec2.begin());
		// std::copy_n(arr.cbegin(), arr.size(), vec2.begin());

		// for(auto idx : arr.extension()) {  // NOLINT(altera-unroll-loops)
		//  vec2[static_cast<std::size_t>(idx)] = arr[idx];
		// }

		BOOST_TEST( vec2[0] == 1 );
		BOOST_TEST( vec2[1] == 2 );
		BOOST_TEST( vec2[2] == 3 );
		BOOST_TEST( vec2[3] == 4 );
		BOOST_TEST( vec2[4] == 5 );

		BOOST_TEST( vec == vec2 );
	}

	{
		multi::array<double, 2> const d2D = {
			{150.0, 16.0, 17.0, 18.0, 19.0},
			{ 30.0,  1.0,  2.0,  3.0,  4.0},
			{100.0, 11.0, 12.0, 13.0, 14.0},
			{ 50.0,  6.0,  7.0,  8.0,  9.0},
		};

		// #if __has_cpp_attribute(no_unique_address) >=201803L and not defined(__NVCC__) and not defined(__PGI)
		//  BOOST_TEST( sizeof(d2D)==sizeof(double*)+7*sizeof(std::size_t) );
		// #endif

		BOOST_TEST( d2D.layout().is_compact() );
		BOOST_TEST( d2D.rotated().layout().is_compact() );
		BOOST_TEST( d2D[3].layout().is_compact() );
		BOOST_TEST( !(d2D.rotated()[2].layout().is_compact()) );
	}
	{
		multi::array<int, 2> d2D({5, 3});
		BOOST_TEST( d2D.layout().is_compact() );
		BOOST_TEST( d2D.rotated().layout().is_compact() );
		BOOST_TEST( d2D[3].layout().is_compact() );
		BOOST_TEST( !d2D.rotated()[2].layout().is_compact() );
	}

	// BOOST_AUTO_TEST_CASE(extensions_layout_to_linear)
	{
		multi::array<double, 3> arr({40, 50, 80});

		auto&& sub = arr({10, 30}, {20, 32}, {60, 75});

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
#endif

		for(int i = 0; i != 10; ++i) {
			for(int j = 0; j != 12; ++j) {
				for(int k = 0; k != 15; ++k) {                                    // NOLINT(altera-unroll-loops)
					BOOST_TEST( &  sub.base()  [sub.layout()(i, j, k)] == &sub(i, j, k) );    // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
					BOOST_TEST( &*(sub.base() + sub.layout()(i, j, k)) == &sub(i, j, k) );  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
				}
			}
		}

#ifdef __clang__
#pragma clang diagnostic pop
#endif
	}

	// BOOST_AUTO_TEST_CASE(extensions_layout_to_linear_2)
	{
		multi::array<double, 3> arr(
#ifdef _MSC_VER  // problem with MSVC 14.3 c++17
			multi::extents_t<3>
#endif
			{40, 50, 80}
		);

		auto&& sub = arr({10, 30}, {20, 32}, {60, 75});

		auto const& rot = sub.rotated();

		auto const [is, js, ks] = rot.extents();

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif

		for(auto const i : is) {
			for(auto const j : js) {
				for(auto const k : ks) {                                          // NOLINT(altera-unroll-loops)
					BOOST_TEST( &  rot.base()  [rot.layout()(i, j, k)] == &rot(i, j, k) );    // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
					BOOST_TEST( &*(rot.base() + rot.layout()(i, j, k)) == &rot(i, j, k) );  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
				}
			}
		}

#ifdef __clang__
#pragma clang diagnostic pop
#endif
	}

	// BOOST_AUTO_TEST_CASE(linearize)
	{
		multi::array<double, 3> const arr(
#ifdef _MSC_VER  // problem with MSVC 14.3 c++17
			multi::extents_t<3>
#endif
			{10, 20, 30}
		);

		BOOST_TEST((  25 % extents(arr) == decltype(  25 % extents(arr)){0, 0, 25} ));
		BOOST_TEST((  55 % extents(arr) == decltype(  55 % extents(arr))(0, 1, 25) ));
		BOOST_TEST(( 655 % extents(arr) == decltype( 655 % extents(arr))(1, 1, 25) ));
		BOOST_TEST((1255 % extents(arr) == decltype(1255 % extents(arr))(2, 1, 25) ));

		auto const point = arr.extents().from_linear(655);
		//  BOOST_TEST( p == std::make_tuple(1, 1, 25) );
		using multi::detail::get;
		BOOST_TEST( get<0>(point) ==  1 );
		BOOST_TEST( get<1>(point) ==  1 );
		BOOST_TEST( get<2>(point) == 25 );
	}

	// BOOST_AUTO_TEST_CASE(layout_tuple_2d)
	{
		multi::extents_t<2> const x1({51, 52});
		multi::extents_t<2> const x2({multi::iextension(0, 51), multi::iextension(0, 52)});

		BOOST_TEST( x1 == x2 );

		multi::extents_t<2> const x3(std::make_tuple(multi::iextension{0, 51}, multi::iextension{0, 52}));
		BOOST_TEST( x1 == x3 );

		multi::extents_t<2> const x4 = std::make_tuple(multi::iextension{0, 51}, multi::iextension{0, 52});
		BOOST_TEST( x1 == x4 );

		multi::extents_t<2> const x5 = std::tuple{
			multi::iextension{0, 51},
			multi::iextension{0, 52},
		};
		BOOST_TEST( x1 == x5 );

		multi::extents_t<2> const x6 = std::tuple{51, 52};
		BOOST_TEST( x1 == x6 );

		multi::extents_t<2> const x7{51, 52};
		BOOST_TEST( x1 == x7 );

		multi::extents_t<2> const x8 = {51, 52};
		BOOST_TEST( x1 == x8 );

		auto const x9 = multi::extents_t<2>{51, 52};
		BOOST_TEST( x1 == x9 );

		// multi::extents_t x10{51, 52, 53};  // TODO(correaa) should it work?
		// BOOST_TEST( x1 == x10 );
	}

	// BOOST_AUTO_TEST_CASE(layout_tuple_3d)
	{
		multi::extents_t<3> const x1({51, 52, 53});
		multi::extents_t<3> const x2({
			multi::iextension{0, 51},
			multi::iextension{0, 52},
			multi::iextension{0, 53},
		});
		BOOST_TEST( x1 == x2 );

		multi::extents_t<3> const x3(std::make_tuple(multi::iextension{0, 51}, multi::iextension{0, 52}, multi::iextension{0, 53}));
		BOOST_TEST( x1 == x3 );

		multi::extents_t<3> const x4 = std::make_tuple(multi::iextension{0, 51}, multi::iextension{0, 52}, multi::iextension{0, 53});
		BOOST_TEST( x1 == x4 );

		multi::extents_t<3> const x5 = std::tuple{
			multi::iextension{0, 51},
			multi::iextension{0, 52},
			multi::iextension{0, 53},
		};
		BOOST_TEST( x1 == x5 );

		multi::extents_t<3> const x6 = std::tuple{51, 52, 53};
		BOOST_TEST( x1 == x6 );

		multi::extents_t<3> const x7{51, 52, 53};
		BOOST_TEST( x1 == x7 );

		// multi::extents_t x8{51, 52, 53};  // TODO(correaa) should it work?
		// BOOST_TEST( x1 == x8 );
	}

	// BOOST_AUTO_TEST_CASE(layout_0)
	{
		multi::array<double, 3> arr(
#ifdef _MSC_VER  // problem with MSVC 14.3 c++17
			multi::extents_t<3>
#endif
			{51, 52, 53}
		);

		using std::size;

		BOOST_TEST( arr.size()  == 51       );
		BOOST_TEST( arr.size() == 51       );

		BOOST_TEST( arr[0].size()  == 52    );
		BOOST_TEST( arr[0].size() == 52    );

		BOOST_TEST( size(arr[0][0])  == 53 );
		BOOST_TEST( arr[0][0].size() == 53 );
	}

	// BOOST_AUTO_TEST_CASE(layout_1)
	{
		// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays): testing feature
		double arr[25][25][25];  // this can overflow the stack: double arr[50][50][50];  50*50*50*8bytes = 1MB

		using multi::size;
		BOOST_TEST( size(arr) == 25 );

		using multi::extent;

		BOOST_TEST(( extent(arr) == multi::index_extension{0, 25} ));
		BOOST_TEST(( extent(arr) == multi::iextension     {0, 25} ));
		// BOOST_TEST(( extension(arr) == multi::irange{0, 25} ));
	}

	// BOOST_AUTO_TEST_CASE(layout_2)
	{
		std::array<std::array<std::array<double, 25>, 25>, 25> const arr{};
		BOOST_TEST( arr.size() == 25 );

		using multi::extent;
		BOOST_TEST(( extent(arr) == multi::index_extension{0, 25} ));
		BOOST_TEST(( extent(arr) == multi::iextension     {0, 25} ));
	}

	// BOOST_AUTO_TEST_CASE(layout_AA)
	{
		multi::array<int, 2> const A2 = {
			{1, 2, 3},
			{4, 5, 6},
			{7, 8, 9},
		};

#if defined(__cpp_lib_ranges) && (__cpp_lib_ranges >= 201911L)
#if !defined(__clang_major__) || ((__clang_major__ < 14) && (__clang_major__ != 10))
#if !defined(__NVCC__)
#if !defined(_MSC_VER)
		static_assert(std::ranges::random_access_range<decltype(A2.extent())>);

		auto tiA2 = std::views::transform(
			// A2.extension(),
			std::views::iota(0L, A2.size()),
			[](auto idx) noexcept { return idx; }
		);
		BOOST_TEST( *tiA2.begin() == 0 );
		BOOST_TEST( tiA2[0] == 0 );
#endif
#endif
#endif
#endif

		BOOST_TEST( size(A2) == 3 );

		multi::array<int, 2> B2({4, 4}, 5);

		BOOST_TEST( size(B2) == 4 );
		B2[3][3] = 99;

		BOOST_TEST( B2[3][3] == 99 );  // cppcheck-suppress knownConditionTrueFalse ; for test

		multi::array<int, 2> B2copy{B2({0, 2}, {0, 2})};

		BOOST_TEST( &B2copy[1][1] != &B2({0, 2}, {0, 2})[1][1] );

		auto B2copy2 = B2({0, 2}, {0, 2}).decay();

		BOOST_TEST( B2copy2 == B2({0, 2}, {0, 2}) );
		BOOST_TEST( B2copy2.base() != B2({0, 2}, {0, 2}).base() );

		// clang-format off
		std::array<std::array<decltype(B2({0, 2}, {0, 2})), 2>, 2> B2blk = {{
			{{B2({0, 2}, {0, 2}), B2({0, 2}, {2, 4})}},
			{{B2({2, 4}, {0, 2}), B2({2, 4}, {2, 4})}},
		}};
		// clang-format on

		BOOST_TEST( &B2blk[1][1][1][1] == &B2[3][3] );
	}

	// BOOST_AUTO_TEST_CASE(layout_BB)
	{
		// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) test legacy type
		double arr[3][4][5] = {};

		using multi::dimensionality;
		static_assert(dimensionality(arr) == 3);

		using multi::extents;
		auto xA = extents(arr);

		using std::get;  // needed for C++17
		using std::size;
		BOOST_TEST( size(get<0>(xA)) == 3 );
		BOOST_TEST( size(get<1>(xA)) == 4 );
		BOOST_TEST( size(get<2>(xA)) == 5 );
	}

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
	//      multi::layout_t<2> lyt(multi::extents_t<2>({
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

	//  boost::multi::extents_t<2> const exts2;

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
	//      multi::layout_t<2> const lyt(multi::extents_t<2>({
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
	//      multi::layout_t<2> const lyt(multi::extents_t<2>({
	//          { 0, 10},
	//          {11, 31},
	//      }));
	//      BOOST_TEST( size(lyt) == 10   );
	//      BOOST_TEST( stride(lyt) == 20 );
	//      BOOST_TEST( offset(lyt) == 0 );
	//  }
	//  {  // this is ambiguous in nvcc
	//      multi::layout_t<2> const lyt(multi::extents_t<2>({
	//          { 0, 10},
	//          {11, 31},
	//      }));
	//      BOOST_TEST( size(lyt) == 10   );
	//      BOOST_TEST( stride(lyt) == 20 );
	//      BOOST_TEST( offset(lyt) == 0 );
	//  }
	//  {
	//      multi::layout_t<2> const lyt(multi::extents_t<2>({
	//          {8, 18},
	//          {0, 20},
	//      }));
	//      BOOST_TEST( size(lyt) == 10 );
	//      BOOST_TEST( stride(lyt) == 20 );
	//  }
	//  // {
	//  //  multi::layout_t<3> const lyt(multi::extents_t<3>({
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
	//  multi::extents_t<1> const exts{11};

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
	//  multi::extents_t<2> const exts{3, 5};
	//  auto [eye, jay] = exts.from_linear(7);

	//  BOOST_TEST_REQUIRE( eye == 1 );
	//  BOOST_TEST_REQUIRE( jay == 2 );
	//  //  BOOST_TEST_REQUIRE( std::apply(l, l.extensions().from_linear(9)) == 9 );
	// }

	// BOOST_AUTO_TEST_CASE(extensions_from_linear_2d_std_get) {
	//  multi::extents_t<2> const exts{3, 5};
	//  auto                         eye = std::get<0>(exts.from_linear(7));
	//  auto                         jay = std::get<1>(exts.from_linear(7));
	//  BOOST_TEST_REQUIRE( eye == 1 );
	//  BOOST_TEST_REQUIRE( jay == 2 );
	// }

	// BOOST_AUTO_TEST_CASE(extensions_from_linear_2d_std_get_using) {
	//  multi::extents_t<2> const exts{3, 5};
	//  using std::get;
	//  auto       fl  = exts.from_linear(7L);
	//  auto const eye = get<0>(fl);
	//  auto const jay = get<1>(fl);
	//  BOOST_TEST_REQUIRE( eye == 1 );
	//  BOOST_TEST_REQUIRE( jay == 2 );
	// }

	// BOOST_AUTO_TEST_CASE(extensions_from_linear_2d_get_using) {
	//  multi::extents_t<2> const exts{3, 5};

	//  using multi::detail::get;

	//  auto eye = get<0>(exts.from_linear(7));
	//  auto jay = get<1>(exts.from_linear(7));
	//  BOOST_TEST_REQUIRE( eye == 1 );
	//  BOOST_TEST_REQUIRE( jay == 2 );
	// }

	// BOOST_AUTO_TEST_CASE(extensions_from_linear_2d) {
	//  multi::extents_t<2> const exts{3, 5};

	//  auto ij = exts.from_linear(7);

	//  using multi::detail::get;

	//  BOOST_TEST_REQUIRE( get<0>(ij) == 1 );
	//  BOOST_TEST_REQUIRE( get<1>(ij) == 2 );

	//  multi::layout_t<2> const lyt{exts};
	//  BOOST_TEST_REQUIRE( lyt[get<0>(ij)][get<1>(ij)] == 7 );
	// }

	// BOOST_AUTO_TEST_CASE(extensions_from_linear_3d_std_get) {
	//  multi::extents_t<3> const exts{11, 13, 17};

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
	//  multi::extents_t<3> const exts{11, 13, 17};

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
	//  multi::extents_t<3> const exts{11, 13, 17};

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
	//  multi::extents_t<3> const exts{11, 13, 17};

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
	//      multi::extents_t<1> const exts(10);
	//      BOOST_TEST_REQUIRE(std::get<0>(exts[0]) == 0);
	//      BOOST_TEST_REQUIRE(std::get<0>(exts[1]) == 1);
	//  }
	//  {
	//      multi::extents_t<1> const exts(multi::iextension{0, 10});
	//      BOOST_TEST_REQUIRE(std::get<0>(exts[0]) == 0);
	//      BOOST_TEST_REQUIRE(std::get<0>(exts[1]) == 1);
	//  }
	// }

	// BOOST_AUTO_TEST_CASE(extensionS_2D_iteration) {
	//  {
	//      multi::extents_t<2> exts({3, 5});
	//      BOOST_TEST_REQUIRE(std::get<0>(exts[0]) == 0);
	//      BOOST_TEST_REQUIRE(std::get<0>(exts[1]) == 1);
	//  }
	//  {
	//      multi::extents_t<2> exts({multi::iextension{0, 3}, multi::iextension{0, 5}});
	//      BOOST_TEST_REQUIRE(std::get<0>(exts[0]) == 0);
	//      BOOST_TEST_REQUIRE(std::get<0>(exts[1]) == 1);
	//  }
	// }

	// BOOST_AUTO_TEST_CASE(layout_1D_iteration) {
	//  multi::layout_t<1> const lyt{multi::extents_t<1>(10)};
	//  BOOST_TEST( lyt[0] == 0 );
	//  BOOST_TEST( lyt[1] == 1 );
	//  BOOST_TEST( lyt[2] == 2 );

	//  //  BOOST_TEST_REQUIRE(std::get<0>(exts[0]) == 0);
	//  //  BOOST_TEST_REQUIRE(std::get<0>(exts[1]) == 1);
	// }

	// BOOST_AUTO_TEST_CASE(layout_2D_iteration) {
	//  multi::layout_t<2> const lyt{multi::extents_t<2>({5, 3})};
	//  BOOST_TEST( lyt[0][0] == 0 );
	//  BOOST_TEST( lyt[0][1] == 1 );
	//  BOOST_TEST( lyt[0][2] == 2 );

	//  BOOST_TEST( lyt[1][0] == 3 );
	//  BOOST_TEST( lyt[1][1] == 4 );
	//  BOOST_TEST( lyt[1][2] == 5 );

	//  //  BOOST_TEST_REQUIRE(std::get<0>(exts[0]) == 0);
	//  //  BOOST_TEST_REQUIRE(std::get<0>(exts[1]) == 1);
	// }

	{
		multi::array<int, 2> arr({6, 8});
		BOOST_TEST( arr.size() == 6 );

		auto&& arrp2 = arr.partitioned(2);
		BOOST_TEST( arrp2.num_elements() == arr.num_elements() );
		BOOST_TEST( arrp2.size() == 2 );
	}
	{
		multi::layout_t<2> const lyt(multi::extents_t<2>{
			{3,  9},
			{0, 15}
		});

		BOOST_TEST( lyt.size() == 6 );

		BOOST_TEST( lyt.extent().front() == 3 );
		BOOST_TEST( lyt.extent().back() == 8 );

		auto sorted_lyt = lyt.sort();

		BOOST_TEST( sorted_lyt == lyt );

		auto lyt_transpose = lyt.transpose();

		auto sorted_lyt2 = lyt_transpose.transpose();

		BOOST_TEST( sorted_lyt2 == sorted_lyt );
	}
	{
		multi::extent_t<int> const ext(5);

		BOOST_TEST( *ext.begin() == 0 );
		BOOST_TEST( *(ext.end() - 1) == 4 );

#if !defined(__clang_major__) || (__clang_major__ > 16)
#if defined(__cpp_lib_ranges) && (__cpp_lib_ranges >= 201911L) && !defined(_MSC_VER)
		BOOST_TEST( *std::ranges::begin(ext) == 0 );
		BOOST_TEST( *(std::ranges::end(ext)-1) == 4 );

		BOOST_TEST( ext[0] == 0 );
		BOOST_TEST( ext[1] == 1 );
		BOOST_TEST( ext[4] == 4 );

		static_assert(std::ranges::range<boost::multi::extent_t<int, int>>);
		static_assert(std::ranges::range<boost::multi::extent_t<int, int> const>);

		// std::ranges::ref_view<const boost::multi::extension_t<int, int>>
		auto rext = ext | std::ranges::views::reverse;

		BOOST_TEST( rext[0] == 4 );
		BOOST_TEST( rext[1] == 3 );
		BOOST_TEST( rext[4] == 0 );
#endif
#endif
	}

	return boost::report_errors();
}
