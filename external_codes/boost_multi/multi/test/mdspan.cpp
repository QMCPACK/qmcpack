// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <boost/core/lightweight_test.hpp>

#if (__cplusplus >= 202302L || (defined(_MSVC_LANG) && _MSVC_LANG > 202002L))
#if __has_include(<mdspan>)
#include <mdspan>
#endif
#endif

namespace multi = boost::multi;

#if defined(__cpp_lib_mdspan) && (__cpp_lib_mdspan >= 202207L)
template<
	class MultiArray,
	typename T                   = std::remove_reference_t<typename std::decay_t<MultiArray>::element_cref>,
	multi::dimensionality_type D = std::decay_t<MultiArray>::dimensionality>
auto to_strided_mdspan(MultiArray&& arr) {
	using std::apply;
	auto shape = apply(
		[](auto... sizes) { return std::dextents<std::size_t, D>{static_cast<std::size_t>(sizes)...}; },
		arr.sizes()
	);

	auto strides = apply(
		[](auto... strds) { return std::array<std::size_t, D>{static_cast<std::size_t>(strds)...}; },
		arr.strides()
	);

	return std::mdspan<T, std::dextents<std::size_t, D>, std::layout_stride>{
		arr.base(), std::layout_stride::mapping{shape, strides}
	};
}

auto fun_const(std::mdspan<int const, std::dextents<std::size_t, 2>, std::layout_stride> mds) -> int const& {
	return mds[0, 0];
}

auto fun(std::mdspan<int, std::dextents<std::size_t, 2>, std::layout_stride> mds) -> int& {
	return mds[0, 0];
}

#endif

auto main() -> int {
	{
		multi::array<int, 2> const carr = {
			{ 1,  2,  3,  4},
			{ 5,  6,  7,  8},
			{ 9, 10, 11, 12},
			{13, 14, 15, 16}
		};

		BOOST_TEST( carr[1][1] == 6 );

#if defined(__cpp_lib_mdspan) && (__cpp_lib_mdspan >= 202207L)
		std::mdspan<int const, std::dextents<std::size_t, 2>, std::layout_stride> mdspan_arr_const = carr;
		BOOST_TEST(( &mdspan_arr_const[0, 0] == &carr[0][0] ));
		fun_const(mdspan_arr_const);
		fun_const(carr);

		auto const& center = carr({1, 3}, {1, 3});
		BOOST_TEST( &center[0][0] == &carr[1][1] );

		std::mdspan<int const, std::dextents<std::size_t, 2>, std::layout_stride> mdspan_center = center;
		BOOST_TEST(( &mdspan_center[0, 0] == &center[0][0] ));
		fun_const(mdspan_center);
		fun_const(center);

		std::mdspan<int const, std::dextents<std::size_t, 2>, std::layout_stride> mdspan_center_rv = carr({1, 3}, {1, 3});
		BOOST_TEST(( &mdspan_center_rv[0, 0] == &carr({1, 3}, {1, 3})[0][0] ));
		fun_const(mdspan_center_rv);
		fun_const(carr({1, 3}, {1, 3}));
#endif
	}
	{
		multi::array<int, 2> arr = {
			{ 1,  2,  3,  4},
			{ 5,  6,  7,  8},
			{ 9, 10, 11, 12},
			{13, 14, 15, 16}
		};

		BOOST_TEST( arr[1][1] == 6 );

#if defined(__cpp_lib_mdspan) && (__cpp_lib_mdspan >= 202207L)
		std::mdspan<int const, std::dextents<std::size_t, 2>, std::layout_stride> mdspan_arr_const = arr;
		BOOST_TEST(( &mdspan_arr_const[0, 0] == &arr[0][0] ));
		fun_const(mdspan_arr_const);
		fun_const(arr);

		std::mdspan<int, std::dextents<std::size_t, 2>, std::layout_stride> mdspan_arr = arr;
		BOOST_TEST(( &mdspan_arr[0, 0] == &arr[0][0] ));
		fun_const(mdspan_arr);
		fun_const(arr);

		fun(mdspan_arr);
		fun(arr);

		auto&& center = arr({1, 3}, {1, 3});
		BOOST_TEST( &center[0][0] == &arr[1][1] );

		std::mdspan<int const, std::dextents<std::size_t, 2>, std::layout_stride> mdspan_center_const = center;
		BOOST_TEST(( &mdspan_center_const[0, 0] == &center[0][0] ));
		fun_const(mdspan_center_const);
		fun_const(center);

		std::mdspan<int, std::dextents<std::size_t, 2>, std::layout_stride> mdspan_center = center;
		BOOST_TEST(( &mdspan_center[0, 0] == &center[0][0] ));
		fun_const(mdspan_center);
		fun_const(center);

		fun(mdspan_center);
		fun(center);

		std::mdspan<int, std::dextents<std::size_t, 2>, std::layout_stride> mdspan_center_rv = arr({1, 3}, {1, 3});
		BOOST_TEST(( &mdspan_center_rv[0, 0] == &arr({1, 3}, {1, 3})[0][0] ));
		fun_const(mdspan_center_rv);
		fun_const(arr({1, 3}, {1, 3}));

		fun(mdspan_center_rv);
		fun(arr({1, 3}, {1, 3}));

		std::mdspan<int const, std::dextents<std::size_t, 2>, std::layout_stride> mdspan_center_crv = arr({1, 3}, {1, 3});
		BOOST_TEST(( &mdspan_center_crv[0, 0] == &arr({1, 3}, {1, 3})[0][0] ));
		fun_const(mdspan_center_crv);
		fun_const(arr({1, 3}, {1, 3}));

		// fun(mdspan_center_crv);  // doesn't compile, good
		fun(arr({1, 3}, {1, 3}));
#endif
	}
	// #if defined(__cpp_lib_mdspan) && (__cpp_lib_mdspan >= 202207L)
	// 	std::mdspan<int const, std::dextents<std::size_t, 2>, std::layout_stride> mdspan_arr_const = arr;
	// 	BOOST_TEST(( &mdspan_arr_const[0, 0] == &arr[0][0] ));
	// 	BOOST_TEST( &fun(arr) == &arr[0][0] );

	// 	// std::mdspan<int const, std::dextents<std::size_t, 1>, std::layout_stride> mds1D = center[0];
	// 	// BOOST_TEST(( &mds1D[0] == &center[0][0] ));

	// 	std::mdspan<int, std::dextents<std::size_t, 2>, std::layout_stride> mdspan_arr = arr;
	// 	BOOST_TEST(( &mdspan_arr[0, 0] == &arr[0][0] ));
	// #endif

	// 	auto const& center = arr({1, 3}, {1, 3});
	// 	BOOST_TEST( &center[0][0] == &arr[1][1] );

	// #if defined(__cpp_lib_mdspan) && (__cpp_lib_mdspan >= 202207L)
	// 	std::mdspan<int const, std::dextents<std::size_t, 2>, std::layout_stride> mds = center;
	// 	BOOST_TEST(( &mds[0, 0] == &center[0][0] ));

	// 	BOOST_TEST( &fun_const(center) == &center[0][0] );

	// 	std::mdspan<int const, std::dextents<std::size_t, 1>, std::layout_stride> mds1D = center[0];
	// 	BOOST_TEST(( &mds1D[0] == &center[0][0] ));
	// #endif

	return boost::report_errors();
}
