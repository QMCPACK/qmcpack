// Copyright 2019-2023 Alfredo A. Correa

#include <boost/test/unit_test.hpp>

#include <multi/array.hpp>
#include <multi/detail/static_allocator.hpp>

#include <vector>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(empty_stride) {
	multi::array<double, 2> ma;
	BOOST_REQUIRE(ma.size() == 0);
	BOOST_REQUIRE(ma.stride() != 0);
	BOOST_REQUIRE(size(ma) == 0);

	multi::array<double, 2> ma0({0, 0}, 0.0);
	BOOST_REQUIRE(ma0.size() == 0);
	BOOST_REQUIRE(ma0.stride() != 0);
	BOOST_REQUIRE(size(ma0) == 0);
}

BOOST_AUTO_TEST_CASE(std_vector_of_arrays) {
	std::vector<multi::array<double, 2>> va;
	std::transform(
		begin(multi::iextension(3)), end(multi::iextension(3)),
		std::back_inserter(va),
		[](auto idx){return multi::array<double, 2>({idx, idx}, static_cast<double>(idx));}
	);

	BOOST_REQUIRE( size(va[0]) == 0 );
	BOOST_REQUIRE( size(va[1]) == 1 );
	BOOST_REQUIRE( size(va[2]) == 2 );
	BOOST_REQUIRE( va[1] [0][0] == 1 );
	BOOST_REQUIRE( va[2] [0][0] == 2 );

	std::vector<multi::array<double, 2>> const wa = {  // testing std::vector of multi:array NOLINT(fuchsia-default-arguments-calls,-warnings-as-errors)
		multi::array<double, 2>({0, 0}, 0.0),
		multi::array<double, 2>({1, 1}, 1.0),
		multi::array<double, 2>({2, 2}, 2.0),
	};
	BOOST_REQUIRE( size(va) == size(wa) );
	BOOST_REQUIRE( va == wa );

	std::vector<multi::array<double, 2>> ua(3, std::allocator<multi::array<double, 2>>{});
	auto iex = multi::iextension(static_cast<multi::size_type>(ua.size()));
	std::transform(
		begin(iex), end(iex),
		begin(ua),
		[](auto idx) {return multi::array<double, 2>({idx, idx}, static_cast<double>(idx));}
	);
	BOOST_REQUIRE( ua == va );
}

// TODO(correaa) make this code work with nvcc compiler (non device function called from device host through adl uninitialized_fill)
#if !(defined(__NVCC__) || defined(__HIP_PLATFORM_NVIDIA__) || defined(__HIP_PLATFORM_AMD__) || defined(__HIPCC__))
BOOST_AUTO_TEST_CASE(array1d_of_arrays2d) {
	multi::array<multi::array<double, 2>, 1> arr(multi::extensions_t<1>(multi::iextension{10}), multi::array<double, 2>{});
	BOOST_REQUIRE( size(arr) == 10 );

	std::transform(
		begin(extension(arr)), end(extension(arr)), begin(arr),
		[](auto idx) {return multi::array<double, 2>({idx, idx}, static_cast<double>(idx));}
	);

	BOOST_REQUIRE( size(arr[0]) == 0 );
	BOOST_REQUIRE( size(arr[1]) == 1 );
	BOOST_REQUIRE( size(arr[8]) == 8 );
	BOOST_REQUIRE( arr[8][4][4] == 8.0 );
}

BOOST_AUTO_TEST_CASE(array_3d_of_array_2d)  {
	multi::array<multi::array<double, 3>, 2> AA({10, 20}, multi::array<double, 3>{});
	std::transform(extension(AA).begin(), extension(AA).end(), AA.begin(), AA.begin(), [](auto idx, auto&& row) -> decltype(row) {
		std::transform(extension(row).begin(), extension(row).end(), row.begin(), [idx](auto jdx) {
			return multi::array<double, 3>({idx + jdx, idx + jdx, idx + jdx}, 99.0);
		});
		return std::forward<decltype(row)>(row);
	});

	BOOST_REQUIRE( size(AA[9][19]) == 9 + 19 );
	BOOST_REQUIRE( AA[9][19][1][1][1] == 99.0 );
}

BOOST_AUTO_TEST_CASE(array_3d_of_array_2d_no_init)  {
	multi::array<multi::array<double, 3>, 2> AA({10, 20});
	std::transform(extension(AA).begin(), extension(AA).end(), AA.begin(), AA.begin(), [](auto idx, auto&& row) -> decltype(row) {
		std::transform(extension(row).begin(), extension(row).end(), row.begin(), [idx](auto jdx) {
			return multi::array<double, 3>({idx + jdx, idx + jdx, idx + jdx}, 99.0);
		});
		return std::forward<decltype(row)>(row);
	});

	BOOST_REQUIRE( size(AA[9][19]) == 9 + 19 );
	BOOST_REQUIRE( AA[9][19][1][1][1] == 99. );
}
#endif

BOOST_AUTO_TEST_CASE(const_elements) {
	auto ptr = std::make_unique<double const>(2.0);
//  *ptr = 3.0;  // ok, can't assign
	BOOST_REQUIRE( *ptr == 2.0 );

//  multi::array<double const, 2, std::allocator<double>> arr({10, 10}, 99.0);
//  BOOST_REQUIRE( arr[1][2] == 99.0 );
}

#if defined(__cpp_lib_memory_resource) and (__cpp_lib_memory_resource >= 201603)
BOOST_AUTO_TEST_CASE(pmr) {
	std::array<char, 13> buffer = {{'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C'}};
	std::pmr::monotonic_buffer_resource pool{std::data(buffer), std::size(buffer)};

	multi::array<char, 2, std::pmr::polymorphic_allocator<char>> Aarr({2, 2}, 'x', &pool);
	Aarr[0][0] = 'x'; Aarr[0][1] = 'y';
	Aarr[1][0] = 'z'; Aarr[1][1] = '&';

	multi::array<char, 2, std::pmr::polymorphic_allocator<char>> Barr({3, 2}, 'o', &pool);

	BOOST_REQUIRE(( buffer != std::array<char, 13>{{'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C'}} ));

#if defined(__GLIBCXX__)
	BOOST_REQUIRE(( buffer == std::array<char, 13>{{'x', 'y', 'z', '&', 'o', 'o', 'o', 'o', 'o', 'o', 'A', 'B', 'C'}} ));
#endif
#if defined(_LIBCPP_VERSION)
	BOOST_REQUIRE(( buffer == std::array<char, 13>{{'0', '1', '2', 'o', 'o', 'o', 'o', 'o', 'o', 'x', 'y', 'z', '&'}} ));
#endif

	BOOST_REQUIRE(Aarr[0][0] == 'x');
	BOOST_REQUIRE(Barr[0][0] == 'o');
}

BOOST_AUTO_TEST_CASE(pmr2) {
	std::array<char, 13> buffer = {{'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X'}};
	std::pmr::monotonic_buffer_resource pool{std::data(buffer), std::size(buffer)};

	multi::pmr::array<char, 2> Aarr({2, 2}, 'a', &pool);
	multi::pmr::array<char, 2> Barr({3, 2}, 'b', &pool);

#if defined(__GLIBCXX__)
	BOOST_REQUIRE(( buffer == std::array<char, 13>{{'a', 'a', 'a', 'a', 'b', 'b', 'b', 'b', 'b', 'b', 'X', 'X', 'X'}} ));
#endif
#if defined(_LIBCPP_VERSION)
	BOOST_REQUIRE(( buffer == std::array<char, 13>{{'X', 'X', 'X', 'b', 'b', 'b', 'b', 'b', 'b', 'a', 'a', 'a', 'a'}} ));
#endif

	BOOST_REQUIRE(Aarr[0][0] == 'a');
	BOOST_REQUIRE(Barr[0][0] == 'b');
}

BOOST_AUTO_TEST_CASE(pmr_double_uninitialized) {
	std::array<double, 12> buffer = {{4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.00, 11.0,  996.0, 997.0, 998.0, 999.0}};
	std::pmr::monotonic_buffer_resource pool{static_cast<void*>(std::data(buffer)), 12*sizeof(double)};

	multi::pmr::array<double, 2> Aarr({2, 2}, &pool);

	BOOST_TEST( buffer[0] == 4.0 );
	BOOST_TEST( buffer[1] == 5.0 );

#if defined(__GLIBCXX__)
	BOOST_REQUIRE(Aarr[0][0] == 4.0);
#endif
#if defined(_LIBCPP_VERSION)
	BOOST_REQUIRE(Aarr[0][0] == 996.0);
#endif
}

BOOST_AUTO_TEST_CASE(static_allocator) {
	multi::detail::static_allocator<double, 32> sa{};
	auto* pp = sa.allocate(10);
	new(std::next(pp, 8)) double{4.2};
	BOOST_REQUIRE( *std::next(pp, 8) == 4.2 );
	// (pp + 8)->~double();
	sa.deallocate(pp, 10);
}

#if defined( __cpp_constexpr) &&  (__cpp_constexpr  > 202306L)
constexpr auto f() {
    std::vector<int> v = {1, 2, 3};
    return v.size();
}

BOOST_AUTO_TEST_CASE(constexpr_allocator_vector) {
	static_assert(f() == 3);
	BOOST_REQUIRE( f() == 3 );
}

constexpr auto g() {
	multi::array<int, 2> arr = {{4, 5, 6}, {1, 2, 3}, {7, 8, 9}};
	std::sort(arr.begin(), arr.end());
	for(auto it = arr.diagonal().begin(); it != arr.diagonal().end(); ++it) {
		*it += 5;
	}
	auto ret = arr[1][1];
	return ret;
}

BOOST_AUTO_TEST_CASE(constexpr_allocator) {
	constexpr auto gg = g();
	static_assert( gg == 10 );
	BOOST_REQUIRE( gg == 10 );
}
#endif

BOOST_AUTO_TEST_CASE(static_allocator_on_vector_double) {
	std::vector<double, multi::detail::static_allocator<double, 32>> vv(10, 4.2);  // NOLINT(fuchsia-default-arguments-calls)
	BOOST_REQUIRE( vv[3] == 4.2 );

	auto ww = vv;
	BOOST_REQUIRE( ww[3] == 4.2 );

	ww[3] = 5.1;
	BOOST_REQUIRE( ww[3] == 5.1 );
	BOOST_REQUIRE( vv[3] == 4.2 );

	auto xx = std::move(ww);
	BOOST_REQUIRE( ww.empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved)
	BOOST_REQUIRE( vv[3] == 4.2 );
	BOOST_REQUIRE( xx[3] == 5.1 );

	// swap(xx, vv);
	// BOOST_REQUIRE( vv[3] == 5.1 );
	// BOOST_REQUIRE( xx[3] == 4.2 );

	{
		std::vector< std::vector<double, multi::detail::static_allocator<double, 32>> > const VV = {vv, xx, vv};  // NOLINT(fuchsia-default-arguments-calls)
		BOOST_REQUIRE( VV.size() == 3 );
		// swap(VV[0], VV[1]);
		// std::sort(VV.begin(), VV.end());
		// BOOST_REQUIRE( std::is_sorted(VV.begin(), VV.end()) );
		// VV.resize(10, xx);
		// std::sort(VV.begin(), VV.end());
		// BOOST_REQUIRE( std::is_sorted(VV.begin(), VV.end()) );
	}
}

BOOST_AUTO_TEST_CASE(static_allocator_on_vector_string) {
	std::string const cat = "catcatcatcatcatcatcatcatcatcatcatcatcatcatcatcatcatcatcatcatcatcatcatcat";  // NOLINT(fuchsia-default-arguments-calls)
	std::string const dog = "dogdogdogdogdogdogdogdogdogdogdogdogdogdogdogdogdogdogdogdogdogdogdogdog";  // NOLINT(fuchsia-default-arguments-calls)

	std::vector<std::string, multi::detail::static_allocator<std::string, 32>> vv(10, cat);  // NOLINT(fuchsia-default-arguments-calls)
	BOOST_REQUIRE( vv[3] == cat );

	auto ww = vv;
	BOOST_REQUIRE( ww[3] == cat );

	ww[3] = dog;
	BOOST_REQUIRE( ww[3] == dog );
	BOOST_REQUIRE( vv[3] == cat );

	auto xx = std::move(ww);
	BOOST_REQUIRE( vv[3] == cat );
	BOOST_REQUIRE( xx[3] == dog );
	BOOST_REQUIRE( ww.empty() );  // NOLINT(bugprone-use-after-move,hicpp-invalid-access-moved)

	// vv.resize(15);

	// swap(xx, vv);
	// BOOST_REQUIRE( vv[3] == dog );
	// BOOST_REQUIRE( xx[3] == cat );

	{
		std::vector< std::vector<std::string, multi::detail::static_allocator<std::string, 32>> > const VV = {vv, xx, vv};  // NOLINT(fuchsia-default-arguments-calls)
		BOOST_REQUIRE( VV.size() == 3 );
		// swap(VV[0], VV[1]);
		// std::sort(VV.begin(), VV.end());
		// BOOST_REQUIRE( std::is_sorted(VV.begin(), VV.end()) );
		// VV.resize(10, xx);
		// std::sort(VV.begin(), VV.end());
		// BOOST_REQUIRE( std::is_sorted(VV.begin(), VV.end()) );
	}
}

template<class T, multi::dimensionality_type D, std::size_t Capacity = 4UL*4UL>
using small_array = multi::static_array<T, D, multi::detail::static_allocator<T, Capacity>>;
// https://godbolt.org/z/d8ozWahna

BOOST_AUTO_TEST_CASE(small_array_double) {
	small_array<double, 2, 4UL*4UL> vv({4, 4}, 4.2);

	BOOST_REQUIRE( vv[3][3] == 4.2 );

	auto ww = vv;
	BOOST_REQUIRE( ww[3][3] == 4.2 );
	BOOST_REQUIRE( ww.base() != vv.base() );
	auto* wwb = ww.base();
	auto* vvb = vv.base();

	ww[3][3] = 5.1;
	BOOST_REQUIRE( ww[3][3] == 5.1 );
	BOOST_REQUIRE( vv[3][3] == 4.2 );

	swap(ww, vv);
	BOOST_REQUIRE( vv[3][3] == 5.1 );
	BOOST_REQUIRE( ww[3][3] == 4.2 );

	BOOST_REQUIRE( ww.base() == wwb );
	BOOST_REQUIRE( vv.base() == vvb );

	auto xx = std::move(ww);
	BOOST_REQUIRE( vv[3][3] == 5.1 );
	BOOST_REQUIRE( xx[3][3] == 4.2 );
	// BOOST_REQUIRE( ww[3][3] == 4.2 );
	BOOST_REQUIRE( xx.base() != vv.base() );
	// BOOST_REQUIRE( ww.empty() );

	small_array<double, 2, 4UL*4UL> yy({4, 4});
	yy = vv;
	BOOST_REQUIRE( yy == vv );
	yy = std::move(vv);
	BOOST_REQUIRE( vv.size() == 4 );  // NOLINT(clang-analyzer-cplusplus.Move,bugprone-use-after-move,hicpp-invalid-access-moved)

	{
		std::vector< small_array<double, 2, 4UL*4UL> > VV = {vv, xx, vv};  // NOLINT(fuchsia-default-arguments-calls)
		BOOST_REQUIRE( VV.size() == 3 );
		swap(VV[0], VV[1]);
		std::sort(VV.begin(), VV.end());
		BOOST_REQUIRE( std::is_sorted(VV.begin(), VV.end()) );
		VV.resize(10, xx);
		std::sort(VV.begin(), VV.end());
		BOOST_REQUIRE( std::is_sorted(VV.begin(), VV.end()) );
	}
}

BOOST_AUTO_TEST_CASE(props_of_static_allocator) {
	{
        std::vector<double> vv(20, 0.1);  // NOLINT(fuchsia-default-arguments-calls)
        std::vector<double> ww = vv;
        BOOST_REQUIRE( ww == vv );

        ww = vv;
        BOOST_REQUIRE( ww == vv );

        ww = std::move(vv);
        BOOST_REQUIRE( vv.size() == 0 );  // NOLINT(readability-container-size-empty,bugprone-use-after-move,hicpp-invalid-access-moved,clang-analyzer-cplusplus.Move)

        std::vector<double> xx(20, 0.2);  // NOLINT(fuchsia-default-arguments-calls)
		swap( ww, xx );
		BOOST_REQUIRE( ww == std::vector<double>(20, 0.2) );  // NOLINT(fuchsia-default-arguments-calls)
    }
    {
        std::vector<double, multi::detail::static_allocator<double, 32>> vv(20, 0.1);  // NOLINT(fuchsia-default-arguments-calls)
        std::vector<double, multi::detail::static_allocator<double, 32>> ww = vv;
        BOOST_REQUIRE( ww == vv );

        ww = vv;
        BOOST_REQUIRE( ww == vv );

        ww = std::move(vv);
        BOOST_REQUIRE( vv.size() == 0 );  // NOLINT(readability-container-size-empty,bugprone-use-after-move,hicpp-invalid-access-moved,clang-analyzer-cplusplus.Move)

		std::vector<double, multi::detail::static_allocator<double, 32>> xx(20, 0.2);  // NOLINT(fuchsia-default-arguments-calls)
		swap( ww, xx );
		BOOST_REQUIRE(( ww == std::vector<double, multi::detail::static_allocator<double, 32>>(20, 0.2) ));  // NOLINT(fuchsia-default-arguments-calls)
    }
}
#endif
