// Copyright 2018-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for array, subarray, static_array

#include <boost/core/lightweight_test.hpp>

#include <array>    // for get, array
#include <cmath>    // for abs  // IWYU pragma: keep
#include <complex>  // for complex, real, operator==, imag
// IWYU pragma: no_include <cstdlib>                          // for abs
#include <iterator>  // for size, begin, end
// IWYU pragma: no_include <memory>       // for allocator
#include <numeric>      // for iota
#include <type_traits>  // for is_same_v
#include <utility>      // for pair

namespace multi = boost::multi;

template<class T> struct complex_dummy {
	T real;
	T imag;
};

#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	// simple class reinterpret
	{
		struct int_class {int val_;};

		multi::array<int_class, 1> arr1d({3}, int_class{});
		arr1d[0].val_ = 5;
		BOOST_TEST( arr1d.reinterpret_array_cast<int>()[0] == 5 );

		multi::array<int_class, 2> arr2d({3, 3}, int_class{});
		arr2d[0][0].val_ = 5;
		BOOST_TEST( arr2d.reinterpret_array_cast<int>()[0][0] == 5 );
	}

	BOOST_AUTO_TEST_CASE(multi_reinterpret_array_cast_struct_to_dimension) {
		struct vec3 {
			double x;
			double y;
			double z;
		};
		multi::array<vec3, 1> arr(multi::extensions_t<1>{multi::iextension{100}});
		arr[8] = {1.0, 2.0, 3.0};
		BOOST_TEST( std::abs( arr[8].y - 2.0 ) < 1E-6 );

#ifndef _MSC_VER  // problems with MSVC 14.3 c++17
		BOOST_TEST( std::abs( arr.reinterpret_array_cast<double>(3)[8][1] - arr[8].y ) < 1E-6 );

		multi::array<double, 2> A2D{arr.reinterpret_array_cast<double>(3)};

		BOOST_TEST( decltype(A2D)::dimensionality == decltype(arr)::dimensionality + 1 );
		BOOST_TEST( dimensionality(A2D) == dimensionality(arr) + 1 );

		BOOST_TEST( size(A2D) == size(arr) );
		BOOST_TEST( std::abs( A2D[8][1] - arr[8].y ) < 1E-6);
		BOOST_TEST( &A2D[8][1] != &arr[8].y );

		BOOST_TEST( & arr[8].x == & arr.reinterpret_array_cast<double>(3)[8][0] );
		BOOST_TEST( & arr[8].y == & arr.reinterpret_array_cast<double>(3)[8][1] );
		BOOST_TEST( & arr[8].z == & arr.reinterpret_array_cast<double>(3)[8][2] );
#endif
	}

	BOOST_AUTO_TEST_CASE(multi_lower_dimension) {
		struct vec3 {
			double x;
			double y;
			double z;

			// [[maybe_unused]] auto operator==(vec3 const& other) const -> bool { return x == other.x && y == other.y && z == other.z; }
		};

		multi::array<double, 2> arr = {
			{0.0, 0.1, 0.2},
			{1.0, 1.1, 1.2},
			{2.0, 2.1, 2.2},
			{3.0, 3.1, 3.2},
		};
		{
			BOOST_TEST( arr.size() == 4 );
			BOOST_TEST( arr.flatted().size() == 12 );
			BOOST_TEST( arr.flatted().strided(3).size() == 4 );
			BOOST_TEST( arr.flatted().strided(3).reinterpret_array_cast<vec3>().size() == 4 );

			auto&& arrvec3 = arr.flatted().strided(3).reinterpret_array_cast<vec3>();

			BOOST_TEST( arr.flatted().size() == arrvec3.size()*3 );
			BOOST_TEST( &arrvec3[2].x == &arr[2][0] );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_lower_dimension_2d) {
		struct vec3 {
			double x;
			double y;
			double z;
		};

		multi::array<double, 2> d2 = {
			{0.0, 0.1, 0.2, 0.0, 0.1, 0.2, 0.0, 0.1, 0.2},
			{1.0, 1.1, 1.2, 1.0, 1.1, 1.2, 1.0, 1.1, 1.2},
			{2.0, 2.1, 2.2, 2.0, 2.1, 2.2, 2.0, 2.1, 2.2},
			{3.0, 3.1, 3.2, 3.0, 3.1, 3.2, 3.0, 3.1, 3.2},
		};

		{
			auto&& d2strided3 = d2.unrotated().strided(3).rotated();
			BOOST_TEST( d2strided3.size() == 4 );
			BOOST_TEST( d2strided3[0].size() == 3 );
			BOOST_TEST( &d2strided3[1][2] == &d2[1][6] );
		}
		{
			auto&& v2view = d2.unrotated().strided(3).rotated().reinterpret_array_cast<vec3>();
			BOOST_TEST( v2view.size() == 4 );
			BOOST_TEST( v2view[0].size() == 3 );
			BOOST_TEST( &v2view[1][2].x == &d2[1][6] );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_lower_dimension_3d) {
		struct vec3 {
			double x;
			double y;
			double z;
		};

		multi::array<double, 3> d3({4, 15, 9}, 0.0);

		{
			auto&& d3strided3 = d3.unrotated().strided(3).rotated();
			BOOST_TEST( d3strided3.size() == 4 );
			BOOST_TEST( d3strided3[0][0].size() == 3 );
			BOOST_TEST( &d3strided3[3][1][2] == &d3[3][1][6] );
		}
		{
			auto&& v3view = d3.unrotated().strided(3).rotated().reinterpret_array_cast<vec3>();
			BOOST_TEST( v3view.size() == 4 );
			BOOST_TEST( v3view[0][0].size() == 3 );
			BOOST_TEST( &v3view[3][1][2].x == &d3[3][1][6] );
			BOOST_TEST( &v3view[3][1][2].y == &d3[3][1][7] );
			BOOST_TEST( &v3view[3][1][2].z == &d3[3][1][8] );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_reinterpret_array_cast_complex_to_real_extra_dimension) {
		using complex = std::complex<double>;
		multi::array<complex, 1> arr(multi::extensions_t<1>{multi::iextension{100}}, complex{1.0, 2.0});
		BOOST_TEST(  size(arr) == 100 );

		{
			complex const arr0 = arr[0];
			BOOST_TEST( std::abs( arr0.real() - 1.0 ) < 1E-6 );
			BOOST_TEST( std::abs( arr0.imag() - 2.0 ) < 1E-6 );
		}

		BOOST_TEST( std::abs( arr[0].real() - 1.0 ) < 1E-6 );
		BOOST_TEST( std::abs( arr[0].imag() - 2.0 ) < 1E-6 );

		BOOST_TEST( std::abs( std::real(arr[0]) - 1.0 ) < 1E-6 );
		BOOST_TEST( std::abs( std::imag(arr[0]) - 2.0 ) < 1E-6 );

		BOOST_TEST( std::abs( real(arr[0]) - 1.0 ) < 1E-6 );
		BOOST_TEST( std::abs( imag(arr[0]) - 2.0 ) < 1E-6 );

		BOOST_TEST( std::abs( arr[0] - complex{1.0, 2.0} ) < 1E-6 );

#ifndef _MSC_VER  // problem with MVSC 14.3 c++17
		multi::array<double, 1> arr2{arr.reinterpret_array_cast<double>()};
		BOOST_TEST( dimensionality(arr2) == dimensionality(arr) );
		BOOST_TEST( std::abs( arr2[0] - 1.0 ) < 1E-6 );
		BOOST_TEST( std::abs( arr2[1] - 1.0 ) < 1E-6 );

		multi::array<double, 2> arr3{arr.reinterpret_array_cast<double>(2)};

		BOOST_TEST(( arr3.sizes() == decltype(arr3.sizes()){100, 2} ));
		BOOST_TEST( std::abs( arr3[5][0] - real(arr[5]) ) < 1E-6 );
		BOOST_TEST( std::abs( arr3[5][1] - imag(arr[5]) ) < 1E-6 );

		arr.reinterpret_array_cast<double>(2)[0][0]   = 99.9;
		arr().reinterpret_array_cast<double>(2)[0][0] = 99.9;

#endif
	}

	BOOST_AUTO_TEST_CASE(multi_reinterpret_array_cast_tuple_as_extra_dimension) {
		using vector3 = std::array<double, 3>;

		vector3 v3d;

#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wunknown-warning-option"
	#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
#endif

		// NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast, cppcoreguidelines-avoid-c-arrays, hicpp-avoid-c-arrays, modernize-avoid-c-arrays): test
		BOOST_TEST( &reinterpret_cast<double(&)[3]>(v3d)[1] == &std::get<1>(v3d) );

#if defined(__clang__)
	#pragma clang diagnostic pop
#endif

#ifndef _MSC_VER  // problem with MVSC 14.3 c++17
		{
			multi::array<vector3, 1> arr(multi::extensions_t<1>{multi::iextension{10}});
			BOOST_TEST( &arr.reinterpret_array_cast<double>(3)[2][1] == &std::get<1>(arr[2]) );
		}
		{
			multi::array<vector3, 2> arr({10, 20});
			BOOST_TEST( &arr.reinterpret_array_cast<double>(3)[5][7][2] == &std::get<2>(arr[5][7]) );
		}
		{
			multi::array<vector3, 2> const arr({4, 5}, vector3{{1.0, 2.0, 3.0}});

			BOOST_TEST( arr.reinterpret_array_cast<double>(3).dimensionality == 3 );
			BOOST_TEST( decltype(arr.reinterpret_array_cast<double>(3))::dimensionality == 3 );
			BOOST_TEST( dimensionality(arr.reinterpret_array_cast<double>(3)) == 3 );

			BOOST_TEST(  arr.reinterpret_array_cast<double>(3).num_elements() == arr.num_elements()*3 );
			BOOST_TEST(  arr.reinterpret_array_cast<double>(3).size() == 4 );
			BOOST_TEST(  arr.reinterpret_array_cast<double>(3)[0].size() == 5 );
			BOOST_TEST(  arr.reinterpret_array_cast<double>(3)[0][0].size() == 3 );
			BOOST_TEST( &arr.reinterpret_array_cast<double>(3)[2][3][0] == &std::get<0>(arr[2][3]) );
			BOOST_TEST( &arr.reinterpret_array_cast<double>(3)[2][3][1] == &std::get<1>(arr[2][3]) );
			BOOST_TEST( &arr.reinterpret_array_cast<double>(3)[2][3][2] == &std::get<2>(arr[2][3]) );

			multi::array<double, 3> const arr2{arr.reinterpret_array_cast<double>(3)};
			BOOST_TEST( std::abs( arr2[2][3][0] - std::get<0>(arr[2][3]) ) < 1E-6 );
			BOOST_TEST( std::abs( arr2[2][3][1] - std::get<1>(arr[2][3]) ) < 1E-6 );
			BOOST_TEST( std::abs( arr2[2][3][2] - std::get<2>(arr[2][3]) ) < 1E-6 );

			auto arr3 = +arr.reinterpret_array_cast<double>(3);
			BOOST_TEST( arr3 == arr2 );
		}
#endif
	}

	BOOST_AUTO_TEST_CASE(multi_reinterpret_array_cast) {
		{
			std::complex<double> cee{1.0, 2.0};

			auto* ptr = reinterpret_cast<complex_dummy<double>*>(&cee);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
			ptr->real = 11.0;
			BOOST_TEST( std::abs( real(cee) - 11.0 ) < 1E-6 );
		}
		{
			multi::array<std::complex<double>, 1> arr(multi::extensions_t<1>{multi::iextension{10}});
			std::iota(begin(arr), end(arr), 1.0);
			BOOST_TEST( std::abs( real(arr[8]) - 9.0 ) < 1E-6);
			auto&& arr2  = arr.reinterpret_array_cast<complex_dummy<double>>();
			arr2[8].real = 1000.0;
			BOOST_TEST( std::abs( real(arr[8]) - 1000.0 ) < 1E-6 );
		}
	}

	BOOST_AUTO_TEST_CASE(multi_reinterpret_array_cast_realcomplex) {
		using complex = std::complex<double>;
		{
			complex cee{1.0, 2.0};
			auto*   conjd_cee = reinterpret_cast<std::array<double, 2>*>(&cee);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
			(*conjd_cee)[0]   = 11.0;
			BOOST_TEST( conjd_cee );
			BOOST_TEST( std::abs( real(cee) - 11.0 ) < 1E-6 );
		}
		{
			complex cee{1.0, 2.0};
			// NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast, cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays): test purposes
			auto* ceePC = reinterpret_cast<double(*)[2]>(&cee);
			(*ceePC)[0] = 11.0;
			BOOST_TEST( ceePC );
			BOOST_TEST( std::abs( real(cee) - 11.0 ) < 1E-6 );
		}
#ifndef _MSC_VER  // problem with MSVC 14.3 c++17
		{
			multi::array<complex, 1> arr(multi::extensions_t<1>{multi::iextension{10}});

			auto&& arr2 = arr.reinterpret_array_cast<double>(2);

			arr2[8][0] = 1000.0;
			arr2[8][1] = 2000.0;

			BOOST_TEST(( arr[8] == std::complex<double>{1000.0, 2000.0} ));
		}
#endif
	}

	BOOST_AUTO_TEST_CASE(multi_reinterpret_array_cast_pair_to_complex) {
		using complex = std::complex<double>;
		using pair    = std::pair<double, double>;
		multi::array<complex, 2> arr({10, 10}, complex{3.0, 4.0});

		multi::array<complex, 2> const& Aconst  = arr;
		auto&&                          A_block = Aconst({0, 5}, {0, 5});

		// const is important below
		auto const& Apair_block = A_block.template reinterpret_array_cast<pair const>();  // cppcheck 1.90 needs `template` to avoid internal bug

		BOOST_TEST( &Apair_block[1][2] == static_cast<void*>(&arr[1][2]) );

#ifndef _MSC_VER  // problems with MSVC 14.3 c++17
		auto&& Adoubles_block = A_block.reinterpret_array_cast<double const>(2);
		BOOST_TEST( &Adoubles_block[1][2][0] == static_cast<void*>(&arr[1][2]) );
#endif
	}

	BOOST_AUTO_TEST_CASE(multi_reinterpret_array_cast_pointer) {
		multi::array<double, 2> arr({10, 10}, 5.0);

		auto&& Aconstcast = arr.reinterpret_array_cast<double, double const*>();
		BOOST_TEST( &arr[0][0] == &Aconstcast[0][0] );
		static_assert(std::is_same_v<decltype(Aconstcast[1][2]), double const&>);
	}

	BOOST_AUTO_TEST_CASE(const_array_cast) {
		multi::array<double, 2> arr({10, 10}, 5.0);  // NOLINT(misc-const-correctness) test const cast

		multi::array<double, 2> const& carr = arr;

		auto&& marr = carr.const_array_cast<double>();

		marr[1][1] = 6.0;

		BOOST_TEST( std::abs( carr[1][1] - 6.0 ) < 1E-6 );
	}

	// test packing 4 doubles
#ifndef _MSC_VER
	{
		using packed_type = std::array<double, 4>;
		multi::array<packed_type, 2> arr_4pc({10, 25});

		BOOST_TEST( arr_4pc[0].size() == 25 );

		arr_4pc[0][0] = packed_type{{1.0, 2.0, 3.0, 4.0}};
		arr_4pc[0][1] = packed_type{{5.0, 6.0, 7.0, 8.0}};

		auto&& arr = arr_4pc.reinterpret_array_cast<double>(4).rotated().flatted().unrotated();

		BOOST_TEST( arr[0].size() == 100 );

		BOOST_TEST( std::abs( arr[0][0] - 1.0 ) < 1E-6 );
		BOOST_TEST( std::abs( arr[0][1] - 2.0 ) < 1E-6 );
		BOOST_TEST( std::abs( arr[0][2] - 3.0 ) < 1E-6 );
		BOOST_TEST( std::abs( arr[0][3] - 4.0 ) < 1E-6 );

		BOOST_TEST( std::abs( arr[0][4] - 5.0 ) < 1E-6 );
		BOOST_TEST( std::abs( arr[0][5] - 6.0 ) < 1E-6 );
		BOOST_TEST( std::abs( arr[0][6] - 7.0 ) < 1E-6 );
		BOOST_TEST( std::abs( arr[0][7] - 8.0 ) < 1E-6 );
	}
#endif

	return boost::report_errors();
}
