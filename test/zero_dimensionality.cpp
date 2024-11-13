// Copyright 2019-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for array_ref, static_array, array_ptr

// IWYU pragma: no_include <algorithm>  // for copy
#include <complex>  // for complex
#include <memory>  // for allocator  // IWYU pragma: keep
// IWYU pragma: no_include <type_traits>  // for remove_reference<>::type
// #include <utility>  // for move
#include <vector>  // for vector, allocator

namespace multi = boost::multi;

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	BOOST_AUTO_TEST_CASE(zero_dimensionality_part1) {
		{
			std::vector<int> v1 = {10, 20, 30};  // NOLINT(fuchsia-default-arguments-calls)

			multi::array_ref<int, 1> m1(v1.data(), multi::extensions_t<1>{multi::iextension{3}});
			BOOST_TEST( size(m1) == 3 );
			BOOST_TEST( &m1[1] == &v1[1] );
			BOOST_TEST( num_elements(m1) == 3 );

			multi::array_ref<int, 0> const m0(v1.data(), {});
			// BOOST_TEST(( &m0 == multi::array_ptr<double, 0>(v1.data(), {}) ));
			BOOST_TEST( m0.data_elements() == v1.data() );
			BOOST_TEST( m0.num_elements() == 1 );

			// m0 = 51;
			// BOOST_TEST( v1[0] == 51 );

			// int const& doub = std::move(m0);
			// BOOST_TEST( doub == 51 );
		}
		{
			// vvv--- TODO(correaa) this might trigger a compiler crash with g++ 7.5 because of operator&() && overloads
			multi::static_array<double, 0> a0 = multi::static_array<double, 0>{45.0};
			BOOST_TEST( num_elements(a0) == 1 );
			BOOST_TEST( a0 == 45.0 );

			a0 = multi::static_array<double, 0>{60.0};
			BOOST_TEST( a0 == 60.0 );
		}
		{
			std::allocator<double> const   alloc;
			multi::static_array<double, 0> a0(45.0, alloc);
			BOOST_TEST( num_elements(a0) == 1 );
			BOOST_TEST( a0 == 45.0 );

			a0 = multi::static_array<double, 0>{60.0};
			BOOST_TEST( a0 == 60.0 );
		}
	}

	BOOST_AUTO_TEST_CASE(zero_dimensionality_part2) {
		{
			multi::array<std::complex<double>, 2> const arr(
#ifdef _MSC_VER  // problem with 14.3 c++17
				multi::extensions_t<2>
#endif
				{1, 2},
				std::allocator<std::complex<double>>{}
			);
			BOOST_TEST( arr.size() == 1 );
		}
		{
			int                      doub = 20;
			multi::array_ref<int, 0> arr(doub);
			int const&               the_doub = static_cast<int&>(arr);
			BOOST_TEST(  the_doub ==  doub );
			BOOST_TEST( &the_doub == &doub );
		}
		{
			int  doub = 20;
			auto dd   = static_cast<int>(multi::array_ref<int, 0>(&doub, {}));

			BOOST_TEST( dd == doub );

			multi::array_ptr<int, 1> const ap1(&doub, multi::extensions_t<1>({0, 1}));
			BOOST_TEST( (*ap1).base() == &doub );
			BOOST_TEST( ap1->base() == &doub );

			// multi::array_ptr<int, 0> const ap0(&doub, {});

			// BOOST_TEST(( ap0 == multi::array_ptr<int, 0>(&doub, {}) ));
			// BOOST_TEST(( ap0 != multi::array_ptr<int, 0>(&dd, {}) ));
			// // BOOST_TEST( ap0->base() == &doub );
			// // BOOST_TEST( (*ap0).base() == &doub );

			// multi::array_ptr<int, 0> const ap0dd{&dd};
			// BOOST_TEST( ap0dd != ap0 );
			// BOOST_TEST( *ap0 == *ap0dd );
			// int d3 = 3141592;
			// BOOST_TEST(( *multi::array_ptr<int, 0>(&d3, {}) == 3141592 ));
		}
	}

	return boost::report_errors();
}
