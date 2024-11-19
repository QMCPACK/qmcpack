// Copyright 2019-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>  // for array, extension_t, static_array

#include <algorithm>  // for fill_n
#include <chrono>     // for high_resolution_clock, operator-  // NOLINT(build/c++11)
#include <cmath>      // for abs  // IWYU pragma: keep
// IWYU pragma: no_include <cstdlib>                          // for abs
#include <cstdint>     // for int64_t
#include <functional>  // for plus  // IWYU pragma: keep
#include <iostream>    // for char_traits, basic_ostream, oper...
#include <iterator>    // for size, data

#if __has_include(<memory_resource>)
	#include <memory_resource>  // for polymorphic_allocator, monotonic...
#endif

#include <numeric>  // for accumulate, transform_reduce

namespace multi = boost::multi;

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	BOOST_AUTO_TEST_CASE(dummy_test) {
	}

#ifdef BOOST_MULTI_HAS_MEMORY_RESOURCE
	BOOST_AUTO_TEST_CASE(pmr_partially_formed) {
		{
	#if defined(__clang__)
		#pragma clang diagnostic push
		#pragma clang diagnostic ignored "-Wunknown-warning-option"
		#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"
	#endif

			// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) use raw memory
			char buffer[] = "0123456789012345678901234567890123456789012345678901234567890123456789";

			std::pmr::monotonic_buffer_resource mbr{std::data(buffer), std::size(buffer)};
			static_assert(std::size(buffer) > 6 * sizeof(double));

			multi::array<double, 2, std::pmr::polymorphic_allocator<double>> const arr({2, 3}, &mbr);
			BOOST_TEST( buffer[ 0] == '0' );  // buffer is intact when initializing without value
			BOOST_TEST( buffer[13] == '3' );

	#if defined(__clang__)
		#pragma clang diagnostic pop
	#endif

			BOOST_TEST( arr.num_elements() == 2*3L );
			//  BOOST_TEST( arr[0][0] != 0.0 );
			//  BOOST_TEST( arr[1][2] != 0.0 );
		}
		{
			// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) use raw memory
			char buffer[] = "0123456789012345678901234567890123456789012345678901234567890123456789";

			std::pmr::monotonic_buffer_resource mbr(std::data(buffer), std::size(buffer));
			static_assert(std::size(buffer) > 6 * sizeof(double));

			multi::array<double, 2, std::pmr::polymorphic_allocator<double>> A({2, 3}, 0.0, &mbr);  // NOLINT(readability-identifier-length)

			BOOST_TEST( std::abs( A[0][0] - 0.0 ) < 1E-6 );
			BOOST_TEST( std::abs( A[1][2] - 0.0 ) < 1E-6 );
		}
		{
			// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) use raw memory
			char buffer[] = "0123456789012345678901234567890123456789012345678901234567890123456789";

			std::pmr::monotonic_buffer_resource mbr(std::data(buffer), std::size(buffer));
			static_assert(std::size(buffer) > 6 * sizeof(double));

			multi::array<double, 2, std::pmr::polymorphic_allocator<double>> arr({2, 3}, {}, &mbr);

			BOOST_TEST( std::abs( arr[0][0] - double{} ) < 1E-6 );
			BOOST_TEST( std::abs( arr[1][2] - double{} ) < 1E-6 );
		}
		{
			// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) use raw memory
			char buffer[] = "0123456789012345678901234567890123456789012345678901234567890123456789";

			std::pmr::monotonic_buffer_resource mbr{std::data(buffer), std::size(buffer)};
			static_assert(std::size(buffer) > 6 * sizeof(double));

			multi::array<double, 2, std::pmr::polymorphic_allocator<double>> arr({2, 3}, 666.0, &mbr);

			BOOST_TEST( std::abs( arr[0][0] - 666.0 ) < 1E-6 );
			BOOST_TEST( std::abs( arr[1][2] - 666.0 ) < 1E-6 );
		}
	}

	#ifndef _MSC_VER  // problems with MSVC 14.3 c++17
	BOOST_AUTO_TEST_CASE(pmr_benchmark) {
		//  auto* resp = std::pmr::unsynchronized_pool_resource(std::pmr::get_default_resource());
		auto* resp = std::pmr::get_default_resource();

		auto count = 50;

		auto start_time = std::chrono::high_resolution_clock::now();

		multi::extension_t const exts{0, count};

		auto acc = std::transform_reduce(
			exts.begin(), exts.end(), int64_t{0},
			std::plus<>{},
			[&resp](auto idx) {
				multi::array<int64_t, 2, std::pmr::polymorphic_allocator<int64_t>> arr(
					multi::extensions_t<2>{1000 - idx % 10, 1000 + idx % 10},  // MSVC needs multi::extensions_t<2>
					resp
				);
				std::fill_n(arr.data_elements(), arr.num_elements(), 1);

				auto const be = arr.elements().begin();
				auto const en = arr.elements().end();
				return std::accumulate(be, en, int64_t{0}, std::plus<int64_t>{});
			}
		);

		auto time = std::chrono::high_resolution_clock::now() - start_time;
		std::cout << time.count() / count << "          " << acc << '\n';
	}
	#endif

#endif
	return boost::report_errors();
}
