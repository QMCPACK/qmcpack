// Copyright 2019-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <numeric>

// Suppress warnings from boost.test
#if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wold-style-cast"
#  pragma clang diagnostic ignored "-Wundef"
#  pragma clang diagnostic ignored "-Wconversion"
#  pragma clang diagnostic ignored "-Wsign-conversion"
#  pragma clang diagnostic ignored "-Wfloat-equal"
#elif defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wold-style-cast"
#  pragma GCC diagnostic ignored "-Wundef"
#  pragma GCC diagnostic ignored "-Wconversion"
#  pragma GCC diagnostic ignored "-Wsign-conversion"
#  pragma GCC diagnostic ignored "-Wfloat-equal"
#elif defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4244)  // narrowing conversion
#endif

#ifndef BOOST_TEST_MODULE
#  define BOOST_TEST_MAIN
#endif

#include <boost/test/unit_test.hpp>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(pmr_dummy) {
}

#ifdef BOOST_MULTI_HAS_MEMORY_RESOURCE
BOOST_AUTO_TEST_CASE(pmr_partially_formed) {
	{
		char buffer[] = "0123456789012345678901234567890123456789012345678901234567890123456789";  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) use raw memory

		std::pmr::monotonic_buffer_resource mbr{std::data(buffer), std::size(buffer)};
		static_assert( std::size(buffer) > 6*sizeof(double) );

		multi::array<double, 2, std::pmr::polymorphic_allocator<double>> const arr({2, 3}, &mbr);
		BOOST_TEST( buffer[ 0] == '0' );  // buffer is intact when initializing without value
		BOOST_TEST( buffer[13] == '3' );

		BOOST_TEST( arr.num_elements() == 2*3 );
		//  BOOST_TEST( arr[0][0] != 0.0 );
		//  BOOST_TEST( arr[1][2] != 0.0 );
	}
	{
		char buffer[] = "0123456789012345678901234567890123456789012345678901234567890123456789";  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) use raw memory

		std::pmr::monotonic_buffer_resource mbr(std::data(buffer), std::size(buffer));
		static_assert( std::size(buffer) > 6*sizeof(double) );

		multi::array<double, 2, std::pmr::polymorphic_allocator<double>> A({2, 3}, 0.0, &mbr);  // NOLINT(readability-identifier-length)
		//  BOOST_TEST( buffer[ 0] != '0' );  // buffer not is intact when initializing with value
		//  BOOST_TEST( buffer[13] != '3' );

		BOOST_TEST( A[0][0] == 0.0 );
		BOOST_TEST( A[1][2] == 0.0 );
	}
	{
		char buffer[] = "0123456789012345678901234567890123456789012345678901234567890123456789";  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) use raw memory

		std::pmr::monotonic_buffer_resource mbr(std::data(buffer), std::size(buffer));
		static_assert( std::size(buffer) > 6*sizeof(double) );

		multi::array<double, 2, std::pmr::polymorphic_allocator<double>> arr({2, 3}, {}, &mbr);
		//  BOOST_TEST( buffer[ 0] != '0' );  // buffer not is intact when initializing with value
		//  BOOST_TEST( buffer[13] != '3' );

		BOOST_TEST( arr[0][0] == double{} );
		BOOST_TEST( arr[1][2] == double{} );
	}
	{
		char buffer[] = "0123456789012345678901234567890123456789012345678901234567890123456789";  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) use raw memory

		std::pmr::monotonic_buffer_resource mbr{std::data(buffer), std::size(buffer)};
		static_assert( std::size(buffer) > 6*sizeof(double) );

		multi::array<double, 2, std::pmr::polymorphic_allocator<double>> arr({2, 3}, 666.0, &mbr);
		//  BOOST_TEST( buffer[ 0] != '0' );  // buffer not is intact when initializing with value
		//  BOOST_TEST( buffer[13] != '3' );

		BOOST_TEST( arr[0][0] == 666.0 );
		BOOST_TEST( arr[1][2] == 666.0 );
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
				multi::extensions_t<2>{1000 - idx%10, 1000 + idx%10},  // MSVC needs multi::extensions_t<2>
				resp
			);
			std::fill_n(arr.data_elements(), arr.num_elements(), 1);
			auto* be = arr.data_elements();
			decltype(be) en = arr.data_elements() + arr.num_elements();
			return std::accumulate(be, en, int64_t{}, std::plus<int64_t>{});
		}
	);

	auto time = std::chrono::high_resolution_clock::now() - start_time;
	std::cout<< time.count() / count <<"          "<< acc << '\n';
}
#endif
#endif
