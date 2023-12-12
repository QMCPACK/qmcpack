// Copyright 2019-2023 Alfredo A. Correa

#include<boost/test/unit_test.hpp>

#include <multi/array.hpp>

#include <numeric>

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(pmr_partially_formed) {
#if defined(__cpp_lib_memory_resource) and (__cpp_lib_memory_resource >= 201603)
	{
		char buffer[] = "0123456789012345678901234567890123456789012345678901234567890123456789";  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) use raw memory

		std::pmr::monotonic_buffer_resource mbr{std::data(buffer), std::size(buffer)};
		static_assert( std::size(buffer) > 6*sizeof(double) );

		multi::array<double, 2, std::pmr::polymorphic_allocator<double>> const A({2, 3}, &mbr);  // NOLINT(readability-identifier-length)
		BOOST_TEST( buffer[ 0] == '0' );  // buffer is intact when initializing without value
		BOOST_TEST( buffer[13] == '3' );

		BOOST_TEST( A.num_elements() == 2*3 );
	//  BOOST_TEST( A[0][0] != 0. );
	//  BOOST_TEST( A[1][2] != 0. );
	}
	{
		char buffer[] = "0123456789012345678901234567890123456789012345678901234567890123456789";  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) use raw memory

		std::pmr::monotonic_buffer_resource mbr(std::data(buffer), std::size(buffer));
		static_assert( std::size(buffer) > 6*sizeof(double) );

		multi::array<double, 2, std::pmr::polymorphic_allocator<double>> A({2, 3}, 0.0, &mbr);  // NOLINT(readability-identifier-length)
	//  BOOST_TEST( buffer[ 0] != '0' );  // buffer not is intact when initializing with value
	//  BOOST_TEST( buffer[13] != '3' );

		BOOST_TEST( A[0][0] == 0. );
		BOOST_TEST( A[1][2] == 0. );
	}
	{
		char buffer[] = "0123456789012345678901234567890123456789012345678901234567890123456789";  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) use raw memory

		std::pmr::monotonic_buffer_resource mbr(std::data(buffer), std::size(buffer));
		static_assert( std::size(buffer) > 6*sizeof(double) );

		multi::array<double, 2, std::pmr::polymorphic_allocator<double>> A({2, 3}, {}, &mbr);  // NOLINT(readability-identifier-length)
	//  BOOST_TEST( buffer[ 0] != '0' );  // buffer not is intact when initializing with value
	//  BOOST_TEST( buffer[13] != '3' );

		BOOST_TEST( A[0][0] == double{} );
		BOOST_TEST( A[1][2] == double{} );
	}
	{
		char buffer[] = "0123456789012345678901234567890123456789012345678901234567890123456789";  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) use raw memory

		std::pmr::monotonic_buffer_resource mbr{std::data(buffer), std::size(buffer)};
		static_assert( std::size(buffer) > 6*sizeof(double) );

		multi::array<double, 2, std::pmr::polymorphic_allocator<double>> A({2, 3}, 666., &mbr);  // NOLINT(readability-identifier-length)
	//  BOOST_TEST( buffer[ 0] != '0' );  // buffer not is intact when initializing with value
	//  BOOST_TEST( buffer[13] != '3' );

		BOOST_TEST( A[0][0] == 666. );
		BOOST_TEST( A[1][2] == 666. );
	}
#endif
}

BOOST_AUTO_TEST_CASE(pmr_benchmark) {
#if defined(__cpp_lib_memory_resource) and (__cpp_lib_memory_resource >= 201603)
//  auto* resp = std::pmr::unsynchronized_pool_resource(std::pmr::get_default_resource());
	auto* resp = std::pmr::get_default_resource();

	auto count = 50;
	auto start_time = std::chrono::high_resolution_clock::now();

	multi::extension_t const exts{0, count};
	auto acc = std::transform_reduce(
		exts.begin(), exts.end(), int64_t{0},
		std::plus<>{},
		[&resp](auto idx) {
			multi::pmr::array<int64_t, 2> arr({1000 - idx%10, 1000 + idx%10}, resp);
			std::fill_n(arr.data_elements(), arr.num_elements(), 1);
			return std::accumulate(arr.data_elements(), arr.data_elements() + arr.num_elements(), 0L);
		}
	);

	auto time = std::chrono::high_resolution_clock::now() - start_time;
	std::cout<< time.count() / count <<"          "<< acc <<std::endl;
#endif
}
