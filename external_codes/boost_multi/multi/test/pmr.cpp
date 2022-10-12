// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2022 Alfredo A. Correa

#define BOOST_TEST_MODULE "C++ Unit Tests for Multi pmr allocators"
#include<boost/test/unit_test.hpp>

#include "multi/array.hpp"

#include <memory_resource>  // for polymorphic memory resource, monotonic buffer

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(pmr_partially_formed) {

	{
		char buffer[] = "0123456789012345678901234567890123456789012345678901234567890123456789";  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) use raw memory

		std::pmr::monotonic_buffer_resource mbr{std::data(buffer), std::size(buffer)};
		static_assert( std::size(buffer) > 6*sizeof(double) );

		multi::array<double, 2, std::pmr::polymorphic_allocator<double>> A({2, 3}, &mbr);  // NOLINT(readability-identifier-length)
		BOOST_TEST( buffer[ 0] == '0' );  // buffer is intact when initializing without value
		BOOST_TEST( buffer[13] == '3' );

		BOOST_TEST( A.num_elements() == 2*3 );
	//  BOOST_TEST( A[0][0] != 0. );
	//  BOOST_TEST( A[1][2] != 0. );
	}
	{
		char buffer[] = "0123456789012345678901234567890123456789012345678901234567890123456789";  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) use raw memory

		std::pmr::monotonic_buffer_resource mbr{std::data(buffer), std::size(buffer)};
		static_assert( std::size(buffer) > 6*sizeof(double) );

		multi::array<double, 2, std::pmr::polymorphic_allocator<double>> A({2, 3}, 0., &mbr);  // NOLINT(readability-identifier-length)
	//  BOOST_TEST( buffer[ 0] != '0' );  // buffer not is intact when initializing with value
	//  BOOST_TEST( buffer[13] != '3' );

		BOOST_TEST( A[0][0] == 0. );
		BOOST_TEST( A[1][2] == 0. );
	}
	{
		char buffer[] = "0123456789012345678901234567890123456789012345678901234567890123456789";  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) use raw memory

		std::pmr::monotonic_buffer_resource mbr{std::data(buffer), std::size(buffer)};
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
}
