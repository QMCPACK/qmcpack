// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2023 Alfredo A. Correa

// #define BOOST_TEST_MODULE "C++ Unit Tests for Multi complex"  // NOLINT(cppcoreguidelines-macro-usage) title
#include <boost/test/unit_test.hpp>

#include <multi/array.hpp>

#include <complex>

#if(MULTI_PROVIDES_PMR_ARRAY)
#include <memory_resource>  // for polymorphic memory resource, monotonic buffer
#endif

#include <vector>

namespace multi = boost::multi;

#ifdef __NVCC__
template<>
inline constexpr bool multi::force_element_trivial_default_construction<std::complex<double>> = true;
template<>
inline constexpr bool multi::force_element_trivial_default_construction<std::complex<float>> = true;
#else
// vvv nvcc (12.1?) doesn't tolerate this kind of customization: "error: expected initializer before ‘<’"
template<class T>
inline constexpr bool multi::force_element_trivial_default_construction<std::complex<T>> = std::is_trivially_default_constructible<T>::value;
#endif

BOOST_AUTO_TEST_CASE(pmr_double) {
	multi::array<std::complex<double>, 2> Aarr({2, 2}, std::complex<double>(4.0, 5.0));
	BOOST_REQUIRE(Aarr[0][0] == std::complex<double>(4.0, 5.0) );
}

#if(MULTI_PROVIDES_PMR_ARRAY)
BOOST_AUTO_TEST_CASE(pmr_double_uninitialized) {
	std::array<double, 12> buffer = {{4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.00, 11.0,  999.9, 999.9, 999.9, 999.9}};
	std::pmr::monotonic_buffer_resource pool{static_cast<void*>(std::data(buffer)), 12*sizeof(double)};

	multi::pmr::array<double, 2> Aarr({2, 2}, &pool);

	BOOST_TEST( buffer[0] == 4.0 );
	BOOST_TEST( buffer[1] == 5.0 );

	BOOST_REQUIRE(Aarr[0][0] == 4.0);
}

BOOST_AUTO_TEST_CASE(pmr_complex_initialized_2) {
	std::array<double, 12> buffer = {{4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.00, 11.0, 999.9, 999.9, 999.9, 999.9}};
	std::pmr::monotonic_buffer_resource pool{static_cast<void*>(std::data(buffer)), 12*sizeof(double)};

	multi::pmr::array<std::complex<double>, 2> Aarr({2, 2}, &pool);

	// BOOST_TEST( buffer[0] == 4.0 );
	// BOOST_TEST( buffer[1] == 5.0 );

	// BOOST_REQUIRE(Aarr[0][0] == std::complex<double>(4.0, 5.0) );

	Aarr[0][0] = std::complex<double>{40.0, 50.0};
	BOOST_TEST( buffer[0] == 40.0 );
	BOOST_TEST( buffer[1] == 50.0 );
}

BOOST_AUTO_TEST_CASE(pmr_complex_initialized_4) {
	std::array<double, 12> buffer = {{4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.00, 11.0, 999.9, 999.9, 999.9, 999.9}};
	std::pmr::monotonic_buffer_resource pool{static_cast<void*>(std::data(buffer)), 12*sizeof(double)};

	multi::pmr::array<std::complex<double>, 2> Aarr({2, 2}, &pool);

	BOOST_REQUIRE(Aarr[0][0] == std::complex<double>(4.0, 5.0) );

	BOOST_TEST( buffer[0] == 4.0 );
	BOOST_TEST( buffer[1] == 5.0 );

	BOOST_TEST( static_cast<void*>(buffer.data()) == static_cast<void*>(&Aarr[0][0]) );
}

BOOST_AUTO_TEST_CASE(pmr_complex_initialized_3) {
	std::array<double, 12> buffer = {{4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.00, 11.0, 999.9, 999.9, 999.9, 999.9}};
	std::pmr::monotonic_buffer_resource pool{static_cast<void*>(std::data(buffer)), 12*sizeof(double)};

	multi::pmr::array<std::complex<double>, 2> const Aarr({2, 2}, std::complex<double>{40.0, 50.0}, &pool);

	BOOST_TEST( Aarr[0][0] == (std::complex<double>{40.0, 50.0}) );

	BOOST_TEST( buffer[0] == 40.0 );
	BOOST_TEST( buffer[1] == 50.0 );
}

BOOST_AUTO_TEST_CASE(pmr_complex_initialized) {
	std::array<double, 12> buffer = {{4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.00, 11.0, 999.9, 999.9, 999.9, 999.9}};
	std::pmr::monotonic_buffer_resource pool{static_cast<void*>(std::data(buffer)), 12*sizeof(double)};

	multi::pmr::array<std::complex<double>, 2> Aarr({2, 2}, &pool);

	if constexpr(multi::force_element_trivial_default_construction<std::complex<double>>) {
		BOOST_TEST( buffer[0] == 4.0 );
		BOOST_TEST( buffer[1] == 5.0 );

		BOOST_REQUIRE(Aarr[0][0] == std::complex<double>(4.0, 5.0) );
	} else {
		BOOST_TEST( buffer[0] == 0.0 );
		BOOST_TEST( buffer[1] == 0.0 );

		BOOST_REQUIRE(Aarr[0][0] == 0.0);
	}
}
#endif
