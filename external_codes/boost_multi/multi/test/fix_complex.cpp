// Copyright 2019-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>
// #include <boost/multi/pmr.hpp>

#include <complex>
#include <vector>

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
#endif

#ifndef BOOST_TEST_MODULE
#  define BOOST_TEST_MAIN
#endif

#include <boost/test/unit_test.hpp>

namespace multi = boost::multi;

#ifdef __NVCC__
template<>
inline constexpr bool multi::force_element_trivial_default_construction<std::complex<double>> = true;
template<>
inline constexpr bool multi::force_element_trivial_default_construction<std::complex<float>> = true;
#else
// vvv nvcc (12.1?) doesn't tolerate this kind of customization: "error: expected initializer before ‘<’"
template<class T>
inline constexpr bool multi::force_element_trivial_default_construction<std::complex<T>> = std::is_trivially_default_constructible_v<T>;
#endif

BOOST_AUTO_TEST_CASE(pmr_double) {
	multi::array<std::complex<double>, 2> Aarr({2, 2}, std::complex<double>(4.0, 5.0));
	BOOST_REQUIRE(Aarr[0][0] == std::complex<double>(4.0, 5.0) );
}

#ifdef BOOST_MULTI_HAS_MEMORY_RESOURCE
BOOST_AUTO_TEST_CASE(pmr_double_uninitialized) {
	{
		std::array<double, 12> buffer = {{4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.00, 11.0,  996.0, 997.0, 998.0, 999.0}};
		std::pmr::monotonic_buffer_resource pool(static_cast<void*>(std::data(buffer)), 12*sizeof(double));

		multi::pmr::array<double, 2> Aarr({2, 2}, &pool);

		BOOST_TEST( buffer[0] == 4.0 );
		BOOST_TEST( buffer[1] == 5.0 );

	#if defined(__GLIBCXX__)
		BOOST_TEST( &Aarr[0][0] == buffer.data() );
		BOOST_TEST( Aarr[0][0] == 4.0);
	#elif defined(_LIBCPP_VERSION)
		BOOST_TEST( &Aarr[0][0] == buffer.data() + (buffer.size() - 4) );  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
		BOOST_TEST( Aarr[0][0] == 996.0);
	#endif
	}
	{
		std::array<double, 12> buffer = {{4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.00, 11.0, 996.0, 997.0, 998.0, 999.0}};
		std::pmr::monotonic_buffer_resource pool(static_cast<void*>(std::data(buffer)), 12*sizeof(double));

		multi::pmr::array<double, 2> Aarr({2, 2}, double{}, &pool);

	#if defined(__GLIBCXX__)
		BOOST_TEST( buffer[0] == 0.0 );
		BOOST_TEST( buffer[1] == 0.0 );
		BOOST_TEST( &Aarr[0][0] == buffer.data() );
	#elif defined(_LIBCPP_VERSION)
		BOOST_TEST( buffer[0] == 4.0 );
		BOOST_TEST( buffer[1] == 5.0 );
		BOOST_TEST( buffer[buffer.size()-4] ==  0.0 );
		BOOST_TEST( buffer[buffer.size()-3] ==  0.0 );
		BOOST_TEST( buffer[buffer.size()-5] == 11.0 );
		BOOST_TEST( &Aarr[0][0] == buffer.data() + (buffer.size() - 4) );  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
	#endif

		BOOST_TEST( Aarr[0][0] == 0.0);
	}
}

BOOST_AUTO_TEST_CASE(pmr_complex_initialized_2) {
	std::array<double, 12> buffer = {{4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.00, 11.0, 996.0, 997.0, 998.0, 999.0}};
	std::pmr::monotonic_buffer_resource pool{static_cast<void*>(std::data(buffer)), 12*sizeof(double)};

	multi::pmr::array<std::complex<double>, 2> Aarr({2, 2}, &pool);

#if defined(__GLIBCXX__)
	BOOST_TEST( buffer[0] == 4.0 );
	BOOST_TEST( buffer[1] == 5.0 );
	BOOST_REQUIRE(Aarr[0][0] == std::complex<double>(4.0, 5.0) );
#elif defined(_LIBCPP_VERSION)
	BOOST_TEST( buffer[buffer.size() - 4] == 996.0 );
	BOOST_TEST( buffer[buffer.size() - 3] == 997.0 );
	BOOST_TEST(Aarr[0][0].real() == 8.0 );
	BOOST_TEST(Aarr[0][0].imag() == 9.0 );
#endif
	Aarr[0][0] = std::complex<double>{40.0, 50.0};

#if defined(__GLIBCXX__)
	BOOST_TEST( buffer[0] == 40.0 );
	BOOST_TEST( buffer[1] == 50.0 );
#elif defined(_LIBCPP_VERSION)
	BOOST_TEST( buffer[buffer.size() - 4] == 996.0 );
	BOOST_TEST( buffer[buffer.size() - 3] == 997.0 );
#endif
}

BOOST_AUTO_TEST_CASE(pmr_complex_initialized_4) {
	std::array<double, 12> buffer = {{4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.00, 11.0, 999.9, 999.9, 999.9, 999.9}};
	std::pmr::monotonic_buffer_resource pool{static_cast<void*>(std::data(buffer)), 12*sizeof(double)};

	multi::pmr::array<std::complex<double>, 2> Aarr({2, 2}, &pool);

#if defined(__GLIBCXX__)
	BOOST_REQUIRE(Aarr[0][0] == std::complex<double>(4.0, 5.0) );
#elif defined(_LIBCPP_VERSION)
	BOOST_REQUIRE(Aarr[0][0] == std::complex<double>(8.0, 9.0) );
#endif

	BOOST_TEST( buffer[0] == 4.0 );
	BOOST_TEST( buffer[1] == 5.0 );

#if defined(__GLIBCXX__)
	BOOST_TEST( static_cast<void*>(buffer.data()) == static_cast<void*>(&Aarr[0][0]) );
#elif defined(_LIBCPP_VERSION)
	BOOST_TEST( static_cast<void*>(buffer.data() + 4) == static_cast<void*>(&Aarr[0][0]) );  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
#endif
}

BOOST_AUTO_TEST_CASE(pmr_complex_initialized_3) {
	std::array<double, 12> buffer = {{4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.00, 11.0, 996.0, 997.0, 998.0, 999.0}};
	std::pmr::monotonic_buffer_resource pool{static_cast<void*>(std::data(buffer)), 12*sizeof(double)};

	multi::pmr::array<std::complex<double>, 2> const Aarr({2, 2}, std::complex<double>{40.0, 50.0}, &pool);

	BOOST_TEST( Aarr[0][0] == (std::complex<double>{40.0, 50.0}) );

#if defined(__GLIBCXX__)
	BOOST_TEST( buffer[0] == 40.0 );
	BOOST_TEST( buffer[1] == 50.0 );
#elif defined(_LIBCPP_VERSION)
	BOOST_TEST( buffer[buffer.size() - 4] == 40.0 );
	BOOST_TEST( buffer[buffer.size() - 3] == 50.0 );
#endif
}

BOOST_AUTO_TEST_CASE(pmr_complex_initialized) {
	std::array<double, 12> buffer = {{4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.00, 11.0, 996.0, 997.0, 998.0, 999.0}};
	std::pmr::monotonic_buffer_resource pool{static_cast<void*>(std::data(buffer)), 12*sizeof(double)};

	multi::pmr::array<std::complex<double>, 2> Aarr({2, 2}, &pool);

	if constexpr(multi::force_element_trivial_default_construction<std::complex<double>>) {
		BOOST_TEST( buffer[0] == 4.0 );
		BOOST_TEST( buffer[1] == 5.0 );

	#if defined(__GLIBCXX__)
		BOOST_REQUIRE(Aarr[0][0] == std::complex<double>(4.0, 5.0) );
	#elif defined(_LIBCPP_VERSION)
		BOOST_TEST(Aarr[0][0].real() == 8.0 );
		BOOST_TEST(Aarr[0][0].imag() == 9.0 );
	#endif
	} else {
		BOOST_TEST( buffer[0] == 0.0 );
		BOOST_TEST( buffer[1] == 0.0 );

		BOOST_REQUIRE(Aarr[0][0] == 0.0);
	}
}
#endif
