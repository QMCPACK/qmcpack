// Copyright 2019-2025 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <boost/core/lightweight_test.hpp>

#include <array>    // for array
#include <cmath>    // for abs  // IWYU pragma: keep
#include <cstddef>  // for _LIBCPP_VERSION  // IWYU pragma: keep  // NOLINT(misc-include-cleaner)
// IWYU pragma: no_include <cstdlib>                          // for abs
#include <complex>   // for complex, operator==
#include <iterator>  // for data

#if __has_include(<version>)
#include <version>  // for __GLIBCXX__  // NOLINT(misc-include-cleaner)
#endif

#ifdef BOOST_MULTI_HAS_MEMORY_RESOURCE
#include <memory_resource>  // for monotonic_buffer_resource
#endif

#include <type_traits>  // for is_trivially_default_cons...  // IWYU pragma: keep

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

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	// BOOST_AUTO_TEST_CASE(pmr_double)
	{
		multi::array<std::complex<double>, 2> Aarr({2, 2}, std::complex<double>(4.0, 5.0));
		BOOST_TEST(Aarr[0][0] == std::complex<double>(4.0, 5.0) );
	}

#ifdef BOOST_MULTI_HAS_MEMORY_RESOURCE
	// BOOST_AUTO_TEST_CASE(pmr_double_uninitialized)
	{
		std::array<double, 12> buffer = {
			{4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.00, 11.0, 996.0, 997.0, 998.0, 999.0}
		};
		std::pmr::monotonic_buffer_resource pool(static_cast<void*>(std::data(buffer)), 12 * sizeof(double));

		multi::pmr::array<double, 2> Aarr({2, 2}, &pool);

		BOOST_TEST( std::abs( buffer[0] - 4.0 ) < 1E-6 );
		BOOST_TEST( std::abs( buffer[1] - 5.0 ) < 1E-6 );

#ifdef __GLIBCXX__
		BOOST_TEST         ( &Aarr[0][0] == buffer.data() );
		BOOST_TEST( std::abs( Aarr[0][0] - 4.0 ) < 1E-6);
#elif defined(_LIBCPP_VERSION)
		BOOST_TEST         ( &Aarr[0][0] == &buffer[buffer.size() - 4] );  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
		BOOST_TEST( std::abs( Aarr[0][0] - 996.0 ) < 1E-6 );
#endif
	}
	{
		std::array<double, 12> buffer = {
			{4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.00, 11.0, 996.0, 997.0, 998.0, 999.0}
		};
		std::pmr::monotonic_buffer_resource pool(static_cast<void*>(std::data(buffer)), 12 * sizeof(double));

		multi::pmr::array<double, 2> Aarr({2, 2}, double{}, &pool);

#ifdef __GLIBCXX__
		BOOST_TEST( std::abs( buffer[0] - 0.0 ) < 1E-6 );
		BOOST_TEST( std::abs( buffer[1] - 0.0 ) < 1E-6 );

		BOOST_TEST( &Aarr[0][0] == buffer.data() );
#elif defined(_LIBCPP_VERSION)
		BOOST_TEST( std::abs( buffer[0] - 4.0 ) < 1E-6);
		BOOST_TEST( std::abs( buffer[1] - 5.0 ) < 1E-6 );
		BOOST_TEST( std::abs( buffer[buffer.size()-4] - 0.0 ) < 1E-6 );
		BOOST_TEST( std::abs( buffer[buffer.size()-3] - 0.0 ) < 1E-6 );
		BOOST_TEST( std::abs( buffer[buffer.size()-5] - 11.0 ) < 1E-6 );

		BOOST_TEST( &Aarr[0][0] == &buffer[buffer.size() - 4] );  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
#endif

		BOOST_TEST( std::abs( Aarr[0][0] - 0.0 ) < 1E-6);
	}

	// BOOST_AUTO_TEST_CASE(pmr_complex_initialized_2)
	{
		std::array<double, 12> buffer = {
			{4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.00, 11.0, 996.0, 997.0, 998.0, 999.0}
		};
		std::pmr::monotonic_buffer_resource pool{static_cast<void*>(std::data(buffer)), 12 * sizeof(double)};

		multi::pmr::array<std::complex<double>, 2> Aarr({2, 2}, &pool);

#ifdef __GLIBCXX__
		BOOST_TEST( std::abs( buffer[0] - 4.0 ) < 1E-6 );
		BOOST_TEST( std::abs( buffer[1] - 5.0 ) < 1E-6 );
		BOOST_TEST( Aarr[0][0] == std::complex<double>(4.0, 5.0) );
#elif defined(_LIBCPP_VERSION)
		BOOST_TEST( std::abs( buffer[buffer.size() - 4] - 996.0 ) < 1E-6 );
		BOOST_TEST( std::abs( buffer[buffer.size() - 3] - 997.0 ) < 1E-6 );
		BOOST_TEST( std::abs( Aarr[0][0].real() - 8.0 ) < 1E-6 );
		BOOST_TEST( std::abs( Aarr[0][0].imag() - 9.0 ) < 1E-6 );
#endif
		Aarr[0][0] = std::complex<double>{40.0, 50.0};

#ifdef __GLIBCXX__
		BOOST_TEST( std::abs( buffer[0] - 40.0 ) < 1E-6 );
		BOOST_TEST( std::abs( buffer[1] - 50.0 ) < 1E-6 );
#elif defined(_LIBCPP_VERSION)
		BOOST_TEST( std::abs( buffer[buffer.size() - 4] - 996.0 ) < 1E-6 );
		BOOST_TEST( std::abs( buffer[buffer.size() - 3] - 997.0 ) < 1E-6 );
#endif
	}

	// BOOST_AUTO_TEST_CASE(pmr_complex_initialized_4)
	{
		std::array<double, 12> buffer = {
			{4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.00, 11.0, 999.9, 999.9, 999.9, 999.9}
		};
		std::pmr::monotonic_buffer_resource pool{static_cast<void*>(std::data(buffer)), 12 * sizeof(double)};

		multi::pmr::array<std::complex<double>, 2> Aarr({2, 2}, &pool);

#ifdef __GLIBCXX__
		BOOST_TEST( std::abs( Aarr[0][0].real() - 4.0 ) < 1E-6 );
		BOOST_TEST( std::abs( Aarr[0][0].imag() - 5.0 ) < 1E-6 );
#elif defined(_LIBCPP_VERSION)
		BOOST_TEST( std::abs( Aarr[0][0].real() - 8.0 ) < 1E-6 );
		BOOST_TEST( std::abs( Aarr[0][0].imag() - 9.0 ) < 1E-6 );
#endif

		BOOST_TEST( std::abs( buffer[0] - 4.0 ) < 1E-6 );
		BOOST_TEST( std::abs( buffer[1] - 5.0 ) < 1E-6 );

#ifdef __GLIBCXX__
		BOOST_TEST( static_cast<void*>(buffer.data()) == static_cast<void*>(&Aarr[0][0]) );
#elif defined(_LIBCPP_VERSION)
		BOOST_TEST( static_cast<void*>(&buffer[4]) == static_cast<void*>(&Aarr[0][0]) );  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
#endif
	}

	// BOOST_AUTO_TEST_CASE(pmr_complex_initialized_3)
	{
		std::array<double, 12> buffer = {
			{4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.00, 11.0, 996.0, 997.0, 998.0, 999.0}
		};
		std::pmr::monotonic_buffer_resource pool{static_cast<void*>(std::data(buffer)), 12 * sizeof(double)};

		multi::pmr::array<std::complex<double>, 2> const Aarr({2, 2}, std::complex<double>{40.0, 50.0}, &pool);

		BOOST_TEST( std::abs( Aarr[0][0].real() - 40.0 ) < 1E-6 );
		BOOST_TEST( std::abs( Aarr[0][0].imag() - 50.0 ) < 1E-6 );

#ifdef __GLIBCXX__
		BOOST_TEST( std::abs( buffer[0] - 40.0 ) < 1E-6 );
		BOOST_TEST( std::abs( buffer[1] - 50.0 ) < 1E-6 );
#elif defined(_LIBCPP_VERSION)
		BOOST_TEST( std::abs( buffer[buffer.size() - 4] - 40.0 ) < 1E-6 );
		BOOST_TEST( std::abs( buffer[buffer.size() - 3] - 50.0 ) < 1E-6 );
#endif
	}

	// BOOST_AUTO_TEST_CASE(pmr_complex_initialized)
	{
		std::array<double, 12> buffer = {
			{4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.00, 11.0, 996.0, 997.0, 998.0, 999.0}
		};
		std::pmr::monotonic_buffer_resource pool{static_cast<void*>(std::data(buffer)), 12 * sizeof(double)};

		multi::pmr::array<std::complex<double>, 2> Aarr({2, 2}, &pool);

		if constexpr(multi::force_element_trivial_default_construction<std::complex<double>>) {
			BOOST_TEST( std::abs( buffer[0] - 4.0 ) < 1E-6 );
			BOOST_TEST( std::abs( buffer[1] - 5.0 ) < 1E-6 );

#ifdef __GLIBCXX__
			BOOST_TEST(Aarr[0][0] == std::complex<double>(4.0, 5.0) );
#elif defined(_LIBCPP_VERSION)
			BOOST_TEST( std::abs( Aarr[0][0].real() - 8.0 ) < 1E-6 );
			BOOST_TEST( std::abs( Aarr[0][0].imag() - 9.0 ) < 1E-6 );
#endif
		} else {
			BOOST_TEST( std::abs( buffer[0] - 0.0 ) < 1E-6 );
			BOOST_TEST( std::abs( buffer[1] - 0.0 ) < 1E-6);

			BOOST_TEST( std::abs( Aarr[0][0] - 0.0 ) < 1E-6);
		}
	}
#endif
	return boost::report_errors();
}
