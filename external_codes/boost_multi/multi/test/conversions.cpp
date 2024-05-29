// Copyright 2019-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

#include <complex>

// Suppress warnings from boost.test
#if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wold-style-cast"
#  pragma clang diagnostic ignored "-Wundef"
#  pragma clang diagnostic ignored "-Wconversion"
#  pragma clang diagnostic ignored "-Wsign-conversion"
// #  pragma clang diagnostic ignored "-Wfloat-equal"
#elif defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wold-style-cast"
#  pragma GCC diagnostic ignored "-Wundef"
#  pragma GCC diagnostic ignored "-Wconversion"
#  pragma GCC diagnostic ignored "-Wsign-conversion"
// #  pragma GCC diagnostic ignored "-Wfloat-equal"
#elif defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4244)
#endif

#ifndef BOOST_TEST_MODULE
#  define BOOST_TEST_MAIN
#endif

#include <boost/test/unit_test.hpp>

#if defined(__clang__)
#  pragma clang diagnostic pop
#elif defined(__GNUC__)
#  pragma GCC diagnostic pop
#elif defined(_MSC_VER)
#  pragma warning(pop)
#  pragma warning(disable : 4244)
#endif

namespace multi = boost::multi;

// NOLINTBEGIN(fuchsia-default-arguments-calls)  // this is a defect in std::complex, not in the library
BOOST_AUTO_TEST_CASE(complex_conversion_float_to_double) {
	std::complex<float> const cee{1.0, 2.0};

	std::complex<double> const zee = cee;

	static_assert(multi::detail::is_explicitly_convertible_v<std::complex<float>, std::complex<double>>);
	static_assert(multi::detail::is_implicitly_convertible_v<std::complex<float>, std::complex<double>>);

	BOOST_CHECK_CLOSE( cee.real(), static_cast<float>(zee.real()), 1E-6 );

	multi::static_array<std::complex<float>, 1> const  CEE1(10, std::complex<float>{});  // NOLINT(fuchsia-default-arguments-calls)
	multi::static_array<std::complex<double>, 1> const ZEE1 = CEE1;
}

BOOST_AUTO_TEST_CASE(complex_conversion_double_to_float) {
	std::complex<double> const zee{1.0, 2.0};

	static_assert( multi::detail::is_explicitly_convertible_v<std::complex<double>, std::complex<float>>);
	static_assert(!multi::detail::is_implicitly_convertible_v<std::complex<double>, std::complex<float>>);

	std::complex<float> const cee{zee};

	BOOST_CHECK_CLOSE( cee.real(), static_cast<float>(zee.real()) , 1E-6);

	multi::static_array<std::complex<double>, 1> const ZEE1(10, std::complex<float>{});
	multi::static_array<std::complex<float>, 1> const  CEE1{ZEE1};
}

BOOST_AUTO_TEST_CASE(double_to_complex_conversion_documentation) {
	// conversions from real to complex is implicit ...
	double const               dee = 5.0;
	std::complex<double> const zee = dee;

	BOOST_REQUIRE_CLOSE( zee.real(), 5.0, 1E-6 );
	BOOST_REQUIRE_CLOSE( zee.imag(), 0.0, 1E-6 );

	// ... therefore from array of reals to arrays of complex is also
	multi::array<double, 2>               DEE({10, 10}, dee);
	multi::array<std::complex<double>, 2> ZEE = DEE;

	BOOST_REQUIRE_CLOSE( ZEE[3][4].real(), 5.0, 1E-6 );
	BOOST_REQUIRE_CLOSE( ZEE[3][4].imag(), 0.0, 1E-6 );

	multi::array<std::complex<double>, 2> ZEE2{DEE};

	BOOST_REQUIRE_CLOSE( ZEE2[3][4].real(), 5.0, 1E-6);
	BOOST_REQUIRE_CLOSE( ZEE2[3][4].imag(), 0.0, 1E-6 );

	// multi::array<double, 2> DEE2{ZEE};  // compilation error, good
}

void fun(multi::array<std::complex<float>, 2> arr);
void fun(multi::array<std::complex<float>, 2> arr) { arr.clear(); }

void gun(multi::array<std::complex<float>, 2> const& /*unused*/);
void gun(multi::array<std::complex<float>, 2> const& /*unused*/) {
	/* no-op */
}

BOOST_AUTO_TEST_CASE(conversion_in_function_call) {
	multi::array<std::complex<double>, 2> ZEE({10, 10});
	fun(multi::array<std::complex<float>, 2>{ZEE});
	gun(multi::array<std::complex<float>, 2>{ZEE});
}

BOOST_AUTO_TEST_CASE(double_to_float) {
	double const dee = 5.0;
	// float const eff{dee};  // -Wc++11-narrowing  // NOLINT(bugprone-narrowing-conversions)
	// float const eff = dee;  // NOLINT(bugprone-narrowing-conversions,cppcoreguidelines-narrowing-conversions)
	// float const eff(dee);  // NOLINT(bugprone-narrowing-conversions,cppcoreguidelines-narrowing-conversions)
	auto const eff = static_cast<float>(dee);

	// BOOST_REQUIRE( eff == 5.0 );  // -Wdouble-promotion
	BOOST_REQUIRE_CLOSE( eff, 5.0F, 1E-6 );

	multi::array<double, 2> const DEE({10, 10}, dee);
	// multi::array<float, 2> const EFF(DEE);
	auto const EFF = static_cast<multi::array<float, 2>>(DEE);  // TODO(correaa) investigate producing intermediate types accessible through interminediate types

	BOOST_REQUIRE_CLOSE( EFF[3][4], 5.0F, 1E-6 );

	// multi::array<float, 2> const EFF = DEE;
}

BOOST_AUTO_TEST_CASE(complex_to_complex_conversion) {
	std::complex<float> const  cee{1.0, 2.0};
	std::complex<double> const zee = cee;

	BOOST_REQUIRE_CLOSE( zee.real(), 1.0, 1E-6 );
	BOOST_REQUIRE_CLOSE( zee.imag(), 2.0, 1E-6 );

	// std::complex<float> cee2 = zee;  // implicit conversion, compilation error
	std::complex<float> const cee2{zee};

	BOOST_REQUIRE_CLOSE( cee2.real(), 1.0F, 1E-6 );
	BOOST_REQUIRE_CLOSE( cee2.imag(), 2.0F, 1E-6 );

	multi::array<std::complex<float>, 2> const  CEE({10, 10}, cee);
	multi::array<std::complex<double>, 2> const ZEE = CEE;

	BOOST_REQUIRE_CLOSE( ZEE[3][4].real(), 1.0, 1E-6);
	BOOST_REQUIRE_CLOSE( ZEE[3][4].imag(), 2.0, 1E-6);

	// multi::array<std::complex<float>, 2> const CEE2 = ZEE;  // implicit conversion, compilation error
	multi::array<std::complex<float>, 2> const CEE2{ZEE};

	BOOST_REQUIRE_CLOSE( CEE2[3][4].real(), 1.0F, 1E-6 );
	BOOST_REQUIRE_CLOSE( CEE2[3][4].imag(), 2.0F, 1E-6 );
}
// NOLINTEND(fuchsia-default-arguments-calls)
