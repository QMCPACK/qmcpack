// Copyright 2019-2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#if defined(_MSC_VER)
	#pragma warning(push)
	#pragma warning(disable : 4244)  // allow conversion from double to float in uninitialized_construct algorithms
#endif

#include <boost/multi/array.hpp>

#include <complex>
// #include <cmath>  // for abs  // IWYU pragma: keep
#include <cstdlib>                          // for abs

namespace multi = boost::multi;

void fun(multi::array<std::complex<float>, 2> arr);
void fun(multi::array<std::complex<float>, 2> arr) { arr.clear(); }

void gun(multi::array<std::complex<float>, 2> const& /*unused*/);
void gun(multi::array<std::complex<float>, 2> const& /*unused*/) {
	/* no-op */
}

#include <boost/core/lightweight_test.hpp>
#define BOOST_AUTO_TEST_CASE(CasenamE) /**/

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
// NOLINTBEGIN(fuchsia-default-arguments-calls)  // std::complex has a constructor with a default argument, not in the library
BOOST_AUTO_TEST_CASE(complex_conversion_float_to_double) {
	std::complex<float> const cee{1.0, 2.0};

	std::complex<double> const zee = cee;

	static_assert(multi::detail::is_explicitly_convertible_v<std::complex<float>, std::complex<double>>);
	static_assert(multi::detail::is_implicitly_convertible_v<std::complex<float>, std::complex<double>>);

	BOOST_TEST( std::abs( cee.real() - static_cast<float>(zee.real())) < 1E-6F );

	multi::static_array<std::complex<float>, 1> const  CEE1(10, std::complex<float>{});  // NOLINT(fuchsia-default-arguments-calls)
	multi::static_array<std::complex<double>, 1> const ZEE1 = CEE1;
}

BOOST_AUTO_TEST_CASE(complex_conversion_double_to_float) {
	std::complex<double> const zee{1.0, 2.0};

	static_assert(multi::detail::is_explicitly_convertible_v<std::complex<double>, std::complex<float>>);
	static_assert(!multi::detail::is_implicitly_convertible_v<std::complex<double>, std::complex<float>>);

	std::complex<float> const cee{zee};

	BOOST_TEST( std::abs( cee.real() - static_cast<float>(zee.real()) ) < 1E-6F );

	multi::static_array<std::complex<double>, 1> const ZEE1(10, std::complex<float>{});
	multi::static_array<std::complex<float>, 1> const  CEE1{ZEE1};
}

BOOST_AUTO_TEST_CASE(double_to_complex_conversion_documentation) {
	// conversions from real to complex is implicit ...
	double const               dee = 5.0;
	std::complex<double> const zee = dee;

	BOOST_TEST( std::abs( zee.real() - 5.0 ) < 1E-6 );
	BOOST_TEST( std::abs( zee.imag() - 0.0 ) < 1E-6 );

	// ... therefore from array of reals to arrays of complex is also
	multi::array<double, 2>               DEE({10, 10}, dee);
	multi::array<std::complex<double>, 2> ZEE = DEE;

	BOOST_TEST( std::abs( ZEE[3][4].real() - 5.0 ) < 1E-6);
	BOOST_TEST( std::abs( ZEE[3][4].imag() - 0.0 ) < 1E-6);

	multi::array<std::complex<double>, 2> ZEE2{DEE};

	BOOST_TEST( std::abs( ZEE2[3][4].real() - 5.0 ) < 1E-6);
	BOOST_TEST( std::abs( ZEE2[3][4].imag() - 0.0 ) < 1E-6);

	// multi::array<double, 2> DEE2{ZEE};  // compilation error, good
}

BOOST_AUTO_TEST_CASE(conversion_in_function_call) {
	multi::array<std::complex<double>, 2> ZEE({10, 10});
	fun(multi::array<std::complex<float>, 2>{ZEE});
	gun(multi::array<std::complex<float>, 2>{ZEE});
}

BOOST_AUTO_TEST_CASE(float_to_double) {
	float const dee = 5.0F;
	// float const eff{dee};  // -Wc++11-narrowing  // NOLINT(bugprone-narrowing-conversions)
	// float const eff = dee;  // NOLINT(bugprone-narrowing-conversions,cppcoreguidelines-narrowing-conversions)
	// float const eff(dee);  // NOLINT(bugprone-narrowing-conversions,cppcoreguidelines-narrowing-conversions)
	auto const eff = static_cast<double>(dee);

	BOOST_TEST( std::abs( eff - 5.0) < 1E-6 );

	multi::array<float, 2> const DEE({10, 10}, dee);
	// multi::array<float, 2> const EFF(DEE);
	auto const EFF = static_cast<multi::array<double, 2>>(DEE);  // TODO(correaa) investigate producing intermediate types accessible through interminediate types

	BOOST_TEST( std::abs( EFF[3][4] - 5.0 ) < 1E-6 );

	// multi::array<float, 2> const EFF = DEE;
}

BOOST_AUTO_TEST_CASE(double_to_float) {
	double const dee = 5.0;
	// float const eff{dee};  // -Wc++11-narrowing  // NOLINT(bugprone-narrowing-conversions)
	// float const eff = dee;  // NOLINT(bugprone-narrowing-conversions,cppcoreguidelines-narrowing-conversions)
	// float const eff(dee);  // NOLINT(bugprone-narrowing-conversions,cppcoreguidelines-narrowing-conversions)
	auto const eff = static_cast<float>(dee);

	BOOST_TEST( std::abs( eff - 5.0F) < 1E-6F );

	multi::array<double, 2> const DEE({10, 10}, dee);

	// multi::array<float, 2> const EFF(DEE);
	auto const EFF = static_cast<multi::array<float, 2>>(DEE);  // TODO(correaa) investigate producing intermediate types accessible through interminediate types

	BOOST_TEST( std::abs( EFF[3][4] - 5.0F ) < 1E-6F );

	// multi::array<float, 2> const EFF = DEE;
}

BOOST_AUTO_TEST_CASE(complex_to_complex_conversion) {
	std::complex<float> const  cee{1.0, 2.0};
	std::complex<double> const zee = cee;

	BOOST_TEST( std::abs( zee.real() - 1.0) < 1E-6);
	BOOST_TEST( std::abs( zee.imag() - 2.0) < 1E-6);

	// std::complex<float> cee2 = zee;  // implicit conversion, compilation error
	std::complex<float> const cee2{zee};

	BOOST_TEST( std::abs( cee2.real() - 1.0F ) < 1E-6F );
	BOOST_TEST( std::abs( cee2.imag() - 2.0F ) < 1E-6F );

	multi::array<std::complex<float>, 2> const  CEE({10, 10}, cee);
	multi::array<std::complex<double>, 2> const ZEE = CEE;

	BOOST_TEST( std::abs( ZEE[3][4].real() - 1.0 ) < 1E-6 );
	BOOST_TEST( std::abs( ZEE[3][4].imag() - 2.0 ) < 1E-6 );

	// multi::array<std::complex<float>, 2> const CEE2 = ZEE;  // implicit conversion, compilation error
	multi::array<std::complex<float>, 2> const CEE2{ZEE};

	BOOST_TEST( std::abs( CEE2[3][4].real() - 1.0F ) < 1E-6F );
	BOOST_TEST( std::abs( CEE2[3][4].imag() - 2.0F ) < 1E-6F );
}
// NOLINTEND(fuchsia-default-arguments-calls)
return boost::report_errors();}

#if defined(_MSC_VER)
	#pragma warning(pop)
#endif
