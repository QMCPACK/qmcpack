// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2019-2023 Alfredo A. Correa

#include<boost/test/unit_test.hpp>

#include <multi/array.hpp>

#include<complex>

namespace multi = boost::multi;

// NOLINTBEGIN(fuchsia-default-arguments-calls)  // this is a defect in std::complex, not in the library
BOOST_AUTO_TEST_CASE(complex_conversion_float_to_double) {
	std::complex<float> const cee{1.0, 2.0};

	std::complex<double> const zee = cee;

	static_assert(    multi::detail::is_explicitly_convertible_v<std::complex<float>, std::complex<double>> );
	static_assert(    multi::detail::is_implicitly_convertible_v<std::complex<float>, std::complex<double>> );

	BOOST_TEST(cee.real() == zee.real());

	multi::static_array<std::complex<float>, 1> const CEE1(10, std::complex<float>{});  // NOLINT(fuchsia-default-arguments-calls)
	multi::static_array<std::complex<double>, 1> const ZEE1 = CEE1;
}

BOOST_AUTO_TEST_CASE(complex_conversion_double_to_float) {
	std::complex<double> const zee{1.0, 2.0};

	static_assert(    multi::detail::is_explicitly_convertible_v<std::complex<double>, std::complex<float>>);
	static_assert(not multi::detail::is_implicitly_convertible_v<std::complex<double>, std::complex<float>>);

	std::complex<float> const cee{zee};

	BOOST_TEST(cee.real() == zee.real());

	multi::static_array<std::complex<double>, 1> const ZEE1(10, std::complex<float>{});
	multi::static_array<std::complex<float>, 1> const CEE1{ZEE1};
}

BOOST_AUTO_TEST_CASE(double_to_complex_conversion_documentation) {
    // conversions from real to complex is implicit ...
    double const dee = 5.0;
    std::complex<double> const zee = dee;

    BOOST_REQUIRE( zee.real() == 5.0 );
    BOOST_REQUIRE( zee.imag() == 0.0 );

    // ... therefore from array of reals to arrays of complex is also
    multi::array<double, 2> DEE({10, 10}, dee);
    multi::array<std::complex<double>, 2> ZEE = DEE;

    BOOST_REQUIRE( ZEE[3][4].real() == 5.0 );
    BOOST_REQUIRE( ZEE[3][4].imag() == 0.0 );

    multi::array<std::complex<double>, 2> ZEE2{DEE};

    BOOST_REQUIRE( ZEE2[3][4].real() == 5.0 );
    BOOST_REQUIRE( ZEE2[3][4].imag() == 0.0 );

    // multi::array<double, 2> DEE2{ZEE};  // compilation error
}

void fun(multi::array<std::complex<float>, 2> /*unused*/);
void fun(multi::array<std::complex<float>, 2> /*unused*/) { }  // NOLINT(performance-unnecessary-value-param)

void gun(multi::array<std::complex<float>, 2> const& /*unused*/);
void gun(multi::array<std::complex<float>, 2> const& /*unused*/) { }

BOOST_AUTO_TEST_CASE(conversion_in_function_call) {
    multi::array<std::complex<double>, 2> ZEE({10, 10});
    fun( multi::array<std::complex<float>, 2>{ZEE} );
    gun( multi::array<std::complex<float>, 2>{ZEE} );
}

BOOST_AUTO_TEST_CASE(double_to_float) {
    double const dee = 5.0;
    // float const eff{dee};  // -Wc++11-narrowing  // NOLINT(bugprone-narrowing-conversions)
    // float const eff = dee;  // NOLINT(bugprone-narrowing-conversions,cppcoreguidelines-narrowing-conversions)
    // float const eff(dee);  // NOLINT(bugprone-narrowing-conversions,cppcoreguidelines-narrowing-conversions)
    auto const eff = static_cast<float>(dee);

    // BOOST_REQUIRE( eff == 5.0 );  // -Wdouble-promotion
    BOOST_REQUIRE( eff == 5.0F );

    multi::array<double, 2> const DEE({10, 10}, dee);
    // multi::array<float, 2> const EFF(DEE);
    auto const EFF = static_cast<multi::array<float, 2>>(DEE);  // TODO(correaa) investigate producing intermediate types accessible through interminediate types

    BOOST_REQUIRE( EFF[3][4] == 5.0F );

    // multi::array<float, 2> const EFF = DEE;
}

BOOST_AUTO_TEST_CASE(complex_to_complex_conversion) {
    std::complex<float> const cee{1.0, 2.0};
    std::complex<double> const zee = cee;

    BOOST_REQUIRE( zee.real() == 1.0 );
    BOOST_REQUIRE( zee.imag() == 2.0 );

    // std::complex<float> cee2 = zee;  // implicit conversion, compilation error
    std::complex<float> const cee2{zee};

    BOOST_REQUIRE( cee2.real() == 1.0F );
    BOOST_REQUIRE( cee2.imag() == 2.0F );

    multi::array<std::complex<float>, 2> const CEE({10, 10}, cee);
    multi::array<std::complex<double>, 2> const ZEE = CEE;

    BOOST_REQUIRE( ZEE[3][4].real() == 1.0 );
    BOOST_REQUIRE( ZEE[3][4].imag() == 2.0 );

    // multi::array<std::complex<float>, 2> const CEE2 = ZEE;  // implicit conversion, compilation error
    multi::array<std::complex<float>, 2> const CEE2{ZEE};

    BOOST_REQUIRE( CEE2[3][4].real() == 1.0F );
    BOOST_REQUIRE( CEE2[3][4].imag() == 2.0F );
}
// NOLINTEND(fuchsia-default-arguments-calls)
