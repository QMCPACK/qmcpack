// Copyright 2023-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// #include <boost/test/unit_test.hpp>

#include <boost/mpl/list.hpp>

#include "../../complex.hpp"
#include <boost/multi/array.hpp>

#include <type_traits>

namespace multi = boost::multi;

using float_types = boost::mpl::list<float, double>;

BOOST_AUTO_TEST_CASE_TEMPLATE(complex_ctors, T, float_types) {
	{
		multi::complex<T> const zeta = T{1.0} + multi::imaginary<T>{T{2.0}};
		BOOST_REQUIRE( zeta.real() == T{1.0});
		BOOST_REQUIRE( zeta.imag() == T{2.0});
	}
	// {
	//  multi::complex<T> zeta = T{1.0} + T{2.0} * multi::imaginary<T>::i;
	//  BOOST_REQUIRE( zeta.real() == T{1.0});
	//  BOOST_REQUIRE( zeta.imag() == T{2.0});
	// }
	//  {
	//      multi::complex<T> zeta = T{1.0} + multi::imaginary{T{2.0}};
	//      BOOST_REQUIRE( zeta.real() == T{1.0});
	//      BOOST_REQUIRE( zeta.imag() == T{2.0});
	//  }
}

BOOST_AUTO_TEST_CASE(double_complex_literals) {
	using multi::literals::operator""_I;
	multi::complex<double> const zeta = 1.0 + 2.0_I;
	//  multi::complex<double> zeta = 1.0 + 2.0i;  // literal i is not standard

	BOOST_REQUIRE( zeta.real() == 1.0 );
	BOOST_REQUIRE( zeta.imag() == 2.0 );
}

BOOST_AUTO_TEST_CASE(imaginary_equal) {
	using multi::literals::operator""_I;
	multi::imaginary<double> const zeta = 2.0_I;

	BOOST_REQUIRE( zeta == multi::imaginary<double>{2.0} );
}

BOOST_AUTO_TEST_CASE(imaginary_assign) {
	using multi::literals::operator""_I;
	multi::imaginary<double> zeta;  // NOLINT(cppcoreguidelines-pro-type-member-init,hicpp-member-init)
	zeta = 2.0_I;

	BOOST_REQUIRE( zeta == multi::imaginary<double>{2.0} );
}

BOOST_AUTO_TEST_CASE(float_complex_literals) {
	using multi::literals::operator""_IF;
	//  multi::complex<float> const zeta = 1.0f + 2.0  _i;  // may induced an undesired or forbidden conversion
	//  multi::complex<float> const zeta = 1.0f + 2.0 f_i;  // literal f_i is not standard
	//  multi::complex<float> const zeta = 1.0f + 2.0_f_i;
	multi::complex<float> const zeta = 1.0F + 2.0_IF;

	BOOST_REQUIRE( zeta.real() == 1.0F );
	BOOST_REQUIRE( zeta.imag() == 2.0F );
}

BOOST_AUTO_TEST_CASE(float_complex_assignment) {
	using multi::literals::operator""_IF;
	multi::complex<float> zeta;  // NOLINT(cppcoreguidelines-pro-type-member-init,hicpp-member-init)

	zeta = 1.0F + 2.0_IF;
	BOOST_REQUIRE( zeta.real() == 1.0F );
	BOOST_REQUIRE( zeta.imag() == 2.0F );

	zeta = 1.0F;
	BOOST_REQUIRE( zeta.real() == 1.0F );
	BOOST_REQUIRE( zeta.imag() == 0.0F );
}

BOOST_AUTO_TEST_CASE(float_complex_aggregate) {
	static_assert(std::is_aggregate_v<multi::complex<float>>);

	// auto const c = multi::complex<float>{._real = 1.0, ._imag = 2.0};

	// BOOST_REQUIRE( real(zeta) == 1.0F );
	// BOOST_REQUIRE( imag(zeta) == 2.0F );
}

BOOST_AUTO_TEST_CASE(double_complex_abs) {
	using multi::literals::operator""_I;
	multi::complex<double> const zeta = 1.0 + 2.0_I;

	BOOST_REQUIRE( abs(zeta) <= std::max(zeta.real(), zeta.imag()) );
}

BOOST_AUTO_TEST_CASE(double_complex_plus_eq) {
	using multi::literals::operator""_I;
	multi::complex<double>       zeta = 1.0 + 2.0_I;
	multi::complex<double> const yeta = 1.0 + 2.0_I;

	zeta += yeta;

	BOOST_REQUIRE( zeta == 2.0 * yeta );
	BOOST_REQUIRE( zeta == yeta / 0.5 );
}

// BOOST_AUTO_TEST_CASE(complex_member_cast) {
//  multi::array<multi::complex<double>, 2> A = {
//      {  {1., 2.}, {3., 4.}},
//      {{22., 33.}, {5., 9.}},
//  };

//  {
//      auto&& Areal = A.member_cast<double>(&multi::complex<double>::re);
//      auto&& Aimag = A.member_cast<double>(&multi::complex<double>::im);

//      BOOST_REQUIRE(Areal[1][0] == 22.);
//      BOOST_REQUIRE(Aimag[1][0] == 33.);
//  }
//  {
//      auto&& Areal = A.member_cast<double>(&multi::complex<double>::re);
//      auto&& Aimag = A.member_cast<double>(&multi::complex<double>::im);

//      BOOST_REQUIRE(Areal[1][0] == 22.);
//      BOOST_REQUIRE(Aimag[1][0] == 33.1);
//  }
// }
