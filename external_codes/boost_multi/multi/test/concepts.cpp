// Copyright 2022-2023 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#include <boost/multi/array.hpp>

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

// #include <boost/mp11.hpp>
#include <boost/mpl/list.hpp>

namespace multi = boost::multi;

// using NDArrays = boost::mp11::mp_list<  // fails with Boost.Test 1.67
using NDArrays = boost::mpl::list<
	multi::array<double, 1>,
	multi::array<double, 2>,
	multi::array<double, 3>
>;

BOOST_AUTO_TEST_CASE_TEMPLATE(convertibles, NDArray, NDArrays)
{
	static_assert( std::is_convertible_v<typename NDArray::      reference, typename NDArray::value_type> );
	static_assert( std::is_convertible_v<typename NDArray::const_reference, typename NDArray::value_type> );

	static_assert( std::is_same_v<typename NDArray::element_type, typename multi::array<double, 1>::value_type>);
	static_assert( std::is_same_v<typename NDArray::element_ref , typename multi::array<double, 1>::reference >);

	using NDRef = typename NDArray::ref;

	static_assert( std::is_convertible_v<NDRef, NDArray> );

	static_assert( std::is_convertible_v<typename NDRef::      reference, typename NDRef::value_type> );
	static_assert( std::is_convertible_v<typename NDRef::const_reference, typename NDRef::value_type> );

	static_assert( std::is_same_v<typename NDRef::element_type, typename multi::array<double, 1>::value_type> );
	static_assert( std::is_same_v<typename NDRef::element_ref , typename multi::array<double, 1>::reference > );
}
