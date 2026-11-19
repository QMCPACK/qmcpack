// Copyright 2022-2026 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wconversion"
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wshadow"
#pragma clang diagnostic ignored "-Wsign-conversion"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif

#ifdef __clang__
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#include <boost/multi/array.hpp>  // for operator!=, implicit...

#include <type_traits>  // for is_same_v, is_convertib...

namespace multi = boost::multi;

#include <boost/core/lightweight_test.hpp>

auto main() -> int {  // NOLINT(readability-function-cognitive-complexity,bugprone-exception-escape)
	// BOOST_AUTO_TEST_CASE(convertibles_1D)
	{
		using NDArray = multi::array<double, 1>;

		NDArray const nda;
		(void)nda;

		static_assert(std::is_same_v<NDArray::element, multi::array<double, 1>::value_type>);
		static_assert(std::is_same_v<NDArray::element_ref, multi::array<double, 1>::reference>);

		using NDRef = multi::array_ref<double, 1>;

		static_assert(std::is_convertible_v<NDRef, NDArray>);

#ifndef __NVCC__
		static_assert(std::is_convertible_v<NDRef::reference, NDRef::value_type>);
		static_assert(std::is_convertible_v<NDRef::const_reference, NDRef::value_type>);
#else
		static_assert(std::is_convertible<NDRef::reference, NDRef::value_type>::value);
		static_assert(std::is_convertible<NDRef::const_reference, NDRef::value_type>::value);
#endif

		static_assert(std::is_same_v<NDRef::element, multi::array<double, 1>::value_type>);
		static_assert(std::is_same_v<NDRef::element_ref, multi::array<double, 1>::reference>);
	}

	// BOOST_AUTO_TEST_CASE(convertibles_2D)
	{
		using NDArray = multi::array<double, 2>;

		NDArray const nda;
		(void)nda;

		static_assert(std::is_same_v<NDArray::element, multi::array<double, 1>::value_type>);
		static_assert(std::is_same_v<NDArray::element_ref, multi::array<double, 1>::reference>);

		using NDRef = multi::array_ref<double, 2>;

		static_assert(std::is_convertible_v<NDRef, NDArray>);

		static_assert(std::is_convertible_v<NDRef::reference, NDRef::value_type>);
		static_assert(std::is_convertible_v<NDRef::const_reference, NDRef::value_type>);

		static_assert(std::is_same_v<NDRef::element, multi::array<double, 1>::value_type>);
		static_assert(std::is_same_v<NDRef::element_ref, multi::array<double, 1>::reference>);
	}

	// BOOST_AUTO_TEST_CASE(convertibles_3D)
	{
		using NDArray = multi::array<double, 3>;

		NDArray const nda;
		(void)nda;

		static_assert(std::is_same_v<NDArray::element, multi::array<double, 1>::value_type>);
		static_assert(std::is_same_v<NDArray::element_ref, multi::array<double, 1>::reference>);

		using NDRef = multi::array_ref<double, 3>;

		static_assert(std::is_convertible_v<NDRef, NDArray>);

		static_assert(std::is_convertible_v<NDRef::reference, NDRef::value_type>);
		static_assert(std::is_convertible_v<NDRef::const_reference, NDRef::value_type>);

		static_assert(std::is_same_v<NDRef::element, multi::array<double, 1>::value_type>);
		static_assert(std::is_same_v<NDRef::element_ref, multi::array<double, 1>::reference>);
	}

	return boost::report_errors();
}
