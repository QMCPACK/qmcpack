// Copyright 2024 Alfredo A. Correa
// Copyright 2024 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Test explicitly calls deprecated function
#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnositc ignored "-Wdeprecated-declarations"
#elif defined(__GNUC__)
	#pragma GCC diagnostic push
	#pragma GCC diagnositc ignored "-Wdeprecated-declarations"
#endif

#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wold-style-cast"
	#pragma clang diagnostic ignored "-Wundef"
	#pragma clang diagnostic ignored "-Wconversion"
	#pragma clang diagnostic ignored "-Wsign-conversion"
	#pragma clang diagnostic ignored "-Wunknown-pragmas"
	#pragma clang diagnostic ignored "-Wunused-variable"
#elif defined(__GNUC__)
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wold-style-cast"
	#pragma GCC diagnostic ignored "-Wundef"
	#pragma GCC diagnostic ignored "-Wconversion"
	#pragma GCC diagnostic ignored "-Wsign-conversion"
	#pragma GCC diagnostic ignored "-Wunknown-pragmas"
	#pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#ifndef BOOST_TEST_MODULE
	#define BOOST_TEST_MAIN
#endif

#include <boost/test/unit_test.hpp>

#if defined(__clang__)
	#pragma clang diagnostic pop
#elif defined(__GNUC__)
	#pragma GCC diagnostic pop
#endif

#include <boost/multi/array.hpp>  // for operator!=, implicit...
// IWYU pragma: no_include <boost/multi_array/subarray.hpp>        // for const_sub_array, sub_array

#include <boost/concept/assert.hpp>              // for BOOST_CONCEPT_ASSERT
#include <boost/concept_check.hpp>               // for Assignable, CopyCons...
#include <boost/multi_array.hpp>                 // for multi_array
#include <boost/multi_array/concept_checks.hpp>  // for ConstMultiArrayConcept

#include <cstddef>  // for ptrdiff_t
#include <vector>   // for vector

BOOST_AUTO_TEST_CASE(concepts_boost_array) {
	using BMA [[maybe_unused]] = boost::multi_array<int, 2>;  // maybe_unused for bug in nvcc 11.8

	BOOST_CONCEPT_ASSERT((boost::multi_array_concepts::ConstMultiArrayConcept<BMA, 2>));
	BOOST_CONCEPT_ASSERT((boost::multi_array_concepts::MutableMultiArrayConcept<BMA, 2>));
}

BOOST_AUTO_TEST_CASE(concepts_boost_array_1D) {
	using BMA = boost::multi_array<int, 1>;

	BOOST_CONCEPT_ASSERT((boost::multi_array_concepts::ConstMultiArrayConcept<BMA, 1>));
	BOOST_CONCEPT_ASSERT((boost::multi_array_concepts::MutableMultiArrayConcept<BMA, 1>));
}

namespace multi = boost::multi;

BOOST_AUTO_TEST_CASE(backwards) {
	multi::array<int, 2> const MA({ 2, 2 });

#ifdef __GNUC__
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
	// BOOST_REQUIRE(A.index_bases()[0] == 0);  // dangles?
	// BOOST_REQUIRE(A.index_bases()[1] == 0);

#ifdef __GNUC__
	#pragma GCC diagnostic pop
#endif

	{
#ifdef __NVCC__
	#pragma nv_diagnostic push
	#pragma nv_diag_suppress = deprecated_entity_with_custom_message  // nvcc #?
#endif

#ifdef __GNUC__
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
	// auto ib = MA.index_bases(); (void)ib;
	// BOOST_REQUIRE(ib[0] == 0);  // dangles?
	// BOOST_REQUIRE(ib[1] == 0);

#ifdef __GNUC__
	#pragma GCC diagnostic pop
#endif

#ifdef __NVCC__
	#pragma nv_diagnostic pop
#endif
	}
	// {
	//  #pragma GCC diagnostic push
	//  #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
	//  std::array<std::ptrdiff_t, 2> ib(MA.index_bases()); (void)ib;
	//  BOOST_REQUIRE(ib[0] == 0);
	//  BOOST_REQUIRE(ib[1] == 0);
	//  #pragma GCC diagnostic pop
	// }
	{
#ifdef __GNUC__
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
	// BOOST_REQUIRE(static_cast<std::ptrdiff_t const*>(MA.index_bases())[0] == 0);  // dangles
	// BOOST_REQUIRE(static_cast<std::ptrdiff_t const*>(MA.index_bases())[1] == 0);

#ifdef __GNUC__
	#pragma GCC diagnostic pop
#endif
	} {
#ifdef __GNUC__
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
	// BOOST_REQUIRE(MA.index_bases()[0] == 0);  // dangles
	// BOOST_REQUIRE(MA.index_bases()[1] == 0);

#ifdef __GNUC__
	#pragma GCC diagnostic pop
#endif
	} {
#ifdef __GNUC__
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
	// std::ptrdiff_t const* ib = MA.index_bases();
	// BOOST_REQUIRE(ib);
	// BOOST_REQUIRE(ib[0] == 0);  // dangles
	// BOOST_REQUIRE(ib[1] == 0);
#ifdef __GNUC__
	#pragma GCC diagnostic pop
#endif
	} {
#ifdef __GNUC__
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
		std::vector<std::ptrdiff_t> const ib(2);
		(void)ib;  // NOLINT(fuchsia-default-arguments-calls)
		           // std::copy_n(static_cast<std::ptrdiff_t const*>(MA.index_bases()), 2, ib.begin());
		           // BOOST_REQUIRE(ib[0] == 0);  // dangles
		           // BOOST_REQUIRE(ib[1] == 0);

#ifdef __GNUC__
	#pragma GCC diagnostic pop
#endif
	}
	// {
	//  #pragma GCC diagnostic push
	//  #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
	//  std::vector<std::ptrdiff_t> ib(2);  // NOLINT(fuchsia-default-arguments-calls)
	//  std::copy_n(MA.index_bases().to_array().data(), 2, ib.begin());
	//  BOOST_REQUIRE(ib[0] == 0);
	//  BOOST_REQUIRE(ib[1] == 0);
	//  #pragma GCC diagnostic pop
	// }
}

BOOST_AUTO_TEST_CASE(concepts_array) {
	using MA = multi::array<int, 2>;

	// BOOST_CONCEPT_ASSERT((boost::multi_array_concepts::ConstMultiArrayConcept<MA, 2>));
	// BOOST_CONCEPT_ASSERT((boost::multi_array_concepts::MutableMultiArrayConcept<MA, 2>));

	BOOST_CONCEPT_ASSERT((boost::Assignable<MA>));
	BOOST_CONCEPT_ASSERT((boost::SGIAssignable<MA>));
	BOOST_CONCEPT_ASSERT((boost::DefaultConstructible<MA>));
	BOOST_CONCEPT_ASSERT((boost::CopyConstructible<MA>));
	BOOST_CONCEPT_ASSERT((boost::EqualityComparable<MA>));
	BOOST_CONCEPT_ASSERT((boost::LessThanComparable<MA>));
	// BOOST_CONCEPT_ASSERT((boost::Comparable<MA>));

	// Function Object Concept Checking Classes
	BOOST_CONCEPT_ASSERT((boost::Generator<MA, boost::multi::subarray<int, 2>>));
	BOOST_CONCEPT_ASSERT((boost::UnaryFunction<MA, MA::reference, MA::index>));
	BOOST_CONCEPT_ASSERT((boost::BinaryFunction<MA, typename MA::element_ref, typename MA::index, typename MA::index>));

	// vvv needs result_type TODO(correaa) add to array?, should result_type be array<T, D-1>? or subarray?
	// BOOST_CONCEPT_ASSERT((boost::AdaptableGenerator<MA, boost::multi::subarray<int, 2>>));
	// BOOST_CONCEPT_ASSERT((boost::AdaptableUnaryFunction<MA, MA::reference, MA::index>));
	// BOOST_CONCEPT_ASSERT((boost::AdaptableBinaryFunction<MA, typename MA::element_ref, typename MA::index, typename MA::index>));

	// Container Concept Checking Classes
	BOOST_CONCEPT_ASSERT((boost::Container<MA>));
	BOOST_CONCEPT_ASSERT((boost::Mutable_Container<MA>));
	BOOST_CONCEPT_ASSERT((boost::ForwardContainer<MA>));
	BOOST_CONCEPT_ASSERT((boost::Mutable_ForwardContainer<MA>));
	// BOOST_CONCEPT_ASSERT((boost::ReversibleContainer<MA>));  // TODO(correaa) make it reversible, `const_reverse_iterator _i = cc.rbegin();`
	// BOOST_CONCEPT_ASSERT((boost::Mutable_ReversibleContainer<MA>));
	// BOOST_CONCEPT_ASSERT((boost::RandomAccessContainer<MA>));
	// BOOST_CONCEPT_ASSERT((boost::Mutable_RandomAccessContainer<MA>));
	// BOOST_CONCEPT_ASSERT((boost::Sequence<MA>));  // TODO(correaa) needs insert and erase, which will not be provided
	BOOST_CONCEPT_ASSERT((boost::Collection<MA>));
}

BOOST_AUTO_TEST_CASE(concepts_array_1D) {
	using MA = multi::array<int, 1>;

	BOOST_CONCEPT_ASSERT((boost::multi_array_concepts::ConstMultiArrayConcept<MA, 1>));
	// BOOST_CONCEPT_ASSERT((boost::multi_array_concepts::MutableMultiArrayConcept<MA, 2>));

	BOOST_CONCEPT_ASSERT((boost::Assignable<MA>));
	BOOST_CONCEPT_ASSERT((boost::SGIAssignable<MA>));
	BOOST_CONCEPT_ASSERT((boost::DefaultConstructible<MA>));
	BOOST_CONCEPT_ASSERT((boost::CopyConstructible<MA>));
	BOOST_CONCEPT_ASSERT((boost::EqualityComparable<MA>));
	BOOST_CONCEPT_ASSERT((boost::LessThanComparable<MA>));
	// BOOST_CONCEPT_ASSERT((boost::Comparable<MA>));

	// Function Object Concept Checking Classes
	BOOST_CONCEPT_ASSERT((boost::Generator<MA, boost::multi::subarray<int, 1>>));
	BOOST_CONCEPT_ASSERT((boost::UnaryFunction<MA, MA::reference, MA::index>));
	// BOOST_CONCEPT_ASSERT((boost::BinaryFunction<MA, typename MA::element_ref, typename MA::index, typename MA::index>));

	// vvv--- needs result_type TODO(correaa) add to array?, should result_type be array<T, D-1>? or subarray?
	// BOOST_CONCEPT_ASSERT((boost::AdaptableGenerator<MA, boost::multi::subarray<int, 2>>));
	// BOOST_CONCEPT_ASSERT((boost::AdaptableUnaryFunction<MA, MA::reference, MA::index>));
	// BOOST_CONCEPT_ASSERT((boost::AdaptableBinaryFunction<MA, typename MA::element_ref, typename MA::index, typename MA::index>));

	// Container Concept Checking Classes
	BOOST_CONCEPT_ASSERT((boost::Container<MA>));
	BOOST_CONCEPT_ASSERT((boost::Mutable_Container<MA>));
	BOOST_CONCEPT_ASSERT((boost::ForwardContainer<MA>));
	BOOST_CONCEPT_ASSERT((boost::Mutable_ForwardContainer<MA>));
	// BOOST_CONCEPT_ASSERT((boost::ReversibleContainer<MA>));  // TODO(correaa) make it reversible, `const_reverse_iterator _i = cc.rbegin();`
	// BOOST_CONCEPT_ASSERT((boost::Mutable_ReversibleContainer<MA>));
	// BOOST_CONCEPT_ASSERT((boost::RandomAccessContainer<MA>));
	// BOOST_CONCEPT_ASSERT((boost::Mutable_RandomAccessContainer<MA>));
	// BOOST_CONCEPT_ASSERT((boost::Sequence<MA>));  // TODO(correaa) needs insert and erase, which will not be provided
	BOOST_CONCEPT_ASSERT((boost::Collection<MA>));
}

BOOST_AUTO_TEST_CASE(concepts_iterator) {
	using MAIt = multi::array<int, 2>::iterator;

	BOOST_CONCEPT_ASSERT((boost::Assignable<MAIt>));
	BOOST_CONCEPT_ASSERT((boost::SGIAssignable<MAIt>));
	BOOST_CONCEPT_ASSERT((boost::DefaultConstructible<MAIt>));
	BOOST_CONCEPT_ASSERT((boost::CopyConstructible<MAIt>));
	BOOST_CONCEPT_ASSERT((boost::EqualityComparable<MAIt>));
	BOOST_CONCEPT_ASSERT((boost::LessThanComparable<MAIt>));

	BOOST_CONCEPT_ASSERT((boost::InputIterator<MAIt>));
	BOOST_CONCEPT_ASSERT((boost::OutputIterator<MAIt, MAIt::reference>));
	BOOST_CONCEPT_ASSERT((boost::OutputIterator<MAIt, MAIt::value_type>));

	// Iterator Concept Checking Classes
	BOOST_CONCEPT_ASSERT((boost::ForwardIterator<MAIt>));
	BOOST_CONCEPT_ASSERT((boost::Mutable_ForwardIterator<MAIt>));
	BOOST_CONCEPT_ASSERT((boost::BidirectionalIterator<MAIt>));
	BOOST_CONCEPT_ASSERT((boost::Mutable_BidirectionalIterator<MAIt>));
	BOOST_CONCEPT_ASSERT((boost::RandomAccessIterator<MAIt>));
	BOOST_CONCEPT_ASSERT((boost::Mutable_RandomAccessIterator<MAIt>));
}

BOOST_AUTO_TEST_CASE(concepts_const_iterator) {
	using MAIt [[maybe_unused]] = multi::array<int, 2>::const_iterator;  // maybe_unused for bug in nvcc 11.8

	BOOST_CONCEPT_ASSERT((boost::Assignable<MAIt>));
	BOOST_CONCEPT_ASSERT((boost::SGIAssignable<MAIt>));
	BOOST_CONCEPT_ASSERT((boost::DefaultConstructible<MAIt>));
	BOOST_CONCEPT_ASSERT((boost::CopyConstructible<MAIt>));
	BOOST_CONCEPT_ASSERT((boost::EqualityComparable<MAIt>));
	BOOST_CONCEPT_ASSERT((boost::LessThanComparable<MAIt>));

	BOOST_CONCEPT_ASSERT((boost::InputIterator<MAIt>));
	// BOOST_CONCEPT_ASSERT((boost::OutputIterator<MAIt, MAIt::reference>));
	//   BOOST_CONCEPT_ASSERT((boost::OutputIterator<MAIt, MAIt::value_type>));

	BOOST_CONCEPT_ASSERT((boost::ForwardIterator<MAIt>));
	//  BOOST_CONCEPT_ASSERT((boost::Mutable_ForwardIterator<MAIt>));
	BOOST_CONCEPT_ASSERT((boost::BidirectionalIterator<MAIt>));
	//  BOOST_CONCEPT_ASSERT((boost::Mutable_BidirectionalIterator<MAIt>));
	BOOST_CONCEPT_ASSERT((boost::RandomAccessIterator<MAIt>));
	//  BOOST_CONCEPT_ASSERT((boost::Mutable_RandomAccessIterator<MAIt>));
}

#if defined(__clang__)
	#pragma clang diagnostic pop
#elif defined(__GNUC__)
	#pragma GCC diagnostic pop
#endif
