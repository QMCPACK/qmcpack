// Copyright 2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_DETAIL_INTERFACES_HPP
#define BOOST_MULTI_DETAIL_INTERFACES_HPP
// #pragma once

#include <cstddef>      // for ptrdiff_t
#include <iterator>     // for random_access_iterator_tag
#include <type_traits>  // for enable_if_t, is_base_of
#include <utility>      // for forward

#ifdef __NVCC__
	#define BOOST_MULTI_HD __host__ __device__
#else
	#define BOOST_MULTI_HD
#endif

namespace boost::multi::detail {

template<class Self> 
class equality_comparable {
	using self_type = Self;
	friend auto operator!=(self_type const& self, self_type const& other) -> bool {
		return !(self == other);
	}

 protected:
	equality_comparable() = default;  // NOLINT(bugprone-crtp-constructor-accessibility)
};

template<class Self>
class regular : equality_comparable<Self> {
 protected:
	regular() = default;  // NOLINT(bugprone-crtp-constructor-accessibility)
};

template<class Self>
class weakly_incrementable_facade {
	using self_type = Self;
 protected:
	weakly_incrementable_facade() = default;  // NOLINT(bugprone-crtp-constructor-accessibility)
	// using difference_type = typename self_type::difference_type;

 public:
	template<int = 0>
	friend auto operator++(self_type& self, int) { ++self; }
};

template<class Self>
class incrementable_facade : public regular<Self>, public weakly_incrementable_facade<Self> {
	using self_type = Self;
 protected:
	incrementable_facade() = default;  // NOLINT(bugprone-crtp-constructor-accessibility)
	
 public:
	friend auto operator++(self_type& self, int) {  // cppcheck-suppress duplInheritedMember;
		auto ret{self}; ++self; return ret;
	}
};

template< class I >
class input_or_output_iterator_facade :
        // requires(I i) {
        //     { *i } -> /*can-reference*/;
        // } &&
        weakly_incrementable_facade<I>
{
 protected:
	input_or_output_iterator_facade() = default;  // NOLINT(bugprone-crtp-constructor-accessibility)
};


template< class I >
class input_iterator_facade
	: input_or_output_iterator_facade<I>
        // std::input_or_output_iterator<I> &&
        // std::indirectly_readable<I> &&
        // requires { typename /*ITER_CONCEPT*/<I>; } &&
        // std::derived_from</*ITER_CONCEPT*/<I>, std::input_iterator_tag>;
{
 protected:
	input_iterator_facade() = default;  // NOLINT(bugprone-crtp-constructor-accessibility)
};

template< class I >
class forward_iterator_facade : 
	public input_iterator_facade<I>,
	// std::derived_from</*ITER_CONCEPT*/<I>, std::forward_iterator_tag> &&
    public incrementable_facade<I>  // std::incrementable<I> &&
    //    std::sentinel_for<I, I>;
{
	protected:
	forward_iterator_facade() = default;  // NOLINT(bugprone-crtp-constructor-accessibility)
};

}  // end namespace boost::multi::detail

#undef BOOST_MULTI_HD

#endif  // BOOST_MULTI_DETAIL_INTERFACES_HPP
