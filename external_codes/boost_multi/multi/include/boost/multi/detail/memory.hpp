// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_DETAIL_MEMORY_HPP
#define BOOST_MULTI_DETAIL_MEMORY_HPP
#pragma once

#include <boost/multi/utility.hpp>  // for has_rank, to_address

#include <iterator>                 // for copy, iterator_traits
#include <memory>                   // for allocator_traits
#include <type_traits>              // for declval, enable_if_t, false_type, is_trivially_default_constructible, true_type, void_t
#include <utility>                  // for addressof, forward

namespace boost::multi {

template<class Alloc>
struct allocator_traits : std::allocator_traits<Alloc> {};

// https://en.cppreference.com/w/cpp/memory/destroy
template<
	class Alloc, class ForwardIt,
	std::enable_if_t<!has_rank<ForwardIt>::value, int> =0  // NOLINT(modernize-use-constraints) TODO(correaa)
>
void destroy(Alloc& alloc, ForwardIt first, ForwardIt last) {
	for(; first != last; ++first) {allocator_traits<Alloc>::destroy(alloc, std::addressof(*first));}  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
}

template<class Alloc, class ForwardIt, std::enable_if_t<has_rank<ForwardIt>::value && ForwardIt::rank_v == 1, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
void destroy(Alloc& alloc, ForwardIt first, ForwardIt last) {
	//  using multi::to_address;
	std::for_each(first, last, [&](auto& elem) {alloc.destroy(addressof(elem));});
	// for(; first != last; ++first) {alloc.destroy(to_address(first));}  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
}

template<class Alloc, class ForwardIt, std::enable_if_t<has_rank<ForwardIt>::value && ForwardIt::rank_v != 1, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
void destroy(Alloc& alloc, ForwardIt first, ForwardIt last) {
	for(; first != last; ++first) {destroy(alloc, begin(*first), end(*first));} // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
}

template<class Alloc, class InputIt, class Size, class ForwardIt>
auto uninitialized_move_n(Alloc& alloc, InputIt first, Size n, ForwardIt dest) -> ForwardIt {
	ForwardIt curr = dest;
//  using std::addressof;
	try {
		for(; n > 0; ++first, ++curr, --n) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			alloc.construct(std::addressof(*curr), std::move(*first));
		}
		return curr;
	} catch(...) {destroy(alloc, dest, curr); throw;}
}

template<class Alloc, class ForwardIt, class Size>
auto uninitialized_default_construct_n(Alloc& alloc, ForwardIt first, Size n) -> ForwardIt {
	ForwardIt current = first;
	try {
		for(; n > 0; ++current, --n) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
		//  allocator_traits<Alloc>::construct(a, to_address(current));
			alloc.construct(to_address(current));
		}
		return current;
	} catch(...) {destroy(alloc, first, current); throw;}
}

template<
	class Alloc, class ForwardIt, class Size,
	typename T = typename std::iterator_traits<ForwardIt>::value_type,
	typename = std::enable_if_t<! std::is_trivially_default_constructible<T>{}>  // NOLINT(modernize-use-constraints) TODO(correaa)
>
auto uninitialized_value_construct_n(Alloc& alloc, ForwardIt first, Size n) -> ForwardIt {
	ForwardIt current = first;  // using std::addressof;
	try {
		for(; n > 0; ++current, --n) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			allocator_traits<Alloc>::construct(alloc, to_address(current), T{});
		}
		//  a.construct(to_address(current), T());  //  a.construct(std::pointer_traits<Ptr>::pointer_to(*current), T());  //   AT::construct(a, to_address(current), T());  // AT::construct(a, addressof(*current), T()); //  a.construct(addressof(*current), T());
		return current;
	} catch(...) {destroy(alloc, first, current); throw;}
}

template<class... Args> auto std_copy(Args&&... args) {
	using std::copy;
	return copy(std::forward<Args>(args)...);
}

namespace xtd {

template<class Alloc, class InputIt, class MIt, typename = std::enable_if_t<! has_rank<MIt>{}> >  // NOLINT(modernize-use-constraints) TODO(correaa)
auto alloc_uninitialized_copy(Alloc& alloc, InputIt first, InputIt last, MIt dest) -> MIt {
	MIt current = dest;
//  using multi::to_address;
	try {
		for(; first != last; ++first, ++current) {alloc.construct(to_address(current), *first);}  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
		return current;
	} catch(...) {destroy(alloc, dest, current); throw;}
}

}  // end namespace xtd

template<class, class = void>
struct is_allocator : std::false_type {};

template<class Alloc>
struct is_allocator<Alloc, std::void_t<decltype(
	std::declval<Alloc const&>() == Alloc{std::declval<Alloc const&>()},
	std::declval<Alloc&>().deallocate(typename Alloc::pointer{std::declval<Alloc&>().allocate(std::declval<typename Alloc::size_type>())}, std::declval<typename Alloc::size_type>())
)>> : std::true_type {};

template<class Alloc> constexpr bool is_allocator_v = is_allocator<Alloc>::value;

// template<dimensionality_type N, class InputIt, class ForwardIt>
// auto uninitialized_copy(InputIt first, InputIt last, ForwardIt dest) {
//  while(first!=last) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
//    uninitialized_copy<N-1>(begin(*first), end(*first), begin(*dest));
//    ++first;
//    ++dest;
//  }
//  return dest;
// }

}  // end namespace boost::multi
#endif  // BOOST_MULTI_DETAIL_MEMORY_HPP
