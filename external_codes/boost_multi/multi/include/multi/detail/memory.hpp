// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2022 Alfredo A. Correa

#ifndef MULTI_DETAIL_MEMORY_HPP
#define MULTI_DETAIL_MEMORY_HPP

#include <memory>  // for std::allocator_traits

namespace boost::multi {

template<class Alloc>
struct allocator_traits : std::allocator_traits<Alloc> {
#if 0
	template<class Ptr, class... Args>
	static auto construct(Alloc& alloc, Ptr p, Args&&... args)  // NOLINT(readability-identifier-length) std naming
	->decltype(alloc.construct(p, std::forward<Args>(args)...)) {
		return alloc.construct(p, std::forward<Args>(args)...); }

	template<class Ptr>
	static auto destroy(Alloc& alloc, Ptr p)  // NOLINT(readability-identifier-length) std naming
	->decltype(alloc.destroy(p)) {
		return alloc.destroy(p); }
#endif
};

// https://en.cppreference.com/w/cpp/memory/destroy
template<class Alloc, class ForwardIt, std::enable_if_t<!has_rank<ForwardIt>::value, int> = 0>
void destroy(Alloc& alloc, ForwardIt first, ForwardIt last) {
	for(; first != last; ++first) {allocator_traits<Alloc>::destroy(alloc, std::addressof(*first));}  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
}

template<class Alloc, class ForwardIt, std::enable_if_t<has_rank<ForwardIt>::value and ForwardIt::rank_v == 1, int> = 0>
void destroy(Alloc& alloc, ForwardIt first, ForwardIt last) {
	//  using multi::to_address;
	for(; first != last; ++first) {alloc.destroy(to_address(first));}  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
}

template<class Alloc, class ForwardIt, std::enable_if_t<has_rank<ForwardIt>::value and ForwardIt::rank_v != 1, int> = 0>
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
	typename = std::enable_if_t<not std::is_trivially_default_constructible<T>{}>
>
auto uninitialized_value_construct_n(Alloc& alloc, ForwardIt first, Size n) -> ForwardIt {
	ForwardIt current = first;  // using std::addressof;
	try {
		for(; n > 0; ++current, --n) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			allocator_traits<Alloc>::construct(alloc, to_address(current), T{});
		}
		//  a.construct(to_address(current), T());  //	a.construct(std::pointer_traits<Ptr>::pointer_to(*current), T());  //	AT::construct(a, to_address(current), T());  //	AT::construct(a, addressof(*current), T()); //	a.construct(addressof(*current), T());
		return current;
	} catch(...) {destroy(alloc, first, current); throw;}
}

template<class... Args> auto std_copy(Args&&... args) {
	using std::copy;
	return copy(std::forward<Args>(args)...);
}

namespace xtd {

template<class Alloc, class InputIt, class MIt, typename = std::enable_if_t<not has_rank<MIt>{}> >
auto alloc_uninitialized_copy(Alloc& alloc, InputIt first, InputIt last, MIt dest) -> MIt {
	MIt current = dest;
//  using multi::to_address;
	try {
		for(; first != last; ++first, ++current) {alloc.construct(to_address(current), *first);}  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
		return current;
	} catch(...) {destroy(alloc, dest, current); throw;}
}

}  // end namespace xtd

// // https://en.cppreference.com/w/cpp/memory/destroy_at
// template<class Alloc, class T, typename AT = std::allocator_traits<Alloc> >
// void destroy_at(Alloc& a, T* p) {AT::destroy(a, p);}

// // https://en.cppreference.com/w/cpp/memory/destroy_n
// template<class Alloc, class ForwardIt, class Size>  // , typename AT = typename std::allocator_traits<Alloc> >
// auto destroy_n(Alloc& a, ForwardIt first, Size n) -> ForwardIt {
// //  using std::addressof;
// 	for(; n > 0; ++first, --n) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
// 		allocator_traits<Alloc>::destroy(a, to_address(first));
// 	}
// 	return first;
// }

template<class AA> class is_allocator {
	template<
		class A,
		class P = typename A::pointer, class S = typename A::size_type,
		typename = decltype(
			std::declval<A const&>() == A{std::declval<A const&>()},
			std::declval<A&>().deallocate(P{std::declval<A&>().allocate(std::declval<S>())}, std::declval<S>())
		)
	>
	static auto  aux(A const&) -> std::true_type;
	static auto  aux(...     ) -> std::false_type;

 public:
	static bool const value = decltype(aux(std::declval<AA>()))::value;
	constexpr explicit operator bool() const {return value;}
};

template<class Alloc> constexpr bool is_allocator_v = is_allocator<Alloc>::value;

template<dimensionality_type N, class InputIt, class ForwardIt>
auto uninitialized_copy(InputIt first, InputIt last, ForwardIt dest) {
	while(first!=last) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
		uninitialized_copy<N-1>(begin(*first), end(*first), begin(*dest));
		++first;
		++dest;
	}
	return dest;
}

}  // end namespace boost::multi

#endif
