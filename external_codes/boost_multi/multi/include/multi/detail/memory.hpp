// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2021 Alfredo A. Correa

#ifndef BOOST_MULTI_DETAIL_MEMORY_HPP
#define BOOST_MULTI_DETAIL_MEMORY_HPP

#include "../detail/adl.hpp"
#include "../utility.hpp"

#include "layout.hpp"

#include <algorithm>  // for copy_n
#include <memory>
#include <utility>

namespace boost::multi {
namespace memory {

template<class Alloc>
struct allocator_traits : std::allocator_traits<Alloc> {
	template<class Ptr, class... Args>
	static auto construct(Alloc& a, Ptr p, Args&&... args)
	->decltype(a.construct(p, std::forward<Args>(args)...)) {
		return a.construct(p, std::forward<Args>(args)...); }

	template<class Ptr>
	static auto destroy(Alloc& a, Ptr p)
	->decltype(a.destroy(p)) {
		return a.destroy(p); }
};

}  // end namespace memory

using memory::allocator_traits;

// https://en.cppreference.com/w/cpp/memory/destroy
template<class Alloc, class ForwardIt, std::enable_if_t<not has_rank<ForwardIt>{}, int> = 0>
void destroy(Alloc& a, ForwardIt f, ForwardIt l) {
	for(; f != l; ++f) {allocator_traits<Alloc>::destroy(a, std::addressof(*f));}
}

template<class Alloc, class ForwardIt, std::enable_if_t<has_rank<ForwardIt>{} and ForwardIt::rank_v == 1, int> = 0>
void destroy(Alloc& a, ForwardIt first, ForwardIt last) {
	//  using multi::to_address;
	for(; first != last; ++first) {a.destroy(to_address(first));}
}

template<class Alloc, class ForwardIt, std::enable_if_t<has_rank<ForwardIt>{} and ForwardIt::rank_v != 1, int> = 0>
void destroy(Alloc& a, ForwardIt first, ForwardIt last) {
	for(; first != last; ++first) {destroy(a, begin(*first), end(*first));}
}

template<class Alloc, class InputIt, class Size, class ForwardIt>
auto uninitialized_move_n(Alloc& a, InputIt f, Size n, ForwardIt d) -> ForwardIt {
	ForwardIt c = d;
//  using std::addressof;
	try {
		for(; n > 0; ++f, ++c, --n) {
			a.construct(std::addressof(*c), std::move(*f));
		}
		return c;
	} catch(...) {destroy(a, d, c); throw;}
}

template<class Alloc, class ForwardIt, class Size>
auto uninitialized_default_construct_n(Alloc& a, ForwardIt first, Size n) -> ForwardIt {
	ForwardIt current = first;
	try {
		for(; n > 0; ++current, --n) {
		//  allocator_traits<Alloc>::construct(a, to_address(current));
			a.construct(to_address(current));
		}
		return current;
	} catch(...) {destroy(a, first, current); throw;}
}

template<
	class Alloc, class ForwardIt, class Size,
	typename T = typename std::iterator_traits<ForwardIt>::value_type,
	typename = std::enable_if_t<not std::is_trivially_default_constructible<T>{}>
>
auto uninitialized_value_construct_n(Alloc& a, ForwardIt first, Size n) -> ForwardIt {
	ForwardIt current = first;  // using std::addressof;
	try {
		for(; n > 0; ++current, --n) {
			allocator_traits<Alloc>::construct(a, to_address(current), T());
		}
		//  a.construct(to_address(current), T());  //	a.construct(std::pointer_traits<Ptr>::pointer_to(*current), T());  //	AT::construct(a, to_address(current), T());  //	AT::construct(a, addressof(*current), T()); //	a.construct(addressof(*current), T());
		return current;
	} catch(...) {destroy(a, first, current); throw;}
}

template<class... Args> auto std_copy(Args&&... args) {
	using std::copy;
	return copy(std::forward<Args>(args)...);
}

namespace xtd {

template<class Alloc, class InputIt, class MIt, typename = std::enable_if_t<not has_rank<MIt>{}> >
auto alloc_uninitialized_copy(Alloc& a, InputIt f, InputIt l, MIt d) -> MIt {
	MIt current = d;
//  using multi::to_address;
	try {
		for(; f != l; ++f, ++current) {a.construct(to_address(current), *f);}
		return current;
	} catch(...) {destroy(a, d, current); throw;}
}

}  // end namespace xtd

// https://en.cppreference.com/w/cpp/memory/destroy_at
template<class Alloc, class T, typename AT = std::allocator_traits<Alloc> >
void destroy_at(Alloc& a, T* p) {AT::destroy(a, p);}

// https://en.cppreference.com/w/cpp/memory/destroy_n
template<class Alloc, class ForwardIt, class Size>  // , typename AT = typename std::allocator_traits<Alloc> >
auto destroy_n(Alloc& a, ForwardIt first, Size n) -> ForwardIt {
//  using std::addressof;
	for(; n > 0; ++first, --n) {
		allocator_traits<Alloc>::destroy(a, to_address(first));
	}
	return first;
}

template<class AA> class is_allocator {
	template<
		class A,
		class P = typename A::pointer, class S = typename A::size_type,
		typename = decltype(
			std::declval<A const&>()==A{std::declval<A const&>()},
			std::declval<A&>().deallocate(P{std::declval<A&>().allocate(std::declval<S>())}, std::declval<S>())
		)
	>
	static auto  aux(A const&) -> std::true_type;
	static auto  aux(...     ) -> std::false_type;

 public:
	static bool const value = decltype(aux(std::declval<AA>()))::value;
	constexpr explicit operator bool() const{return value;}
};

template<dimensionality_type N, class InputIt, class ForwardIt>
auto uninitialized_copy(InputIt first, InputIt last, ForwardIt dest) {
	while(first!=last) {
		uninitialized_copy<N-1>(begin(*first), end(*first), begin(*dest));
		++first;
		++dest;
	}
	return dest;
}

template<dimensionality_type N> struct recursive_fill_aux;

template<dimensionality_type N, class Out, class T>
void recursive_fill(Out f, Out l, T const& value) {
	return recursive_fill_aux<N>::call(f, l, value);
}

template<dimensionality_type N>
struct recursive_fill_aux {
	template<class Out, class T>
	static auto call(Out first, Out last, T const& value) {
		using std::begin; using std::end;
		for(; first != last; ++first) {
			recursive_fill<N-1>(begin(*first), end(*first), value);  // (*first).begin() instead of first->begin() to make it work with T[][]
		}
	}
};

template<> struct recursive_fill_aux<1> {
	template<class O, class T>  static auto call(O f, O l, T const& v) {
		using std::fill; return fill(f, l, v);
	}
};

#if 0
template<dimensionality_type N> struct recursive_uninitialized_fill_aux;

template<dimensionality_type N, class Alloc, class Out, class T>
void recursive_uninitialized_fill(Alloc& a, Out f, Out l, T const& value) {
	return recursive_uninitialized_fill_aux<N>::call(a, f, l, value);
}

template<dimensionality_type N>
struct uninitialized_fill_aux{
	template<class Alloc, class Out, class T>
	static auto call(Alloc& a, Out first, Out last, T const& v){
		using std::begin; using std::end;
		for(; first != last; ++first)
			recursive_uninitialized_fill<N-1>(a, begin(*first), end(*first), v);  // (*first).begin() instead of first->begin() to make it work with T[][]
	}
};

template<> struct uninitialized_fill_aux<1>{template<class Alloc, class O, class T>
	static auto call(Alloc& a, O f, O l, T const& v){return uninitialized_fill(a, f, l, v);}
};
#endif

}  // end namespace boost::multi

//#if defined(__INCLUDE_LEVEL__) and not __INCLUDE_LEVEL__

//#include<vector>

//namespace multi = boost::multi;

//template<class T> void what(T&&);

//int main() {
//	static_assert(multi::is_allocator<std::allocator<double>>{}, "!");
//	static_assert(not multi::is_allocator<double>{}, "!");
//	static_assert(not multi::is_allocator<std::vector<double>>{}, "!");

//	 {
//		double* p = nullptr;
//		auto a = multi::default_allocator_of(p);
//		static_assert(std::is_same<decltype(a), std::allocator<double>>{}, "!");
//	//  what(typename std::iterator_traits<double*>::value_type{});
//	}
//#if 0
//	 {
//		std::vector<double>::iterator it;
//		auto a = multi::default_allocator_of(it);
//	//  what(typename std::iterator_traits<decltype(it)>::value_type{});
//		static_assert(std::is_same<decltype(a), std::allocator<double>>{}, "!");
//	}
//#endif
//}
//#endif
#endif
