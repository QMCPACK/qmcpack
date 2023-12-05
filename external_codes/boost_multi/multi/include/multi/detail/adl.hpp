// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2020-2022 Alfredo A. Correa

#ifndef MULTI_DETAIL_ADL_HPP
#define MULTI_DETAIL_ADL_HPP

#include<cstddef>      // std::size_t
#include<type_traits>  // std::conditional_t
#include<utility>

#include "multi/detail/memory.hpp"

#if defined(__NVCC__)
#include<thrust/copy.h>
#include<thrust/equal.h>
#include<thrust/detail/memory_algorithms.h>
#include<thrust/uninitialized_copy.h>
#endif

#include<algorithm>  // for std::copy, std::copy_n, std::equal, etc
#include<iterator>   // for begin, end
#include<memory>     // for uninitialized_copy, etc

#define BOOST_MULTI_DEFINE_ADL(FuN)  /*NOLINT(cppcoreguidelines-macro-usage) TODO(correaa) consider replacing for all ADL'd operations*/ \
namespace boost { \
namespace multi { \
namespace adl { \
	namespace custom {template<class...> struct FuN##_t;} 	__attribute__((unused))  \
	static constexpr class FuN##_t { \
		template<class... As> [[deprecated]] auto _(priority<0>,        As&&... args) const = delete; \
		template<class... As>          auto _(priority<1>,        As&&... args) const DECLRETURN(std::FuN(std::forward<As>(args)...)) \
		template<class... As>          auto _(priority<2>,        As&&... args) const DECLRETURN(     FuN(std::forward<As>(args)...)) \
		template<class T, class... As> auto _(priority<3>, T&& t, As&&... args) const DECLRETURN(std::forward<T>(t).FuN(std::forward<As>(args)...))     \
		template<class... As>          auto _(priority<4>,        As&&... args) const DECLRETURN(custom::FuN##_t<As&&...>::_(std::forward<As>(args)...)) \
	public: \
		template<class... As> auto operator()(As&&... args) const-> decltype(_(priority<4>{}, std::forward<As>(args)...)) {return _(priority<4>{}, std::forward<As>(args)...);} \
	} (FuN); \
}  /* end namespace adl   */ \
}  /* end namespace multi */ \
}  /* end namespace boost */

namespace boost::multi {

template<std::size_t I> struct priority : std::conditional_t<I==0, std::true_type, priority<I-1>> {};

#define DECLRETURN(ExpR) ->decltype(ExpR) {return ExpR;}  // NOLINT(cppcoreguidelines-macro-usage) saves a lot of typing
#define JUSTRETURN(ExpR)                  {return ExpR;}  // NOLINT(cppcoreguidelines-macro-usage) saves a lot of typing

constexpr class adl_copy_n_t {
	template<class... As>          constexpr auto _(priority<0>/**/,          As&&... args) const DECLRETURN(std::                copy_n(                      std::forward<As>(args)...))
#if defined(__NVCC__)
	template<class... As> 		   constexpr auto _(priority<1>/**/,          As&&... args) const DECLRETURN(::thrust::           copy_n(                      std::forward<As>(args)...))
#endif
	template<class... As>          constexpr auto _(priority<2>/**/,          As&&... args) const DECLRETURN(                     copy_n(                      std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& arg, As&&... args) const DECLRETURN(std::decay_t<T>::    copy_n(std::forward<T>(arg), std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<4>/**/, T&& arg, As&&... args) const DECLRETURN(std::forward<T>(arg).copy_n(                      std::forward<As>(args)...))

 public:
	template<class... As> constexpr auto operator()(As&&... args) const DECLRETURN(_(priority<4>{}, std::forward<As>(args)...))
} adl_copy_n;

constexpr class adl_move_t {
	template<class... As>           constexpr auto _(priority<0>/**/,                      As&&... args) const DECLRETURN(              std::    move(                      std::forward<As>(args)...))
#if defined(__NVCC__)  // there is no thrust::move algorithm
	template<class It, class... As> constexpr auto _(priority<1>/**/, It first, It last, As&&... args) const DECLRETURN(           thrust::copy(std::make_move_iterator(first), std::make_move_iterator(last), std::forward<As>(args)...))
#endif
	template<class... As>           constexpr auto _(priority<2>/**/,                      As&&... args) const DECLRETURN(                     move(                      std::forward<As>(args)...))
	template<class T, class... As>  constexpr auto _(priority<3>/**/, T&& arg,             As&&... args) const DECLRETURN(std::decay_t<T>::    move(std::forward<T>(arg), std::forward<As>(args)...))
	template<class T, class... As>  constexpr auto _(priority<4>/**/, T&& arg,             As&&... args) const DECLRETURN(std::forward<T>(arg).move(                      std::forward<As>(args)...))

 public:
	template<class... As> constexpr auto operator()(As&&... args) const DECLRETURN(_(priority<4>{}, std::forward<As>(args)...))
} adl_move;

constexpr class adl_fill_n_t {
	template<         class... As> constexpr auto _(priority<0>/**/,          As&&... args) const DECLRETURN(              std::  fill_n              (std::forward<As>(args)...))
#if defined(__NVCC__)
	template<         class... As> constexpr auto _(priority<1>/**/,          As&&... args) const DECLRETURN(           thrust::  fill_n              (std::forward<As>(args)...))
#endif
	template<         class... As> constexpr auto _(priority<2>/**/,          As&&... args) const DECLRETURN(                     fill_n              (std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& arg, As&&... args) const DECLRETURN(std::decay_t<T>::    fill_n(std::forward<T>(arg), std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<4>/**/, T&& arg, As&&... args) const DECLRETURN(std::forward<T>(arg).fill_n              (std::forward<As>(args)...))

 public:
	template<class... As> constexpr auto operator()(As&&... args) const DECLRETURN(_(priority<4>{}, std::forward<As>(args)...))
} adl_fill_n;

constexpr class adl_equal_t {
	template<         class...As> /*[[gnu::pure]]*/ constexpr auto _(priority<1>/**/,          As&&...args) const DECLRETURN(               std::  equal(                      std::forward<As>(args)...))
#if defined(__NVCC__)
	template<         class...As>               constexpr auto _(priority<2>/**/,          As&&...args) const DECLRETURN(          ::thrust::  equal(                      std::forward<As>(args)...))
#endif
	template<         class...As> /*[[gnu::pure]]*/ constexpr auto _(priority<3>/**/,          As&&...args) const DECLRETURN(                      equal(                      std::forward<As>(args)...))
	template<class T, class...As>               constexpr auto _(priority<4>/**/, T&& arg, As&&...args) const DECLRETURN( std::decay_t<T>::    equal(std::forward<T>(arg), std::forward<As>(args)...))
	template<class T, class...As>               constexpr auto _(priority<5>/**/, T&& arg, As&&...args) const DECLRETURN( std::forward<T>(arg).equal(                      std::forward<As>(args)...))

 public:
	template<class...As> /*[[gnu::pure]]*/          constexpr auto operator()(As&&...args) const DECLRETURN(_(priority<5>{}, std::forward<As>(args)...))
} adl_equal;

template<class... Args> struct adl_custom_copy;

constexpr class adl_copy_t {
	class Copy;
	template<class InputIt, class OutputIt,
		class=std::enable_if_t<std::is_assignable_v<typename std::iterator_traits<OutputIt>::reference, typename std::iterator_traits<InputIt>::reference>>
	>
	                               constexpr auto _(priority<1>/**/, InputIt first, InputIt last, OutputIt d_first)	const DECLRETURN(std::copy(first, last, d_first))
#if defined(__NVCC__)
	template<class... As> 		   constexpr auto _(priority<2>/**/,          As&&... args) const DECLRETURN(           thrust::copy(std::forward<As>(args)...))
#endif
	template<         class... As> constexpr auto _(priority<3>/**/,          As&&... args) const DECLRETURN(                   copy(std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<4>/**/, T&& arg, As&&... args) const DECLRETURN(  std::decay_t<T>::copy(std::forward<T>(arg), std::forward<As>(args)...))
//	template<class... As         > constexpr auto _(priority<5>/**/,          As&&... args) const DECLRETURN(boost::multi::adl_custom_copy<std::decay_t<As>...>::copy(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<6>/**/, T&& arg, As&&... args) const DECLRETURN(std::forward<T>(arg).copy(std::forward<As>(args)...))

 public:
	template<class... As> constexpr auto operator()(As&&... args) const DECLRETURN( _(priority<6>{}, std::forward<As>(args)...) ) \
} adl_copy;

namespace adl {
	namespace custom {template<class...> struct fill_t;}
	static constexpr class fill_t {
		template<class... As>          auto _(priority<1>/**/,          As&&... args) const DECLRETURN(              std::  fill              (std::forward<As>(args)...))
		template<class... As>          auto _(priority<2>/**/,          As&&... args) const DECLRETURN(                     fill              (std::forward<As>(args)...))
		template<class T, class... As> auto _(priority<3>/**/, T&& arg, As&&... args) const DECLRETURN(std::forward<T>(arg).fill              (std::forward<As>(args)...))
		template<class... As>          auto _(priority<4>/**/,          As&&... args) const DECLRETURN(custom::             fill_t<As&&...>::_(std::forward<As>(args)...))
	
	 public:
		template<class... As> auto operator()(As&&... args) const DECLRETURN(_(priority<5>{}, std::forward<As>(args)...))
	} fill [[maybe_unused]];
}  // end namespace adl

template<class Alloc>
struct alloc_construct_elem_t {
	Alloc* palloc_;
	template<class T> auto operator()(T&& ptr) const
	->decltype(std::allocator_traits<Alloc>::construct(*palloc_, std::addressof(ptr))) {
		return std::allocator_traits<Alloc>::construct(*palloc_, std::addressof(ptr)); }
};

namespace xtd {

template<class T>  // this one goes last!!!
constexpr auto to_address(const T& ptr) noexcept;

template<class T>
constexpr auto me_to_address(priority<0>/**/, const T& ptr) noexcept
->decltype(to_address(ptr.operator->())) {
	return to_address(ptr.operator->()); }

template<class T>
constexpr auto me_to_address(priority<1>/**/, const T& ptr) noexcept
->decltype(std::pointer_traits<T>::to_address(ptr)) {
	return std::pointer_traits<T>::to_address(ptr); }

template<class T, std::enable_if_t<std::is_pointer<T>{}, int> = 0>
constexpr auto me_to_address(priority<2>/**/, T const& ptr) noexcept -> T {
    static_assert(not std::is_function_v<T>, "!");
    return ptr;
}

template<class T>  // this one goes last!!!
constexpr auto to_address(T const& ptr) noexcept
->decltype(me_to_address(priority<2>{}/**/, ptr)) {
	return me_to_address(priority<2>{}    , ptr); }

template<class Alloc, class ForwardIt, class Size, typename Value = typename std::iterator_traits<ForwardIt>::value_type, typename = decltype(std::addressof(*ForwardIt{})), typename = decltype(Value())>
auto alloc_uninitialized_value_construct_n(Alloc& alloc, ForwardIt first, Size count) -> ForwardIt {
// ->std::decay_t<decltype(std::allocator_traits<Alloc>::construct(alloc, std::addressof(*first), Value()), first)>
	ForwardIt current = first;
	try {
		for (; count > 0 ; ++current, --count) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			std::allocator_traits<Alloc>::construct(alloc, std::addressof(*current), Value());  // !!!!!!!!!!!!!! if you are using std::complex type consider making complex default constructible (e.g. by type traits)
		}
		//  ::new (static_cast<void*>(std::addressof(*current))) Value();
		return current;
	} catch(...) {
		for(; current != first; ++first) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			std::allocator_traits<Alloc>::destroy(alloc, std::addressof(*first));
		}
		throw;
	}
}

template<class Alloc, class ForwardIt, class Size, class T = typename std::iterator_traits<ForwardIt>::value_type>
auto alloc_uninitialized_default_construct_n(Alloc& alloc, ForwardIt first, Size count)
->std::decay_t<decltype(std::allocator_traits<Alloc>::construct(alloc, std::addressof(*first)), first)> {
	ForwardIt current = first;
	if constexpr(std::is_trivially_default_constructible_v<T>) {
		std::advance(current, count);
	} else {
		using _ = std::allocator_traits<Alloc>;
		try {
			for(; count > 0; ++current, --count) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
				_::construct(alloc, std::addressof(*current));
			}
		} catch(...) {
			for(; current != first; ++first) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
				_::destroy(alloc, std::addressof(*first));
			}
			throw;
		}
	}
	return current;
}

template<class ForwardIt, class Size>
auto uninitialized_default_construct_n(ForwardIt first, Size count) -> ForwardIt {
    using T = typename std::iterator_traits<ForwardIt>::value_type;
    ForwardIt current = first;
    try {
        for (; count > 0 ; (void) ++current, --count) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
            ::new (static_cast<void*>(std::addressof(*current))) T;
        }
        return current;
    }  catch (...) {assert(0);
//        std::destroy(first, current);
        throw;
    }
}

}  // end namespace xtd

template<class Alloc> struct alloc_destroy_elem_t {
	Alloc* palloc_;
	template<class T> constexpr auto operator()(T&& ptr) const {  // ->decltype(std::allocator_traits<Alloc>::construct(*palloc_, std::forward<T>(t)...)){
		return std::allocator_traits<Alloc>::destroy(*palloc_, std::addressof(ptr));
	}
};

template<class BidirIt, class Size, class T = typename std::iterator_traits<BidirIt>::value_type>
constexpr auto destroy_n(BidirIt first, Size count)
->std::decay_t<decltype(std::addressof(*(first - 1)), first)> {
	first += count;
	for(; count != 0; --first, --count) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
		std::addressof(*(first-1))->~T();
	}
	return first;
}

template<class Alloc, class BidirIt, class Size, class T = typename std::iterator_traits<BidirIt>::value_type>
constexpr auto alloc_destroy_n(Alloc& alloc, BidirIt first, Size count)
->std::decay_t<decltype(std::addressof(*(first-1)), first)> {
	first += count;
	for (; count != 0; --first, --count) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
		std::allocator_traits<Alloc>::destroy(alloc, std::addressof(*(first - 1)));
	}
	return first;
}

constexpr class adl_uninitialized_copy_t {
	template<class InIt, class FwdIt, class=decltype(std::addressof(*FwdIt{}))>  // sfinae friendy std::uninitialized_copy
	[[nodiscard]]                  constexpr auto _(priority<1>/**/, InIt first, InIt last, FwdIt d_first) const DECLRETURN(       std::uninitialized_copy(first, last, d_first))
	template<class... As>          constexpr auto _(priority<2>/**/,        As&&... args)       const DECLRETURN(                       uninitialized_copy(std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& arg, As&&... args)       const DECLRETURN(  std::decay_t<T>::  uninitialized_copy(std::forward<T>(arg), std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<4>/**/, T&& arg, As&&... args)       const DECLRETURN(std::forward<T>(arg).uninitialized_copy(std::forward<As>(args)...))

 public:
	template<class... As> constexpr auto operator()(As&&... args) const DECLRETURN(_(priority<5>{}, std::forward<As>(args)...))
} adl_uninitialized_copy;

namespace xtd {

template<class InputIt, class Size, class ForwardIt, class Value = typename std::iterator_traits<ForwardIt>::value_type>
auto uninitialized_copy_n(InputIt first, Size count, ForwardIt d_first)
->std::decay_t<decltype(::new (static_cast<void*>(std::addressof(*d_first))) Value(*first), d_first)> {
	ForwardIt current = d_first;
	try {
		for (; count > 0; ++first, (void) ++current, --count) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			::new (static_cast<void*>(std::addressof(*current))) Value(*first);
		}
	} catch(...) {
		for(; d_first != current; ++d_first) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			d_first->~Value();
		}
		throw;
	}
	return current;
}

template<class InputIt, class Size, class ForwardIt, class Value = typename std::iterator_traits<ForwardIt>::value_type>
auto uninitialized_move_n(InputIt first, Size count, ForwardIt d_first)
->std::decay_t<decltype(::new (static_cast<void*>(std::addressof(*d_first))) Value(std::move(*first)), d_first)> {
	ForwardIt current = d_first;
	try {
		for (; count > 0; ++first, (void) ++current, --count) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			::new (static_cast<void*>(std::addressof(*current))) Value(std::move(*first));
		}
	} catch(...) {
		for(; d_first != current; ++d_first) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			d_first->~Value();
		}
		throw;
	}
	return current;
}

}  // end namespace xtd

constexpr class adl_uninitialized_copy_n_t {
	template<class... As>          constexpr auto _(priority<1>/**/,        As&&... args) const DECLRETURN(                  std::uninitialized_copy_n(std::forward<As>(args)...))
#if defined(__NVCC__)
	template<class... As> 		   constexpr auto _(priority<2>/**/,        As&&... args) const DECLRETURN(           thrust::copy_n(                    std::forward<As>(args)...))
#endif
	template<class... As>          constexpr auto _(priority<3>/**/,        As&&... args) const DECLRETURN(                       uninitialized_copy_n(std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<4>/**/, T&& arg, As&&... args) const DECLRETURN(std::decay_t<T>::    uninitialized_copy_n(std::forward<T>(arg), std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<5>/**/, T&& arg, As&&... args) const DECLRETURN(std::forward<T>(arg).uninitialized_copy_n(std::forward<As>(args)...))

 public:
	template<class... As> constexpr auto operator()(As&&... args) const {return _(priority<5>{}, std::forward<As>(args)...);}  // TODO(correaa) this might trigger a compiler crash with g++ 7.5 because of operator&() && overloads
} adl_uninitialized_copy_n;

constexpr class adl_uninitialized_move_n_t {
	template<class... As>          constexpr auto _(priority<1>/**/,          As&&... args) const DECLRETURN(              xtd::  uninitialized_move_n(std::forward<As>(args)...))
	template<class... As>          constexpr auto _(priority<2>/**/,          As&&... args) const DECLRETURN(                     uninitialized_move_n(std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& arg, As&&... args) const DECLRETURN(std::decay_t<T>::    uninitialized_move_n(std::forward<T>(arg), std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<4>/**/, T&& arg, As&&... args) const DECLRETURN(std::forward<T>(arg).uninitialized_move_n(std::forward<As>(args)...))

 public:
	template<class... As> constexpr auto operator()(As&&... args) const {return _(priority<4>{}, std::forward<As>(args)...);}
} adl_uninitialized_move_n;

namespace xtd {

template<class T, class InputIt, class Size, class ForwardIt>
constexpr auto alloc_uninitialized_copy_n(std::allocator<T>& /*alloc*/, InputIt first, Size count, ForwardIt d_first) {
	return adl_uninitialized_copy_n(first, count, d_first);}

template<class T, class InputIt, class Size, class ForwardIt>
constexpr auto alloc_uninitialized_move_n(std::allocator<T>& /*alloc*/, InputIt first, Size count, ForwardIt d_first) {
	return adl_uninitialized_move_n(first, count, d_first);}

template<class Alloc, class InputIt, class Size, class ForwardIt, class = decltype(std::addressof(*ForwardIt{}))>
auto alloc_uninitialized_copy_n(Alloc& alloc, InputIt first, Size count, ForwardIt d_first) {
	ForwardIt current = d_first;
	try {
		for(; count > 0; ++first, ++current, --count) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			std::allocator_traits<Alloc>::construct(alloc, std::addressof(*current), *first);
		}
		return current;
	} catch(...) {
		for(; d_first != current; ++d_first) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			std::allocator_traits<Alloc>::destroy(alloc, std::addressof(*d_first));
		}
		throw;
	}
}

template<class Alloc, class InputIt, class Size, class ForwardIt>
auto alloc_uninitialized_move_n(Alloc& alloc, InputIt first, Size count, ForwardIt d_first) {
	ForwardIt current = d_first;
	try {
		for(; count > 0; ++first, ++current, --count) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			std::allocator_traits<Alloc>::construct(alloc, std::addressof(*current), std::move(*first));
		}
		return current;
	} catch(...) {
		for(; d_first != current; ++d_first) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			std::allocator_traits<Alloc>::destroy(alloc, std::addressof(*d_first));
		}
		throw;
	}
}

template<class T, class InputIt, class ForwardIt>
constexpr auto alloc_uninitialized_copy(std::allocator<T>&/*allocator*/, InputIt first, InputIt last, ForwardIt d_first) {
	return adl_uninitialized_copy(first, last, d_first);
}

template<class Alloc, class InputIt, class ForwardIt, class=decltype(std::addressof(*std::declval<ForwardIt>())), class=std::enable_if_t<std::is_constructible_v<typename std::iterator_traits<ForwardIt>::value_type, typename std::iterator_traits<InputIt>::reference>>>
auto alloc_uninitialized_copy(Alloc& alloc, InputIt first, InputIt last, ForwardIt d_first) {
// ->std::decay_t<decltype(a.construct(std::addressof(*d_first), *first), d_first)> // problematic in clang-11 + gcc-9
	ForwardIt current = d_first;
	try {
		for(; first != last; ++first, (void)++current) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			std::allocator_traits<std::decay_t<Alloc>>::construct(alloc, std::addressof(*current), *first);
		}
		return current;
	} catch(...) {
		for(; d_first != current; ++d_first) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			std::allocator_traits<std::decay_t<Alloc>>::destroy(alloc, std::addressof(*d_first));
		}
		throw;
	}
}

template<class Alloc, class ForwardIt, class Size, class T>
auto alloc_uninitialized_fill_n(Alloc& alloc, ForwardIt first, Size n, T const& value)
->std::decay_t<decltype(std::allocator_traits<Alloc>::construct(alloc, std::addressof(*first), value), first)> {
	ForwardIt current = first;  // using std::to_address;
	try {
		for(; n > 0; ++current, --n) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			std::allocator_traits<Alloc>::construct(alloc, std::addressof(*current), value);
		}
		return current;
	} catch(...) {
		for(; first != current; ++first) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			std::allocator_traits<Alloc>::destroy(alloc, std::addressof(*first));
		}
		throw;
	}
}
}  // end namespace xtd

constexpr class adl_distance_t {
	template<class... As>          constexpr auto _(priority<1>/**/,        As&&... args) const DECLRETURN(              std::    distance(std::forward<As>(args)...))
	template<class... As>          constexpr auto _(priority<2>/**/,        As&&... args) const DECLRETURN(                       distance(std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& arg, As&&... args) const DECLRETURN(  std::decay_t<T>::  distance(std::forward<T>(arg), std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<4>/**/, T&& arg, As&&... args) const DECLRETURN(std::forward<T>(arg).distance(std::forward<As>(args)...))

 public:
	template<class... As> constexpr auto operator()(As&&... args) const DECLRETURN(_(priority<4>{}, std::forward<As>(args)...))
} adl_distance;

constexpr class adl_begin_t {
	template<class... As>          constexpr auto _(priority<1>/**/,          As&&... args) const DECLRETURN(                std::begin(std::forward<As>(args)...))
	template<class... As>          constexpr auto _(priority<2>/**/,          As&&... args) const DECLRETURN(                     begin(std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& arg, As&&... args) const DECLRETURN(    std::decay_t<T>::begin(std::forward<T>(arg), std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<4>/**/, T&& arg, As&&... args) const DECLRETURN(std::forward<T>(arg).begin(std::forward<As>(args)...))

 public:
	template<class... As> [[nodiscard]] constexpr auto operator()(As&&... args) const DECLRETURN(_(priority<4>{}, std::forward<As>(args)...))
} adl_begin;

constexpr class adl_end_t {
	template<class... As>          constexpr auto _(priority<1>/**/,          As&&... args) const DECLRETURN(              std::  end(std::forward<As>(args)...))
	template<class... As>          constexpr auto _(priority<2>/**/,          As&&... args) const DECLRETURN(                     end(std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& arg, As&&... args) const DECLRETURN(  std::decay_t<T>::  end(std::forward<T>(arg), std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<4>/**/, T&& arg, As&&... args) const DECLRETURN(std::forward<T>(arg).end(std::forward<As>(args)...))

 public:
	template<class... As> constexpr auto operator()(As&&... args) const DECLRETURN(_(priority<4>{}, std::forward<As>(args)...))
} adl_end;

constexpr class adl_swap_ranges_t {
	template<class... As>          constexpr auto _(priority<1>/**/,          As&&... args) const DECLRETURN(              std::  swap_ranges(std::forward<As>(args)...))
	template<class... As>          constexpr auto _(priority<2>/**/,          As&&... args) const DECLRETURN(                     swap_ranges(std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& arg, As&&... args) const DECLRETURN(  std::decay_t<T>::  swap_ranges(std::forward<T>(arg), std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<4>/**/, T&& arg, As&&... args) const DECLRETURN(std::forward<T>(arg).swap_ranges(std::forward<As>(args)...))

 public:
	template<class... As> constexpr auto operator()(As&&... args) const DECLRETURN(_(priority<4>{}, std::forward<As>(args)...))
} adl_swap_ranges;

constexpr class adl_lexicographical_compare_t {
	template<class... As>          /*[[gnu::pure]]*/ constexpr auto _(priority<1>/**/,          As&&... args) const DECLRETURN(              std::  lexicographical_compare(std::forward<As>(args)...))
	template<class... As>          constexpr auto _(priority<2>/**/,          As&&... args) const DECLRETURN(                     lexicographical_compare(std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& arg, As&&... args) const DECLRETURN(  std::decay_t<T>::  lexicographical_compare(std::forward<T>(arg), std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<4>/**/, T&& arg, As&&... args) const DECLRETURN(std::forward<T>(arg).lexicographical_compare(std::forward<As>(args)...))

 public:
	template<class... As> /*[[gnu::pure]]*/ constexpr	auto operator()(As&&... args) const DECLRETURN(_(priority<4>{}, std::forward<As>(args)...))
} adl_lexicographical_compare;

constexpr class adl_uninitialized_value_construct_n_t {
	template<class... As>              constexpr auto _(priority<1>/**/,                      As&&... args) const DECLRETURN(              std::  uninitialized_value_construct_n(std::forward<As>(args)...))  // TODO(correaa) use boost alloc_X functions?
	template<class... As>              constexpr auto _(priority<2>/**/,                      As&&... args) const DECLRETURN(                     uninitialized_value_construct_n(std::forward<As>(args)...))
	template<class T, class... As>     constexpr auto _(priority<3>/**/, T&& arg,             As&&... args) const DECLRETURN(  std::decay_t<T>::uninitialized_value_construct_n(std::forward<T>(arg), std::forward<As>(args)...))
	template<class T, class... As>     constexpr auto _(priority<4>/**/, T&& arg,             As&&... args) const DECLRETURN(std::forward<T>(arg).uninitialized_value_construct_n(std::forward<As>(args)...))

 public:
	template<class... As> constexpr auto operator()(As&&... args) const {return (_(priority<4>{}, std::forward<As>(args)...));}
} adl_uninitialized_value_construct_n;

constexpr class adl_alloc_uninitialized_value_construct_n_t {
	template<class Alloc, class... As> constexpr auto _(priority<1>/**/, Alloc&& /*alloc*/, As&&... args) const DECLRETURN(                       adl_uninitialized_value_construct_n(std::forward<As>(args)...))
	template<class... As>              constexpr auto _(priority<2>/**/,                    As&&... args) const DECLRETURN(              xtd::  alloc_uninitialized_value_construct_n(std::forward<As>(args)...))  // TODO(correaa) use boost alloc_X functions?
	template<class... As>              constexpr auto _(priority<3>/**/,                    As&&... args) const DECLRETURN(                     alloc_uninitialized_value_construct_n(std::forward<As>(args)...))
	template<class T, class... As>     constexpr auto _(priority<4>/**/, T&& arg,           As&&... args) const DECLRETURN(  std::decay_t<T>::  alloc_uninitialized_value_construct_n(std::forward<T>(arg), std::forward<As>(args)...))
	template<class T, class... As>     constexpr auto _(priority<5>/**/, T&& arg,           As&&... args) const DECLRETURN(std::forward<T>(arg).alloc_uninitialized_value_construct_n(std::forward<As>(args)...))

 public:
	template<class... As> constexpr auto operator()(As&&... args) const {return (_(priority<5>{}, std::forward<As>(args)...));}
} adl_alloc_uninitialized_value_construct_n;

constexpr class adl_uninitialized_default_construct_n_t {
	template<class... As>          constexpr auto _(priority<1>/**/,          As&&... args) const {return                  xtd::  uninitialized_default_construct_n(                      std::forward<As>(args)...);}
	template<class... As>          constexpr auto _(priority<2>/**/,          As&&... args) const DECLRETURN(                     uninitialized_default_construct_n(                      std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& arg, As&&... args) const DECLRETURN(  std::decay_t<T>::  uninitialized_default_construct_n(std::forward<T>(arg), std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<4>/**/, T&& arg, As&&... args) const DECLRETURN(std::forward<T>(arg).uninitialized_default_construct_n(                      std::forward<As>(args)...))

 public:
		template<class... As> constexpr auto operator()(As&&... args) const {return (_(priority<4>{}, std::forward<As>(args)...));}
} adl_uninitialized_default_construct_n;

constexpr class adl_alloc_uninitialized_default_construct_n_t {
	template<class Alloc, class... As> constexpr auto _(priority<1>/**/, Alloc&&/*unused*/, As&&... args) const JUSTRETURN(         adl_uninitialized_default_construct_n(std::forward<As>(args)...))
//  #if defined(__NVCC__)
//  template<class Alloc, class... As> constexpr auto _(priority<3>/**/,        As&&... as) const DECLRETURN(          thrust::         uninitialized_construct_n_with_allocator(               std::forward<As>(as)...))
//  #endif
	template<class... As>          constexpr auto _(priority<2>/**/,          As&&... args) const DECLRETURN(              xtd::  alloc_uninitialized_default_construct_n(                      std::forward<As>(args)...))  // TODO(correaa) use boost alloc_X functions?
	template<class... As>          constexpr auto _(priority<4>/**/,          As&&... args) const DECLRETURN(                     alloc_uninitialized_default_construct_n(                      std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<5>/**/, T&& arg, As&&... args) const DECLRETURN(  std::decay_t<T>::  alloc_uninitialized_default_construct_n(std::forward<T>(arg), std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<6>/**/, T&& arg, As&&... args) const DECLRETURN(std::forward<T>(arg).alloc_uninitialized_default_construct_n(                      std::forward<As>(args)...))
 public:
	template<class... As> constexpr auto operator()(As&&... args) const {return (_(priority<6>{}, std::forward<As>(args)...));}
} adl_alloc_uninitialized_default_construct_n;

constexpr class destroy_n_t {
	template<class... As>          constexpr auto _(priority<1>/**/,          As&&... args) const DECLRETURN(            multi::  destroy_n              (std::forward<As>(args)...))
	template<class... As>          constexpr auto _(priority<2>/**/,          As&&... args) const DECLRETURN(                     destroy_n              (std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& arg, As&&... args) const DECLRETURN(  std::decay_t<T>::  destroy_n(std::forward<T>(arg), std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<4>/**/, T&& arg, As&&... args) const DECLRETURN(std::forward<T>(arg).destroy_n              (std::forward<As>(args)...))
public:
	template<class... As> constexpr auto operator()(As&&... args) const DECLRETURN(_(priority<4>{}, std::forward<As>(args)...))
} adl_destroy_n;

constexpr class alloc_destroy_n_t {
	template<class Alloc, class... As> constexpr auto _(priority<1>/**/, Alloc&&/*unused*/, As&&... args) const DECLRETURN(             adl_destroy_n              (std::forward<As>(args)...))
	template<             class... As> constexpr auto _(priority<2>/**/,          As&&... args) const DECLRETURN(multi::              alloc_destroy_n              (std::forward<As>(args)...))  // TODO(correaa) use boost alloc_X functions?
	template<             class... As> constexpr auto _(priority<3>/**/,          As&&... args) const DECLRETURN(                     alloc_destroy_n              (std::forward<As>(args)...))
	template<class T,     class... As> constexpr auto _(priority<4>/**/, T&& arg, As&&... args) const DECLRETURN(std::decay_t<T>::    alloc_destroy_n(std::forward<T>(arg), std::forward<As>(args)...))
	template<class T,     class... As> constexpr auto _(priority<5>/**/, T&& arg, As&&... args) const DECLRETURN(std::forward<T>(arg).alloc_destroy_n              (std::forward<As>(args)...))
 public:
	template<class... As> constexpr auto operator()(As&&... args) const DECLRETURN(_(priority<5>{}, std::forward<As>(args)...))
} adl_alloc_destroy_n;

constexpr class adl_alloc_uninitialized_copy_t {
	template<class Alloc, class... As> constexpr auto _(priority<1>/**/, Alloc&&/*alloc*/, As&&... args) const DECLRETURN(                     adl_uninitialized_copy(std::forward<As>(args)...))
// TODO(correaa) : remove T from below?
	template<class T,     class... As> constexpr auto _(priority<2>/**/, T&& arg, As&&... args) const DECLRETURN(                xtd::alloc_uninitialized_copy(std::forward<T>(arg), std::forward<As>(args)...))
	template<class T,     class... As> constexpr auto _(priority<3>/**/, T&& arg, As&&... args) const DECLRETURN(                     alloc_uninitialized_copy(std::forward<T>(arg), std::forward<As>(args)...))
	template<class T,     class... As> constexpr auto _(priority<4>/**/, T&& arg, As&&... args) const DECLRETURN(    std::decay_t<T>::alloc_uninitialized_copy(std::forward<T>(arg), std::forward<As>(args)...))
	template<class T,     class... As> constexpr auto _(priority<5>/**/, T&& arg, As&&... args) const DECLRETURN(std::forward<T>(arg).alloc_uninitialized_copy(std::forward<As>(args)...))
 public:
	template<class... As> constexpr auto operator()(As&&... args) const DECLRETURN(_(priority<5>{}, std::forward<As>(args)...))
} adl_alloc_uninitialized_copy;

constexpr class alloc_uninitialized_copy_n_t {
	template<class Alloc, class... As> constexpr auto _(priority<1>/**/, Alloc&&/*alloc*/, As&&... args) const DECLRETURN(                         uninitialized_copy_n(std::forward<As>(args)...))
	template<class... As>              constexpr auto _(priority<2>/**/,                   As&&... args) const DECLRETURN(              xtd::alloc_uninitialized_copy_n(std::forward<As>(args)...))
#if defined(__NVCC__)
	template<class Alloc, class... As> constexpr auto _(priority<3>/**/, Alloc&&/*alloc*/, As&&... as) const DECLRETURN(                 thrust::uninitialized_copy_n(std::forward<As>(as)...))
#endif
// TODO(correaa) revise
	template<class Alloc, class... As> constexpr auto _(priority<4>/**/, Alloc&&/*alloc*/, As&&... args) const DECLRETURN(                           uninitialized_copy_n(std::forward<As>(args)...))
	template<class... As>              constexpr auto _(priority<5>/**/,                   As&&... args) const DECLRETURN(                     alloc_uninitialized_copy_n(std::forward<As>(args)...))
	template<class T, class... As>     constexpr auto _(priority<6>/**/, T&& arg,          As&&... args) const DECLRETURN(std::forward<T>(arg).alloc_uninitialized_copy_n(std::forward<As>(args)...))

 public:
	template<class... As> constexpr auto operator()(As&&... args) const {return _(priority<6>{}, std::forward<As>(args)...);}
} adl_alloc_uninitialized_copy_n;

constexpr class alloc_uninitialized_move_n_t {
// TODO(correaa) : fallback to no alloc version
	template<class... As>          constexpr auto _(priority<1>/**/,          As&&... args) const {return(                 xtd::  alloc_uninitialized_move_n(std::forward<As>(args)...));}
	template<class... As>          constexpr auto _(priority<2>/**/,          As&&... args) const DECLRETURN(                     alloc_uninitialized_move_n(std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& arg, As&&... args) const DECLRETURN(std::forward<T>(arg).alloc_uninitialized_move_n(std::forward<As>(args)...))

 public:
	template<class... As> constexpr auto operator()(As&&... args) const {return _(priority<3>{}, std::forward<As>(args)...);} \
} adl_alloc_uninitialized_move_n;

constexpr class uninitialized_fill_n_t {
	template<class... As>          constexpr auto _(priority<1>/**/,          As&&... args) const DECLRETURN(               std::    uninitialized_fill_n(std::forward<As>(args)...))
	template<class... As>          constexpr auto _(priority<2>/**/,          As&&... args) const DECLRETURN(                        uninitialized_fill_n(std::forward<As>(args)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& arg, As&&... args) const DECLRETURN( std::forward<T>(arg).uninitialized_fill_n(std::forward<As>(args)...))

 public:
	template<class T1, class... As> constexpr auto operator()(T1&& arg, As&&... args) const DECLRETURN(_(priority<3>{}, arg, std::forward<As>(args)...))
} adl_uninitialized_fill_n;

constexpr class alloc_uninitialized_fill_n_t {
	template<             class... As> constexpr auto _(priority<1>/**/,                   As&&... args) const DECLRETURN(                       xtd::alloc_uninitialized_fill_n(std::forward<As>(args)...))
	template<class Alloc, class... As> constexpr auto _(priority<2>/**/, Alloc&&/*alloc*/, As&&... args) const DECLRETURN(                              adl_uninitialized_fill_n(std::forward<As>(args)...))
	template<             class... As> constexpr auto _(priority<3>/**/,                   As&&... args) const DECLRETURN(                            alloc_uninitialized_fill_n(std::forward<As>(args)...))
	template<class Alloc, class... As> constexpr auto _(priority<4>/**/, Alloc&&  alloc  , As&&... args) const DECLRETURN( std::forward<Alloc>(alloc).alloc_uninitialized_fill_n(std::forward<As>(args)...))
 public:
	template<class T1, class... As> constexpr auto operator()(T1&& arg, As&&... args) const DECLRETURN(_(priority<4>{}, arg, std::forward<As>(args)...))
} adl_alloc_uninitialized_fill_n;

template<dimensionality_type N>
struct recursive {
	template<class Alloc, class InputIt, class ForwardIt>
	static constexpr auto alloc_uninitialized_copy(Alloc& alloc, InputIt first, InputIt last, ForwardIt dest){
		using std::begin; using std::end;
		while(first!=last) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			recursive<N-1>::alloc_uninitialized_copy(alloc, begin(*first), end(*first), begin(*dest));
			++first;
			++dest;
		}
		return dest;
	}
};

template<> struct recursive<1> {
	template<class Alloc, class InputIt, class ForwardIt>
	static auto alloc_uninitialized_copy(Alloc& alloc, InputIt first, InputIt last, ForwardIt dest){
		return adl_alloc_uninitialized_copy(alloc, first, last, dest);
	}
};

}  // end namespace boost::multi
#endif
