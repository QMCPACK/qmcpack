// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2020-2021 Alfredo A. Correa

#ifndef MULTI_DETAIL_ADL_HPP
#define MULTI_DETAIL_ADL_HPP

#include<cstddef>      // std::size_t
#include<type_traits>  // std::conditional_t
#include<utility>

#include "../config/MAYBE_UNUSED.hpp"
#include "../config/NODISCARD.hpp"
#include "../detail/memory.hpp"

#if defined(__NVCC__)
#include<thrust/copy.h>
#include<thrust/equal.h>
#include<thrust/detail/memory_algorithms.h>
#include<thrust/uninitialized_copy.h>
#endif

#include<algorithm>  // for std::copy, std::copy_n, std::equal, etc
#include<iterator>   // for begin, end
#include<memory>     // for uninitialized_copy, etc

#define BOOST_MULTI_DEFINE_ADL(FuN) /*NOLINT(cppcoreguidelines-macro-usage) TODO(correaa) consider replacing for all ADL'd operations*/ \
namespace boost { \
namespace multi { \
namespace adl { \
	namespace custom {template<class...> struct FuN##_t;} 	__attribute__((unused))  \
	static constexpr class FuN##_t { \
		template<class... As> [[deprecated]] auto _(priority<0>,        As&&... as) const = delete; \
		template<class... As>          auto _(priority<1>,        As&&... as) const DECLRETURN(std::FuN(std::forward<As>(as)...)) \
		template<class... As>          auto _(priority<2>,        As&&... as) const DECLRETURN(     FuN(std::forward<As>(as)...)) \
		template<class T, class... As> auto _(priority<3>, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).FuN(std::forward<As>(as)...))     \
		template<class... As>          auto _(priority<4>,        As&&... as) const DECLRETURN(custom::FuN##_t<As&&...>::_(std::forward<As>(as)...)) \
	public: \
		template<class... As> auto operator()(As&&... as) const->decltype(_(priority<4>{}, std::forward<As>(as)...)) {return _(priority<4>{}, std::forward<As>(as)...);} \
	} (FuN); \
}  /* end namespace adl   */ \
}  /* end namespace multi */ \
}  /* end namespace boost */

namespace boost::multi {

template<std::size_t I> struct priority : std::conditional_t<I==0, std::true_type, priority<I-1>> {};

#define DECLRETURN(ExpR) ->decltype(ExpR) {return ExpR;}  // NOLINT(cppcoreguidelines-macro-usage) saves a lot of typing
#define JUSTRETURN(ExpR)                  {return ExpR;}  // NOLINT(cppcoreguidelines-macro-usage) saves a lot of typing

constexpr class adl_copy_n_t {
	template<class... As>          constexpr auto _(priority<0>/**/,        As&&... as) const DECLRETURN(              std::copy_n(                    std::forward<As>(as)...))
#if defined(__NVCC__)
	template<class... As> 		   constexpr auto _(priority<1>/**/,        As&&... as) const DECLRETURN(         ::thrust::copy_n(                    std::forward<As>(as)...))
#endif
	template<class... As>          constexpr auto _(priority<2>/**/,        As&&... as) const DECLRETURN(                   copy_n(                    std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& t, As&&... as) const DECLRETURN(std::decay_t<T>::  copy_n(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>/**/, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).copy_n(                    std::forward<As>(as)...))

 public:
	template<class... As> constexpr auto operator()(As&&... as) const DECLRETURN(_(priority<4>{}, std::forward<As>(as)...))
} adl_copy_n;

constexpr class adl_move_t {
	template<class... As>           constexpr auto _(priority<0>/**/,                    As&&... as) const DECLRETURN(              std::move(                    std::forward<As>(as)...))
#if defined(__NVCC__)  // there is no thrust::move algorithm
	template<class It, class... As> constexpr auto _(priority<1>/**/, It first, It last, As&&... as) const DECLRETURN(           thrust::copy(std::make_move_iterator(first), std::make_move_iterator(last), std::forward<As>(as)...))
#endif
	template<class... As>           constexpr auto _(priority<2>/**/,                    As&&... as) const DECLRETURN(                   move(                    std::forward<As>(as)...))
	template<class T, class... As>  constexpr auto _(priority<3>/**/, T&& t,             As&&... as) const DECLRETURN(std::decay_t<T>::  move(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As>  constexpr auto _(priority<4>/**/, T&& t,             As&&... as) const DECLRETURN(std::forward<T>(t).move(                    std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const DECLRETURN(_(priority<4>{}, std::forward<As>(as)...))
} adl_move;

constexpr class adl_fill_n_t {
	template<         class... As> constexpr auto _(priority<0>/**/,        As&&... as) const DECLRETURN(              std::fill_n              (std::forward<As>(as)...))
#if defined(__NVCC__)
	template<         class... As> constexpr auto _(priority<1>/**/,        As&&... as) const DECLRETURN(           thrust::fill_n              (std::forward<As>(as)...))
#endif
	template<         class... As> constexpr auto _(priority<2>/**/,        As&&... as) const DECLRETURN(                   fill_n              (std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& t, As&&... as) const DECLRETURN(std::decay_t<T>::  fill_n(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>/**/, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).fill_n              (std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const DECLRETURN(_(priority<4>{}, std::forward<As>(as)...))
} adl_fill_n;

constexpr class adl_equal_t {
	template<         class...As> constexpr auto _(priority<1>/**/,        As&&...as) const DECLRETURN(               std::equal(                    std::forward<As>(as)...))
#if defined(__NVCC__)
	template<         class...As> constexpr auto _(priority<2>/**/,        As&&...as) const DECLRETURN(          ::thrust::equal(                    std::forward<As>(as)...))
#endif
	template<         class...As> constexpr auto _(priority<3>/**/,        As&&...as) const DECLRETURN(                    equal(                    std::forward<As>(as)...))
	template<class T, class...As> constexpr auto _(priority<4>/**/, T&& t, As&&...as) const DECLRETURN( std::decay_t<T>::  equal(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class...As> constexpr auto _(priority<5>/**/, T&& t, As&&...as) const DECLRETURN( std::forward<T>(t).equal(                    std::forward<As>(as)...))
public:
	template<class...As> constexpr auto operator()(As&&...as) const DECLRETURN(_(priority<5>{}, std::forward<As>(as)...))
} adl_equal;

template<class... Args> struct adl_custom_copy;

constexpr class adl_copy_t {
	class Copy;
	template<class InputIt, class OutputIt,
		class=std::enable_if_t<std::is_assignable<typename std::iterator_traits<OutputIt>::reference, typename std::iterator_traits<InputIt>::reference>{}>
	>
	NODISCARD("")                  constexpr auto _(priority<1>/**/, InputIt first, InputIt last, OutputIt d_first)	                                                                                const DECLRETURN(              std::copy(first, last, d_first))
#if defined(__NVCC__)
	template<class... As> 		   constexpr auto _(priority<2>/**/,        As&&... as) const DECLRETURN(           thrust::copy(std::forward<As>(as)...))
#endif
	template<         class... As> constexpr auto _(priority<3>/**/,        As&&... as) const DECLRETURN(                   copy(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>/**/, T&& t, As&&... as) const DECLRETURN(  std::decay_t<T>::copy(std::forward<T>(t), std::forward<As>(as)...))
//	template<class... As         > constexpr auto _(priority<5>/**/,        As&&... as) const DECLRETURN(boost::multi::adl_custom_copy<std::decay_t<As>...>::copy(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<6>/**/, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).copy(std::forward<As>(as)...))

 public:
	template<class... As> constexpr auto operator()(As&&... as) const DECLRETURN( _(priority<6>{}, std::forward<As>(as)...) ) \
} adl_copy;

namespace adl {
	namespace custom {template<class...> struct fill_t;}  // __attribute__((unused))
	static constexpr class fill_t {
		template<class... As>          auto _(priority<1>/**/,        As&&... as) const DECLRETURN(              std::fill              (std::forward<As>(as)...))
		template<class... As>          auto _(priority<2>/**/,        As&&... as) const DECLRETURN(                   fill              (std::forward<As>(as)...))
		template<class T, class... As> auto _(priority<3>/**/, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).fill              (std::forward<As>(as)...))
		template<class... As>          auto _(priority<4>/**/,        As&&... as) const DECLRETURN(custom::           fill_t<As&&...>::_(std::forward<As>(as)...))
	public:
		template<class... As> auto operator()(As&&... as) const DECLRETURN(_(priority<5>{}, std::forward<As>(as)...))
	} fill [[maybe_unused]];
}  // end namespace adl

template<class Alloc>
struct alloc_construct_elem_t {
	Alloc* palloc_;
	template<class T> auto operator()(T&& p) const
	->decltype(std::allocator_traits<Alloc>::construct(*palloc_, std::addressof(p))) {
		return std::allocator_traits<Alloc>::construct(*palloc_, std::addressof(p)); }
};

namespace xtd {

template<class T>  // this one goes last!!!
constexpr auto to_address(const T& p) noexcept;

template<class T>
constexpr auto me_to_address(priority<0>/**/, const T& p) noexcept
->decltype(to_address(p.operator->())) {
	return to_address(p.operator->()); }

template<class T>
constexpr auto me_to_address(priority<1>/**/, const T& p) noexcept
->decltype(std::pointer_traits<T>::to_address(p)) {
	return std::pointer_traits<T>::to_address(p); }

template<class T, std::enable_if_t<std::is_pointer<T>{}, int> = 0>
constexpr auto me_to_address(priority<2>/**/, T const& p) noexcept -> T{
    static_assert(!std::is_function<T>{}, "!");
    return p;
}

template<class T>  // this one goes last!!!
constexpr auto to_address(const T& p) noexcept
->decltype(me_to_address(priority<2>{}/**/, p)) {
	return me_to_address(priority<2>{}, p);
}

template<class Alloc, class ForwardIt, class Size, typename Value = typename std::iterator_traits<ForwardIt>::value_type, typename = decltype(std::addressof(*ForwardIt{}))>
auto alloc_uninitialized_value_construct_n(Alloc& alloc, ForwardIt first, Size n) -> ForwardIt {
// ->std::decay_t<decltype(std::allocator_traits<Alloc>::construct(alloc, std::addressof(*first), Value()), first)>
	ForwardIt current = first;
	try {
		for (; n > 0 ; ++current, --n) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
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
auto alloc_uninitialized_default_construct_n(Alloc& alloc, ForwardIt first, Size n)
->std::decay_t<decltype(std::allocator_traits<Alloc>::construct(alloc, std::addressof(*first)), first)> {
	ForwardIt current = first;
	if(std::is_trivially_default_constructible<T>{}) {
		std::advance(current, n);
	} else {
		using _ = std::allocator_traits<Alloc>;
		try {
			for(; n > 0; ++current, --n) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
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
auto uninitialized_default_construct_n( ForwardIt first, Size n ) -> ForwardIt {
    using T = typename std::iterator_traits<ForwardIt>::value_type;
    ForwardIt current = first;
    try {
        for (; n > 0 ; (void) ++current, --n) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
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
	template<class T> constexpr auto operator()(T&& p) const {  // ->decltype(std::allocator_traits<Alloc>::construct(*palloc_, std::forward<T>(t)...)){
		return std::allocator_traits<Alloc>::destroy(*palloc_, std::addressof(p));
	}
};

template<class BidirIt, class Size, class T = typename std::iterator_traits<BidirIt>::value_type>
constexpr auto destroy_n(BidirIt first, Size n)
->std::decay_t<decltype(std::addressof(*(first-1)), first)> {
	first += n;
	for(; n != 0; --first, --n) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
		std::addressof(*(first-1))->~T();
	}
	return first;
}

template<class Alloc, class BidirIt, class Size, class T = typename std::iterator_traits<BidirIt>::value_type>
constexpr auto alloc_destroy_n(Alloc& a, BidirIt first, Size n)
->std::decay_t<decltype(std::addressof(*(first-1)), first)> {
	first += n;
	for (; n != 0; --first, --n) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
		std::allocator_traits<Alloc>::destroy(a, std::addressof(*(first-1)));
	}
	return first;
}

constexpr class adl_uninitialized_copy_t {
	template<class InIt, class FwdIt, class=decltype(std::addressof(*FwdIt{}))>  // sfinae friendy std::uninitialized_copy
	NODISCARD("")                  constexpr auto _(priority<1>/**/, InIt f, InIt l, FwdIt d) const DECLRETURN(              std::uninitialized_copy(f, l, d))
	template<class... As>          constexpr auto _(priority<2>/**/,        As&&... as)       const DECLRETURN(                   uninitialized_copy(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& t, As&&... as)       const DECLRETURN(  std::decay_t<T>::uninitialized_copy(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>/**/, T&& t, As&&... as)       const DECLRETURN(std::forward<T>(t).uninitialized_copy(std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const DECLRETURN(_(priority<5>{}, std::forward<As>(as)...))
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
	template<class... As>          constexpr auto _(priority<1>/**/,        As&&... as) const DECLRETURN(              std::uninitialized_copy_n(std::forward<As>(as)...))
#if defined(__NVCC__)
	template<class... As> 		   constexpr auto _(priority<2>/**/,        As&&... as) const DECLRETURN(           thrust::copy_n(                    std::forward<As>(as)...))
#endif
	template<class... As>          constexpr auto _(priority<3>/**/,        As&&... as) const DECLRETURN(                   uninitialized_copy_n(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>/**/, T&& t, As&&... as) const DECLRETURN(std::decay_t<T>::  uninitialized_copy_n(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<5>/**/, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).uninitialized_copy_n(std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const {return _(priority<5>{}, std::forward<As>(as)...);}  // TODO(correaa) this might trigger a compiler crash with g++ 7.5 because of operator&() && overloads
} adl_uninitialized_copy_n;

constexpr class adl_uninitialized_move_n_t {
	template<class... As>          constexpr auto _(priority<1>/**/,        As&&... as) const DECLRETURN(              xtd::uninitialized_move_n(std::forward<As>(as)...))
	template<class... As>          constexpr auto _(priority<2>/**/,        As&&... as) const DECLRETURN(                   uninitialized_move_n(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& t, As&&... as) const DECLRETURN(std::decay_t<T>::  uninitialized_move_n(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>/**/, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).uninitialized_move_n(std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const {return _(priority<4>{}, std::forward<As>(as)...);}
} adl_uninitialized_move_n;

namespace xtd {

template<class T, class InputIt, class Size, class ForwardIt>
constexpr auto alloc_uninitialized_copy_n(std::allocator<T>&/*alloc*/, InputIt f, Size n, ForwardIt d) {
	return adl_uninitialized_copy_n(f, n, d);}

template<class T, class InputIt, class Size, class ForwardIt>
constexpr auto alloc_uninitialized_move_n(std::allocator<T>&/*alloc*/, InputIt f, Size n, ForwardIt d) {
	return adl_uninitialized_move_n(f, n, d);}

template<class Alloc, class InputIt, class Size, class ForwardIt, class = decltype(std::addressof(*ForwardIt{}))>
auto alloc_uninitialized_copy_n(Alloc& a, InputIt f, Size n, ForwardIt d) {
// ->std::decay_t<decltype(a.construct(std::addressof(*d), *f), d)>
	ForwardIt c = d;
	try {
		for(; n > 0; ++f, ++c, --n) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			std::allocator_traits<Alloc>::construct(a, std::addressof(*c), *f);
		}
		return c;
	} catch(...) {
		for(; d != c; ++d) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			std::allocator_traits<Alloc>::destroy(a, std::addressof(*d));
		}
		throw;
	}
}

template<class Alloc, class InputIt, class Size, class ForwardIt>
auto alloc_uninitialized_move_n(Alloc& a, InputIt f, Size n, ForwardIt d) {
// ->std::decay_t<decltype(a.construct(std::addressof(*d), *f), d)>
	ForwardIt c = d;
	try {
		for(; n > 0; ++f, ++c, --n) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			std::allocator_traits<Alloc>::construct(a, std::addressof(*c), std::move(*f));
		}
		return c;
	} catch(...) {
		for(; d != c; ++d) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			std::allocator_traits<Alloc>::destroy(a, std::addressof(*d));
		}
		throw;
	}
}

template<class T, class InputIt, class ForwardIt>
constexpr auto alloc_uninitialized_copy(std::allocator<T>&/*allocator*/, InputIt f, InputIt l, ForwardIt d) {
	return adl_uninitialized_copy(f, l, d);
}

template<class Alloc, class InputIt, class ForwardIt, class=decltype(std::addressof(*std::declval<ForwardIt>())), class=std::enable_if_t<std::is_constructible<typename std::iterator_traits<ForwardIt>::value_type, typename std::iterator_traits<InputIt>::reference>{}> >
auto alloc_uninitialized_copy(Alloc& a, InputIt first, InputIt last, ForwardIt d_first) {
// ->std::decay_t<decltype(a.construct(std::addressof(*d_first), *first), d_first)> // problematic in clang-11 + gcc-9
	ForwardIt current = d_first;
	try {
		for(; first != last; ++first, (void)++current) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			std::allocator_traits<std::decay_t<Alloc>>::construct(a, std::addressof(*current), *first);
		}
		return current;
	} catch(...) {
		for(; d_first != current; ++d_first) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			std::allocator_traits<std::decay_t<Alloc>>::destroy(a, std::addressof(*d_first));
		}
		throw;
	}
}

template<class Alloc, class ForwardIt, class Size, class T>
auto alloc_uninitialized_fill_n(Alloc& a, ForwardIt first, Size n, T const& v)
->std::decay_t<decltype(std::allocator_traits<Alloc>::construct(a, std::addressof(*first), v), first)> {
	ForwardIt current = first;  // using std::to_address;
	try {
		for(; n > 0; ++current, --n) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			std::allocator_traits<Alloc>::construct(a, std::addressof(*current), v);
		}
		return current;
	} catch(...) {
		for(; first != current; ++first) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			std::allocator_traits<Alloc>::destroy(a, std::addressof(*first));
		}
		throw;
	}
}
}  // end namespace xtd

constexpr class adl_distance_t {
	template<class... As>          constexpr auto _(priority<1>/**/,        As&&... as) const DECLRETURN(              std::distance(std::forward<As>(as)...))
	template<class... As>          constexpr auto _(priority<2>/**/,        As&&... as) const DECLRETURN(                   distance(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& t, As&&... as) const DECLRETURN(  std::decay_t<T>::distance(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>/**/, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).distance(std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const DECLRETURN(_(priority<4>{}, std::forward<As>(as)...))
} adl_distance;

constexpr class adl_begin_t {
	template<class... As>          NODISCARD("") constexpr auto _(priority<1>/**/,        As&&... as) const DECLRETURN(              std::begin(std::forward<As>(as)...))
	template<class... As>          NODISCARD("") constexpr auto _(priority<2>/**/,        As&&... as) const DECLRETURN(                   begin(std::forward<As>(as)...))
	template<class T, class... As>               constexpr auto _(priority<3>/**/, T&& t, As&&... as) const DECLRETURN(  std::decay_t<T>::begin(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> NODISCARD("") constexpr auto _(priority<4>/**/, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).begin(std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const DECLRETURN(_(priority<4>{}, std::forward<As>(as)...))
} adl_begin;

constexpr class adl_end_t {
	template<class... As>                        constexpr auto _(priority<1>/**/,        As&&... as) const DECLRETURN(              std::end(std::forward<As>(as)...))
	template<class... As>          NODISCARD("") constexpr auto _(priority<2>/**/,        As&&... as) const DECLRETURN(                   end(std::forward<As>(as)...))
	template<class T, class... As>               constexpr auto _(priority<3>/**/, T&& t, As&&... as) const DECLRETURN(  std::decay_t<T>::end(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> NODISCARD("") constexpr auto _(priority<4>/**/, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).end(std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const DECLRETURN(_(priority<4>{}, std::forward<As>(as)...))
} adl_end;

constexpr class adl_swap_ranges_t {
	template<class... As>          constexpr auto _(priority<1>/**/,        As&&... as) const DECLRETURN(              std::swap_ranges(std::forward<As>(as)...))
	template<class... As>          constexpr auto _(priority<2>/**/,        As&&... as) const DECLRETURN(                   swap_ranges(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& t, As&&... as) const DECLRETURN(  std::decay_t<T>::swap_ranges(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>/**/, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).swap_ranges(std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const DECLRETURN(_(priority<4>{}, std::forward<As>(as)...))
} adl_swap_ranges;

constexpr class adl_lexicographical_compare_t {
	template<class... As>          constexpr auto _(priority<1>/**/,        As&&... as) const DECLRETURN(              std::lexicographical_compare(std::forward<As>(as)...))
	template<class... As>          constexpr auto _(priority<2>/**/,        As&&... as) const DECLRETURN(                   lexicographical_compare(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& t, As&&... as) const DECLRETURN(  std::decay_t<T>::lexicographical_compare(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>/**/, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).lexicographical_compare(std::forward<As>(as)...))
public:
	template<class... As> constexpr	auto operator()(As&&... as) const DECLRETURN(_(priority<4>{}, std::forward<As>(as)...))
} adl_lexicographical_compare;

constexpr class adl_alloc_uninitialized_value_construct_n_t {
	template<class... As>          constexpr auto _(priority<1>/**/,        As&&... as) const DECLRETURN(              xtd::alloc_uninitialized_value_construct_n(std::forward<As>(as)...))  // TODO(correaa) use boost alloc_X functions?
	template<class... As>          constexpr auto _(priority<2>/**/,        As&&... as) const DECLRETURN(                   alloc_uninitialized_value_construct_n(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& t, As&&... as) const DECLRETURN(  std::decay_t<T>::alloc_uninitialized_value_construct_n(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>/**/, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).alloc_uninitialized_value_construct_n(std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const {return (_(priority<4>{}, std::forward<As>(as)...));}
} adl_alloc_uninitialized_value_construct_n;

constexpr class adl_uninitialized_default_construct_n_t {
	template<class... As>          constexpr auto _(priority<1>/**/,        As&&... as) const {return                  xtd::uninitialized_default_construct_n              (std::forward<As>(as)...);}
	template<class... As>          constexpr auto _(priority<2>/**/,        As&&... as) const DECLRETURN(                   uninitialized_default_construct_n(                    std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& t, As&&... as) const DECLRETURN(  std::decay_t<T>::uninitialized_default_construct_n(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>/**/, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).uninitialized_default_construct_n              (std::forward<As>(as)...))

 public:
		template<class... As> constexpr auto operator()(As&&... as) const {return (_(priority<4>{}, std::forward<As>(as)...));}
} adl_uninitialized_default_construct_n;

constexpr class adl_alloc_uninitialized_default_construct_n_t {
	template<class Alloc, class... As> constexpr auto _(priority<1>/**/, Alloc&&/*unused*/, As&&... as) const JUSTRETURN(adl_uninitialized_default_construct_n(std::forward<As>(as)...))
//  #if defined(__NVCC__)
//  template<class Alloc, class... As> constexpr auto _(priority<3>/**/,        As&&... as) const DECLRETURN(          thrust::uninitialized_construct_n_with_allocator(                    std::forward<As>(as)...))
//  #endif
	template<class... As>          constexpr auto _(priority<2>/**/,        As&&... as) const DECLRETURN(              xtd::alloc_uninitialized_default_construct_n(                    std::forward<As>(as)...))  // TODO(correaa) use boost alloc_X functions?
	template<class... As>          constexpr auto _(priority<4>/**/,        As&&... as) const DECLRETURN(                   alloc_uninitialized_default_construct_n(                    std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<5>/**/, T&& t, As&&... as) const DECLRETURN(  std::decay_t<T>::alloc_uninitialized_default_construct_n(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<6>/**/, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).alloc_uninitialized_default_construct_n(                    std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const {return (_(priority<6>{}, std::forward<As>(as)...));}
} adl_alloc_uninitialized_default_construct_n;

constexpr class destroy_n_t {
	template<class... As>          constexpr auto _(priority<1>/**/,        As&&... as) const DECLRETURN(            multi::destroy_n              (std::forward<As>(as)...))
	template<class... As>          constexpr auto _(priority<2>/**/,        As&&... as) const DECLRETURN(                   destroy_n              (std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& t, As&&... as) const DECLRETURN(  std::decay_t<T>::destroy_n(std::forward<T>(t), std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<4>/**/, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).destroy_n              (std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const DECLRETURN(_(priority<4>{}, std::forward<As>(as)...))
} adl_destroy_n;

constexpr class alloc_destroy_n_t {
	template<class Alloc, class... As> constexpr auto _(priority<1>/**/, Alloc&&/*unused*/, As&&... as) const DECLRETURN(                     adl_destroy_n              (std::forward<As>(as)...))
	template<             class... As> constexpr auto _(priority<2>/**/,        As&&... as) const DECLRETURN(multi::            alloc_destroy_n              (std::forward<As>(as)...))  // TODO(correaa) use boost alloc_X functions?
	template<             class... As> constexpr auto _(priority<3>/**/,        As&&... as) const DECLRETURN(                   alloc_destroy_n              (std::forward<As>(as)...))
	template<class T,     class... As> constexpr auto _(priority<4>/**/, T&& t, As&&... as) const DECLRETURN(std::decay_t<T>::  alloc_destroy_n(std::forward<T>(t), std::forward<As>(as)...))
	template<class T,     class... As> constexpr auto _(priority<5>/**/, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).alloc_destroy_n              (std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const DECLRETURN(_(priority<5>{}, std::forward<As>(as)...))
} adl_alloc_destroy_n;

constexpr class adl_alloc_uninitialized_copy_t {
	template<class Alloc, class... As> constexpr auto _(priority<1>/**/, Alloc&&/*alloc*/, As&&... as) const DECLRETURN(                     adl_uninitialized_copy(std::forward<As>(as)...))
// TODO(correaa) : remove T from below?
	template<class T,     class... As> constexpr auto _(priority<2>/**/, T&& t, As&&... as) const DECLRETURN(              xtd::alloc_uninitialized_copy(std::forward<T>(t), std::forward<As>(as)...))
	template<class T,     class... As> constexpr auto _(priority<3>/**/, T&& t, As&&... as) const DECLRETURN(                   alloc_uninitialized_copy(std::forward<T>(t), std::forward<As>(as)...))
	template<class T,     class... As> constexpr auto _(priority<4>/**/, T&& t, As&&... as) const DECLRETURN(  std::decay_t<T>::alloc_uninitialized_copy(std::forward<T>(t), std::forward<As>(as)...))
	template<class T,     class... As> constexpr auto _(priority<5>/**/, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).alloc_uninitialized_copy(std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const DECLRETURN(_(priority<5>{}, std::forward<As>(as)...))
} adl_alloc_uninitialized_copy;

constexpr class alloc_uninitialized_copy_n_t {
	template<class Alloc, class... As> constexpr auto _(priority<1>/**/, Alloc&&/*alloc*/, As&&... as) const DECLRETURN(                         uninitialized_copy_n(std::forward<As>(as)...))
	template<class... As>              constexpr auto _(priority<2>/**/,                   As&&... as) const DECLRETURN(              xtd::alloc_uninitialized_copy_n(std::forward<As>(as)...))
#if defined(__NVCC__)
	template<class Alloc, class... As> constexpr auto _(priority<3>/**/, Alloc&&/*alloc*/, As&&... as) const DECLRETURN(                 thrust::uninitialized_copy_n(std::forward<As>(as)...))
#endif
// TODO(correaa) revise
	template<class Alloc, class... As> constexpr auto _(priority<4>/**/, Alloc&&/*alloc*/, As&&... as) const DECLRETURN(                         uninitialized_copy_n(std::forward<As>(as)...))
	template<class... As>              constexpr auto _(priority<5>/**/,                   As&&... as) const DECLRETURN(                   alloc_uninitialized_copy_n(std::forward<As>(as)...))
	template<class T, class... As>     constexpr auto _(priority<6>/**/, T&& t,            As&&... as) const DECLRETURN(std::forward<T>(t).alloc_uninitialized_copy_n(std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const {return _(priority<6>{}, std::forward<As>(as)...);}
} adl_alloc_uninitialized_copy_n;

constexpr class alloc_uninitialized_move_n_t {
// TODO(correaa) : fallback to no alloc version
	template<class... As>          constexpr auto _(priority<1>/**/,        As&&... as) const {return(                  xtd::alloc_uninitialized_move_n(std::forward<As>(as)...));}
	template<class... As>          constexpr auto _(priority<2>/**/,        As&&... as) const DECLRETURN(                   alloc_uninitialized_move_n(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& t, As&&... as) const DECLRETURN(std::forward<T>(t).alloc_uninitialized_move_n(std::forward<As>(as)...))
public:
	template<class... As> constexpr auto operator()(As&&... as) const {return _(priority<3>{}, std::forward<As>(as)...);} \
} adl_alloc_uninitialized_move_n;

constexpr class uninitialized_fill_n_t {
	template<class... As>          constexpr auto _(priority<1>/**/,        As&&... as) const DECLRETURN(               std::uninitialized_fill_n(std::forward<As>(as)...))
	template<class... As>          constexpr auto _(priority<2>/**/,        As&&... as) const DECLRETURN(                    uninitialized_fill_n(std::forward<As>(as)...))
	template<class T, class... As> constexpr auto _(priority<3>/**/, T&& t, As&&... as) const DECLRETURN( std::forward<T>(t).uninitialized_fill_n(std::forward<As>(as)...))
public:
	template<class T1, class... As> constexpr auto operator()(T1&& t1, As&&... as) const DECLRETURN(_(priority<3>{}, t1, std::forward<As>(as)...))
} adl_uninitialized_fill_n;

constexpr class alloc_uninitialized_fill_n_t {
	template<class Alloc, class... As> constexpr auto _(priority<1>/**/, Alloc&&/*alloc*/, As&&... as) const DECLRETURN(                              adl_uninitialized_fill_n(std::forward<As>(as)...))
	template<             class... As> constexpr auto _(priority<2>/**/,                   As&&... as) const DECLRETURN(                       xtd::alloc_uninitialized_fill_n(std::forward<As>(as)...))
	template<             class... As> constexpr auto _(priority<3>/**/,                   As&&... as) const DECLRETURN(                            alloc_uninitialized_fill_n(std::forward<As>(as)...))
	template<class Alloc, class... As> constexpr auto _(priority<4>/**/, Alloc&&  alloc  , As&&... as) const DECLRETURN( std::forward<Alloc>(alloc).alloc_uninitialized_fill_n(std::forward<As>(as)...))
public:
	template<class T1, class... As> constexpr auto operator()(T1&& t1, As&&... as) const DECLRETURN(_(priority<4>{}, t1, std::forward<As>(as)...))
} adl_alloc_uninitialized_fill_n;

template<dimensionality_type N>
struct recursive {
	template<class Alloc, class InputIt, class ForwardIt>
	static constexpr auto alloc_uninitialized_copy(Alloc& a, InputIt first, InputIt last, ForwardIt dest){
		using std::begin; using std::end;
		while(first!=last) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
			recursive<N-1>::alloc_uninitialized_copy(a, begin(*first), end(*first), begin(*dest));
			++first;
			++dest;
		}
		return dest;
	}
};

template<> struct recursive<1> {
	template<class Alloc, class InputIt, class ForwardIt>
	static auto alloc_uninitialized_copy(Alloc& a, InputIt first, InputIt last, ForwardIt dest){
		return adl_alloc_uninitialized_copy(a, first, last, dest);
	}
};

}  // end namespace boost::multi
#endif
