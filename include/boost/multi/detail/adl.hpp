// Copyright 2020-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_DETAIL_ADL_HPP
#define BOOST_MULTI_DETAIL_ADL_HPP

// super fancy block to detect Thrust
#if defined(__CUDACC__) ||              /* nvcc, or clang --cuda (device compilation) */                                        \
	defined(__NVCC__) ||                /* nvcc legacy macro                            */                                      \
	defined(__CUDA__) ||                /* clang --cuda (host-side TU)                 */                                       \
	defined(__HIP_PLATFORM_AMD__) ||    /* hipcc targeting AMD, or clang+HIP on AMD    */                                       \
	defined(__HIP_PLATFORM_NVIDIA__) || /* hipcc targeting NVIDIA                      */                                       \
	defined(__HIPCC__) ||               /* hipcc generic                                */                                      \
	defined(__HIP__) ||                 /* clang --offload=hip                         */                                       \
	defined(THRUST_DEVICE_SYSTEM) ||    /* explicitly configured by CMake (e.g. NVIDIA thrust_create_target propagates this) */ \
	defined(BOOST_MULTI_HAS_THRUST) ||  /* explicit opt-in, e.g. set by Multi's CMake when find_package(Thrust) succeeds    */ \
	(defined(__has_include) &&          /* __has_include-based auto-detection:                                              */ \
	 __has_include(<thrust/version.h>))           /* backed by AMD ROCm                             */
#define BOOST_MULTI_ADL_HAS_THRUST 1
#endif

#ifdef BOOST_MULTI_ADL_HAS_THRUST
#if !defined(__NVCC__) && !defined(__CUDACC__) && !defined(__CUDA__) && \
    !defined(__HIPCC__) && !defined(__HIP_PLATFORM_AMD__) && !defined(__HIP_PLATFORM_NVIDIA__) && \
    !(defined(__has_include) && (__has_include(<cuda_runtime.h>) || __has_include(<hip/hip_runtime.h>)))
#if !defined(THRUST_DEVICE_SYSTEM)
#define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_CPP  // fallback: Thrust present but no GPU runtime detected
#endif
#endif

#ifdef __NVCC__
#pragma nv_diagnostic push
#pragma nv_diag_suppress = 20011  // deep inside Thrust: calling a __host__ function("std::vector<double, ::std::allocator<double> > ::vector(const ::std::vector<double, ::std::allocator<double> > &)") from a __host__ __device__ function("thrust::system::detail::generic::detail::uninitialized_copy_functor<    ::std::vector<double, ::std::allocator<double> > ,     ::std::vector<double, ::std::allocator<double> > > ::operator ()< ::thrust::detail::tuple_of_iterator_references<    ::std::vector<double, ::std::allocator<double> >  &,     ::std::vector<double, ::std::allocator<double> >  & > > ") is not allowed
#pragma nv_diag_suppress = 20014  // deep inside Thrust: calling a __host__ function from a __host__ __device__ function is not allowed
#pragma nv_diag_suppress = 20015  // deep inside Thrust: calling a constexpr __host__ function from a __host__ __device__ function is not allowed
#endif

#include <thrust/copy.h>
#include <thrust/detail/allocator/destroy_range.h>
#include <thrust/detail/memory_algorithms.h>
#include <thrust/equal.h>
#include <thrust/swap.h>
#include <thrust/uninitialized_copy.h>

#include <exception>  // for std::terminate, fixes a bug un Thrust 2

namespace thrust {
template<class... Ts> void thrust_involved(Ts const&... /*unused*/) {}
}  // namespace thrust

#if defined(__NVCC__) || defined(__HIP_PLATFORM_NVIDIA__) || defined(__HIP_PLATFORM_AMD__) || defined(__HIPCC__)
#if THRUST_VERSION >= 300102
#include <cuda/std/iterator>
#endif
#endif

#ifdef __NVCC__
#pragma nv_diagnostic pop  // nv_diagnostics pop
#endif

#endif

#include <algorithm>    // for for_each, copy_n, fill, fill_n, lexicographical_compare, swap_ranges  // IWYU pragma: keep  // bug in iwyu 0.18
#include <cstddef>      // for size_t
#include <functional>   // for equal_to
#include <iterator>     // for iterator_traits, distance, size
#include <memory>       // for allocator_traits, allocator, pointer_traits
#include <type_traits>  // for decay_t, enable_if_t, conditional_t, declval, is_pointer, true_type
#include <utility>      // for forward, addressof

#ifdef _MULTI_FORCE_TRIVIAL_STD_COMPLEX
#include <complex>
#endif

#define BOOST_MULTI_DEFINE_ADL(FuN) /*NOLINT(cppcoreguidelines-macro-usage) TODO(correaa) consider replacing for all ADL'd operations*/                                                                                                                                                                                                                                                                                                                                                                                                                            \
	namespace boost {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              \
	namespace multi {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              \
	namespace adl {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                \
	namespace custom {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             \
	template<class...> struct FuN##_t;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             \
	}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              \
	__attribute__((unused)) static constexpr class FuN##_t {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       \
		template<class... As> [[deprecated]] static auto _(priority<0>, As&&... args) = delete;                                                                                                                                                                                                                                                                                                                                                                                                                                                                    \
	 template<class... As> static auto                   _(priority<1>, As&&... args) BOOST_MULTI_DECLRET(std::FuN(std::forward<As>(args)...)) template<class... As> static auto _(priority<2>, As&&... args) BOOST_MULTI_DECLRET(FuN(std::forward<As>(args)...)) template<class T, class... As> static auto _(priority<3>, T&& t, As&&... args) BOOST_MULTI_DECLRET(std::forward<T>(t).FuN(std::forward<As>(args)...)) template<class... As> static auto _(priority<4>, As&&... args) BOOST_MULTI_DECLRET(custom::FuN##_t<As&&...>::_(std::forward<As>(args)...)) \
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   \
		 public : template<class... As>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            \
				  auto operator()(As&&... args) const -> decltype(_(priority<4>(), std::forward<As>(args)...)) { return _(priority<4>(), std::forward<As>(args)...); }                                                                                                                                                                                                                                                                                                                                                                                             \
	}(FuN);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        \
	} /* end namespace adl   */                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    \
	} /* end namespace multi */                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    \
	} /* end namespace boost */

#define BOOST_MULTI_DECLRET(ExpR) \
	->decltype(ExpR) { return ExpR; }  // NOLINT(cppcoreguidelines-macro-usage) saves a lot of typing
#define BOOST_MULTI_JUSTRET(ExpR) \
	{                             \
		return ExpR;              \
	}  // NOLINT(cppcoreguidelines-macro-usage) saves a lot of typing

namespace boost::multi {

template<std::size_t N> struct priority : std::conditional_t<N == 0, std::true_type, priority<N - 1>> {};

// clang-format off
class adl_copy_n_t {
	template<class... As>          constexpr static auto _(priority<0>/**/,          As&&... args) BOOST_MULTI_DECLRET(std::                copy_n(                      std::forward<As>(args)...))
#ifdef BOOST_MULTI_ADL_HAS_THRUST
	template<class... As>          constexpr static auto _(priority<1>/**/,          As&&... args) BOOST_MULTI_DECLRET(::thrust::           copy_n(                      std::forward<As>(args)...))
#endif
	template<class... As>          constexpr static auto _(priority<2>/**/,          As&&... args) BOOST_MULTI_DECLRET(                     copy_n(                      std::forward<As>(args)...))
	template<class T, class... As> constexpr static auto _(priority<3>/**/, T&& arg, As&&... args) BOOST_MULTI_DECLRET(std::decay_t<T>::    copy_n(std::forward<T>(arg), std::forward<As>(args)...))
	template<class T, class... As> constexpr static auto _(priority<4>/**/, T&& arg, As&&... args) BOOST_MULTI_DECLRET(std::forward<T>(arg).copy_n(                      std::forward<As>(args)...))

 public:
	template<class... As> constexpr auto operator()(As&&... args) const BOOST_MULTI_DECLRET(_(priority<4>(), std::forward<As>(args)...))
};
inline constexpr adl_copy_n_t adl_copy_n;
// clang-format on

// there is no move_n (std::move_n), use copy_n(std::make_move_iterator(first), count) instead

class adl_move_t {
 template<class... As> static constexpr auto _(priority<0> /**/, As&&... args) BOOST_MULTI_DECLRET(std::move(std::forward<As>(args)...))
#ifdef BOOST_MULTI_ADL_HAS_THRUST
	 template<class It, class... As>
	 static constexpr auto _(priority<1> /**/, It first, It last, As&&... args) BOOST_MULTI_DECLRET(thrust::copy(std::make_move_iterator(first), std::make_move_iterator(last), std::forward<As>(args)...))
#endif
		 template<class... As>
		 static constexpr auto _(priority<2> /**/, As&&... args) BOOST_MULTI_DECLRET(move(std::forward<As>(args)...)) template<class T, class... As> static constexpr auto _(priority<3> /**/, T&& arg, As&&... args) BOOST_MULTI_DECLRET(std::decay_t<T>::move(std::forward<T>(arg), std::forward<As>(args)...)) template<class T, class... As> static constexpr auto _(priority<4> /**/, T&& arg, As&&... args) BOOST_MULTI_DECLRET(std::forward<T>(arg).move(std::forward<As>(args)...))

			 public : template<class... As>
					  constexpr auto operator()(As&&... args) const BOOST_MULTI_DECLRET(_(priority<4>(), std::forward<As>(args)...))
};
inline constexpr adl_move_t adl_move;

class adl_fill_n_t {
 template<class... As> static constexpr auto _(priority<0> /**/, As&&... args) BOOST_MULTI_DECLRET(std::fill_n(std::forward<As>(args)...))
#ifdef BOOST_MULTI_ADL_HAS_THRUST
	 template<class... As>
	 static constexpr auto _(priority<1> /**/, As&&... args) BOOST_MULTI_DECLRET(thrust::fill_n(std::forward<As>(args)...))
#endif
		 template<class... As>
		 static constexpr auto _(priority<2> /**/, As&&... args) BOOST_MULTI_DECLRET(fill_n(std::forward<As>(args)...)) template<class T, class... As> static constexpr auto _(priority<3> /**/, T&& arg, As&&... args) BOOST_MULTI_DECLRET(std::decay_t<T>::fill_n(std::forward<T>(arg), std::forward<As>(args)...)) template<class T, class... As> static constexpr auto _(priority<4> /**/, T&& arg, As&&... args) BOOST_MULTI_DECLRET(std::forward<T>(arg).fill_n(std::forward<As>(args)...))

			 public : template<class... As>
					  constexpr auto operator()(As&&... args) const BOOST_MULTI_DECLRET(_(priority<4>(), std::forward<As>(args)...))
};
inline constexpr adl_fill_n_t adl_fill_n;

class adl_fill_t {
 template<class... As> static constexpr auto _(priority<0> /**/, As&&... args) BOOST_MULTI_DECLRET(std::fill(std::forward<As>(args)...))
#ifdef BOOST_MULTI_ADL_HAS_THRUST
	 template<class... As>
	 static constexpr auto _(priority<1> /**/, As&&... args) BOOST_MULTI_DECLRET(thrust::fill(std::forward<As>(args)...))
#endif
		 template<class... As>
		 static constexpr auto _(priority<2> /**/, As&&... args) BOOST_MULTI_DECLRET(fill(std::forward<As>(args)...)) template<class T, class... As> static constexpr auto _(priority<3> /**/, T&& arg, As&&... args) BOOST_MULTI_DECLRET(std::decay_t<T>::fill(std::forward<T>(arg), std::forward<As>(args)...)) template<class T, class... As> static constexpr auto _(priority<4> /**/, T&& arg, As&&... args) BOOST_MULTI_DECLRET(std::forward<T>(arg).fill(std::forward<As>(args)...))

			 public : template<class... As>
					  constexpr auto operator()(As&&... args) const BOOST_MULTI_DECLRET(_(priority<4>(), std::forward<As>(args)...))
};
inline constexpr adl_fill_t adl_fill;

class adl_equal_t {
 template<class... As> constexpr auto _(priority<1> /**/, As&&... args) const BOOST_MULTI_DECLRET(std::equal(std::forward<As>(args)...))
#ifdef BOOST_MULTI_ADL_HAS_THRUST
	 template<class... As>
	 constexpr auto _(priority<2> /**/, As&&... args) const BOOST_MULTI_DECLRET(::thrust::equal(std::forward<As>(args)...))
#endif
		 template<class... As>
		 constexpr auto _(priority<3> /**/, As&&... args) const BOOST_MULTI_DECLRET(equal(std::forward<As>(args)...)) template<class... As> constexpr auto _(priority<4> /**/, As&&... args) const BOOST_MULTI_DECLRET(equal(std::forward<As>(args)..., std::equal_to<>{}))  // WORKAROUND makes syntax compatible with boost::ranges::equal if, for some reason, it is included.
	 template<class T, class... As>
	 constexpr auto _(priority<5> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRET(std::decay_t<T>::equal(std::forward<T>(arg), std::forward<As>(args)...)) template<class T, class... As> constexpr auto _(priority<6> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRET(std::forward<T>(arg).equal(std::forward<As>(args)...))

		 public : template<class... As>
				  constexpr auto operator()(As&&... args) const BOOST_MULTI_DECLRET(_(priority<6>(), std::forward<As>(args)...))
};
inline constexpr adl_equal_t adl_equal;

#ifndef _MSC_VER
template<class... As, class = std::enable_if_t<sizeof...(As) == 0>> void copy(As...) = delete;  // NOLINT(modernize-use-constraints) TODO(correaa)
#endif

class adl_copy_t {
 template<class InputIt, class OutputIt, class = std::enable_if_t<std::is_assignable_v<typename std::iterator_traits<OutputIt>::reference, typename std::iterator_traits<InputIt>::reference>>  // NOLINT(modernize-use-constraints) TODO(correaa)
		  >
 constexpr auto _(priority<1> /**/, InputIt first, InputIt last, OutputIt d_first) const BOOST_MULTI_DECLRET(std::copy(first, last, d_first))  // cppcheck-suppress functionStatic ; TODO(correaa) consider making these functions static
#ifdef BOOST_MULTI_ADL_HAS_THRUST
	 template<class... As>
	 constexpr auto _(priority<2> /**/, As&&... args) const BOOST_MULTI_DECLRET(::thrust::copy(std::forward<As>(args)...))
#endif
		 template<class... As>
		 constexpr auto _(priority<3> /**/, As&&... args) const BOOST_MULTI_DECLRET(copy(std::forward<As>(args)...)) template<class T, class... As> constexpr auto _(priority<4> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRET(std::decay_t<T>::copy(std::forward<T>(arg), std::forward<As>(args)...))
	 //  template<class... As         > constexpr auto _(priority<5>/**/,          As&&... args) const BOOST_MULTI_DECLRET(boost::multi::adl_custom_copy<std::decay_t<As>...>::copy(std::forward<As>(as)...))
	 template<class T, class... As>
	 constexpr auto _(priority<6> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRET(std::forward<T>(arg).copy(std::forward<As>(args)...))

		 public : template<class... As>
				  constexpr auto operator()(As&&... args) const BOOST_MULTI_DECLRET(_(priority<6>(), std::forward<As>(args)...))
	// If constexpr std::copy/std::sort (C++20/26) hit a non-constexpr thrust::copy in the priority cascade,
	// route to std::copy during constant evaluation. Keep the trailing return type for SFINAE-friendliness.
	// public: template<class... As>
	// constexpr auto operator()(As&&... args) const -> decltype(_(priority<6>(), std::forward<As>(args)...)) {
	// #if defined(__cpp_lib_is_constant_evaluated) && (__cpp_lib_is_constant_evaluated >= 201811L)
	//   if(std::is_constant_evaluated()) { return std::copy(std::forward<As>(args)...); }
	// #endif
	//   return _(priority<6>(), std::forward<As>(args)...);
	// }
};
inline constexpr adl_copy_t adl_copy;

namespace adl {
// namespace custom {template<class...> struct fill_t;}
class fill_t {
 template<class... As> auto _(priority<1> /**/, As&&... args) const BOOST_MULTI_DECLRET(std::fill(std::forward<As>(args)...)) template<class... As> auto _(priority<2> /**/, As&&... args) const BOOST_MULTI_DECLRET(fill(std::forward<As>(args)...)) template<class T, class... As> auto _(priority<3> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRET(std::forward<T>(arg).fill(std::forward<As>(args)...))
	 // template<class... As>          auto _(priority<4>/**/,          As&&... args) const BOOST_MULTI_DECLRET(custom::             fill_t<As&&...>::_(std::forward<As>(args)...))

	 public : template<class... As>
			  auto operator()(As&&... args) const BOOST_MULTI_DECLRET(_(priority<5>{}, std::forward<As>(args)...))
};
inline constexpr fill_t fill;

class transform_t {
//  template<class InputIt, class OutputIt, class = std::enable_if_t<std::is_assignable_v<typename std::iterator_traits<OutputIt>::reference, typename std::iterator_traits<InputIt>::reference>>  // NOLINT(modernize-use-constraints) TODO(correaa)
// 		  >
//  constexpr auto _(priority<1> /**/, InputIt first, InputIt last, OutputIt d_first) const BOOST_MULTI_DECLRET(std::copy(first, last, d_first))  // cppcheck-suppress functionStatic ; TODO(correaa) consider making these functions static
	template<         class... As> static constexpr auto _(priority<1> /**/,          As&&... args) BOOST_MULTI_DECLRET(              ::std::transform(                      std::forward<As>(args)...))
#ifdef BOOST_MULTI_ADL_HAS_THRUST
	template<         class... As> static constexpr auto _(priority<2> /**/,          As&&... args) BOOST_MULTI_DECLRET(           ::thrust::transform(                      std::forward<As>(args)...))
#endif
	template<         class... As> static constexpr auto _(priority<3> /**/,          As&&... args) BOOST_MULTI_DECLRET(                     transform(                      std::forward<As>(args)...))
	template<class T, class... As> static constexpr auto _(priority<4> /**/, T&& arg, As&&... args) BOOST_MULTI_DECLRET(    std::decay_t<T>::transform(std::forward<T>(arg), std::forward<As>(args)...))
	template<class T, class... As> static constexpr auto _(priority<5> /**/, T&& arg, As&&... args) BOOST_MULTI_DECLRET(std::forward<T>(arg).transform(                      std::forward<As>(args)...))

 public : template<class... As>
	constexpr auto operator()(As&&... args) const BOOST_MULTI_DECLRET(_(priority<5>(), std::forward<As>(args)...))
};
inline constexpr transform_t transform;


}  // end namespace adl

namespace xtd {

template<class T>  // this one goes last!!!
constexpr auto to_address(T const& ptr) noexcept;

template<class T>
constexpr auto me_to_address(priority<0> /**/, T const& ptr) noexcept
	-> decltype(to_address(ptr.operator->())) {
	return to_address(ptr.operator->());
}

template<class T>
constexpr auto me_to_address(priority<1> /**/, T const& ptr) noexcept
	-> decltype(std::pointer_traits<T>::to_address(ptr)) {
	return std::pointer_traits<T>::to_address(ptr);
}

template<class T, std::enable_if_t<std::is_pointer<T>{}, int> = 0>  // NOLINT(modernize-use-constraints) TODO(correaa)
constexpr auto me_to_address(priority<2> /**/, T const& ptr) noexcept -> T {
	static_assert(!std::is_function_v<T>);
	return ptr;
}

template<class T>  // this one goes last!!!
constexpr auto to_address(T const& ptr) noexcept
	-> decltype(me_to_address(priority<2>(), ptr)) {
	return me_to_address(priority<2>(), ptr);
}

template<class Alloc, class ForwardIt, class Size, typename Value = typename std::iterator_traits<ForwardIt>::value_type, typename = decltype(std::addressof(*ForwardIt{})), typename = decltype(Value())>
auto alloc_uninitialized_value_construct_n(Alloc& alloc, ForwardIt first, Size count) -> ForwardIt {
	// ->std::decay_t<decltype(std::allocator_traits<Alloc>::construct(alloc, std::addressof(*first), Value()), first)>
	ForwardIt current = first;
	try {
		for(; count > 0; ++current, --count) {                                                  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
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

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
#endif

template<class Alloc, class ForwardIt, class Size, class T = typename std::iterator_traits<ForwardIt>::value_type>
auto alloc_uninitialized_default_construct_n(Alloc& alloc, ForwardIt first, Size count)
	-> std::decay_t<decltype(std::allocator_traits<Alloc>::construct(alloc, std::addressof(*first)), first)> {
	if(std::is_trivially_default_constructible_v<T>) {
		std::advance(first, count);
		return first;
	}
	using alloc_traits = std::allocator_traits<Alloc>;
	ForwardIt current  = first;

	try {
		//  return std::for_each_n(first, count, [&](T& elem) { alloc_traits::construct(alloc, std::addressof(elem)); ++current; });
		//  workadoung for gcc 8.3.1 in Lass
		std::for_each(first, first + count, [&](T& elem) -> void { alloc_traits::construct(alloc, std::addressof(elem)); ++current; });  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
		return first + count;                                   // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
	} catch(...) {
		// LCOV_EXCL_START  // TODO(correaa) add test
		std::for_each(first, current, [&](T& elem) -> void { alloc_traits::destroy(alloc, std::addressof(elem)); });
		throw;
		// LCOV_EXCL_STOP
	}

	// return current;
}

#ifdef __clang__
#pragma clang diagnostic pop
#endif

}  // end namespace xtd

template<class BidirIt, class Size, class T = typename std::iterator_traits<BidirIt>::value_type>
constexpr auto destroy_n(BidirIt first, Size count)
	-> std::decay_t<decltype(std::addressof(*(first - 1)), first)> {  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
	first += count;
	for(; count != 0; --first, --count) {    // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
		std::addressof(*(first - 1))->~T();  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
	}
	return first;
}

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
#endif

template<class Alloc, class BidirIt, class Size, class T = typename std::iterator_traits<BidirIt>::value_type>
constexpr auto alloc_destroy_n(Alloc& alloc, BidirIt first, Size count)
	-> std::decay_t<decltype(std::addressof(*(first - 1)), first)> {                 // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
	first += count;                                                                  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
	for(; count != 0; --first, --count) {                                            // NOLINT(altera-unroll-loops,cppcoreguidelines-pro-bounds-pointer-arithmetic) TODO(correaa) consider using an algorithm
		std::allocator_traits<Alloc>::destroy(alloc, std::addressof(*(first - 1)));  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
	}
	return first;
}

#ifdef __clang__
#pragma clang diagnostic pop
#endif

class adl_uninitialized_copy_t {
	template<class InIt, class FwdIt, class = decltype(std::addressof(*FwdIt{}))>                 // sfinae friendy std::uninitialized_copy
	[[nodiscard]] constexpr auto _(priority<1> /**/, InIt first, InIt last, FwdIt d_first) const  // N_O_L_I_N_T(performance-unnecessary-value-param) bug in clang-tidy
	// BOOST_MULTI_DECLRET(       std::uninitialized_copy(first, last, d_first))
	{
#if (__cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L) && !defined(__CUDA_ARCH__))
		using ValueType = typename std::iterator_traits<FwdIt>::value_type;
		if(
			std::is_constant_evaluated() && (std::is_trivially_default_constructible_v<ValueType> || multi::force_element_trivial_default_construction<ValueType>)
		) {
			return std::copy(std::move(first), std::move(last), std::move(d_first));
		} else  // NOLINT(llvm-else-after-return,readability-else-after-return,-warnings-as-errors)
#endif
		{
			return std::uninitialized_copy(std::move(first), std::move(last), std::move(d_first));
		}
	}
#ifdef BOOST_MULTI_ADL_HAS_THRUST
	template<class InIt, class FwdIt, class ValueType = typename std::iterator_traits<FwdIt>::value_type>
	constexpr auto _(priority<2> /**/, InIt first, InIt last, FwdIt d_first) const -> decltype(::thrust::uninitialized_copy(first, last, d_first))  // doesn't work with culang 17, cuda 12 ?
	{
		if constexpr(std::is_trivially_default_constructible_v<ValueType> || multi::force_element_trivial_default_construction<ValueType>) {
			return ::thrust::copy(first, last, d_first);
		} else {
			return ::thrust::uninitialized_copy(first, last, d_first);
		}
	}
#endif
 template<class TB, class... As> constexpr auto _(priority<3> /**/, TB&& first, As&&... args) const BOOST_MULTI_DECLRET(uninitialized_copy(std::forward<TB>(first), std::forward<As>(args)...)) template<class TB, class TE, class DB> constexpr auto _(priority<4> /**/, TB&& first, TE&& last, DB&& d_first) const BOOST_MULTI_DECLRET(std::decay_t<DB>::uninitialized_copy(std::forward<TB>(first), std::forward<TE>(last), std::forward<DB>(d_first))) template<class TB, class... As> constexpr auto _(priority<5> /**/, TB&& first, As&&... args) const BOOST_MULTI_DECLRET(std::decay_t<TB>::uninitialized_copy(std::forward<TB>(first), std::forward<As>(args)...)) template<class TB, class... As> constexpr auto _(priority<6> /**/, TB&& first, As&&... args) const BOOST_MULTI_DECLRET(std::forward<TB>(first).uninitialized_copy(std::forward<As>(args)...))

	 public : template<class... As>
			  constexpr auto operator()(As&&... args) const BOOST_MULTI_DECLRET(_(priority<6>(), std::forward<As>(args)...))
};
inline constexpr adl_uninitialized_copy_t adl_uninitialized_copy;

#ifdef __NVCC__
#pragma nv_diagnostic push
#pragma nv_diag_suppress = implicit_return_from_non_void_function
#endif

class adl_uninitialized_copy_n_t {
	template<class... As> constexpr auto _(priority<1> /**/, As&&... args) const BOOST_MULTI_DECLRET(std::uninitialized_copy_n(std::forward<As>(args)...)) template<class... As> constexpr auto _(priority<2> /**/, As&&... args) const BOOST_MULTI_DECLRET(uninitialized_copy_n(std::forward<As>(args)...))
#ifdef BOOST_MULTI_ADL_HAS_THRUST
	template<
		class It, class Size, class ItFwd,
		class ValueType = typename std::iterator_traits<ItFwd>::value_type,
		class = std::enable_if_t<! std::is_rvalue_reference_v<typename std::iterator_traits<It>::reference> >  // NOLINT(modernize-use-constraints) for C++20
	>
	constexpr auto _(priority<3>/**/, It first, Size count, ItFwd d_first) const -> decltype(::thrust::uninitialized_copy_n(first, count, d_first)) {  // NOLINT(performance-unnecessary-value-param)
		if constexpr(std::is_trivially_default_constructible_v<ValueType> || multi::force_element_trivial_default_construction<ValueType>) {
			return count?::thrust::copy_n(first, count, d_first):d_first;  // condition fixes a bug in Thrust 2, github.com/NVIDIA/thrust/issues/939, fixed later
		} else {
			return ::thrust::uninitialized_copy_n(first, count, d_first);
		}
	}
#endif
 template<class T, class... As> constexpr auto _(priority<4> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRET(std::decay_t<T>::uninitialized_copy_n(std::forward<T>(arg), std::forward<As>(args)...)) template<class T, class... As> constexpr auto _(priority<5> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRET(std::forward<T>(arg).uninitialized_copy_n(std::forward<As>(args)...))

	 public : template<class... As>
			  constexpr auto operator()(As&&... args) const BOOST_MULTI_DECLRET(_(priority<5>{}, std::forward<As>(args)...))  // TODO(correaa) this might trigger a compiler crash with g++ 7.5 because of operator&() && overloads
};
inline constexpr adl_uninitialized_copy_n_t adl_uninitialized_copy_n;

#ifdef __NVCC__
#pragma nv_diagnostic pop
#endif

class adl_uninitialized_move_n_t {
 template<class... As> constexpr auto _(priority<1> /**/, As&&... args) const BOOST_MULTI_DECLRET(std::uninitialized_move_n(std::forward<As>(args)...)) template<class... As> constexpr auto _(priority<2> /**/, As&&... args) const BOOST_MULTI_DECLRET(uninitialized_move_n(std::forward<As>(args)...)) template<class T, class... As> constexpr auto _(priority<3> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRET(std::decay_t<T>::uninitialized_move_n(std::forward<T>(arg), std::forward<As>(args)...)) template<class T, class... As> constexpr auto _(priority<4> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRET(std::forward<T>(arg).uninitialized_move_n(std::forward<As>(args)...))

	 public : template<class... As>
			  constexpr auto operator()(As&&... args) const { return _(priority<4>(), std::forward<As>(args)...); }
};
inline constexpr auto adl_uninitialized_move_n = adl_uninitialized_move_n_t{};

namespace xtd {

template<class T, class InputIt, class Size, class ForwardIt>
constexpr auto alloc_uninitialized_copy_n(std::allocator<T>& /*alloc*/, InputIt first, Size count, ForwardIt d_first) {
	// a std::allocator unambiguously denotes host memory, so never route through thrust's
	// __host__ __device__ path (which fails to compile for host-only types like std::string)
	return std::uninitialized_copy_n(first, count, d_first);
}

// template<class T, class InputIt, class Size, class ForwardIt>
// constexpr auto alloc_uninitialized_move_n(std::allocator<T>& /*alloc*/, InputIt first, Size count, ForwardIt d_first) {
//  return adl_uninitialized_move_n(first, count, d_first);}

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

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-warning-option"
#pragma clang diagnostic ignored "-Wunsafe-buffer-usage"  // TODO(correaa) use checked span
#endif

template<class Alloc, class InputIt, class Size, class ForwardIt>
// [[deprecated("check")]]
auto alloc_uninitialized_move_n(Alloc& alloc, InputIt first, Size count, ForwardIt d_first) {
	ForwardIt current = d_first;  // mull-ignore: cxx_init_const
	try {
		// NOLINTNEXTLINE(altera-unroll-loops,cppcoreguidelines-pro-bounds-pointer-arithmetic) TODO(correaa) consider using an algorithm
		for(; count > 0; ++first, ++current, --count) {  // mull-ignore: cxx_gt_to_ge
			std::allocator_traits<Alloc>::construct(alloc, std::addressof(*current), std::move(*first));
		}
		return current;
	} catch(...) {
		for(; d_first != current; ++d_first) {  // NOLINT(altera-unroll-loops,altera-id-dependent-backward-branch,cppcoreguidelines-pro-bounds-pointer-arithmetic) TODO(correaa) consider using an algorithm
			std::allocator_traits<Alloc>::destroy(alloc, std::addressof(*d_first));
		}
		throw;
	}
}

#ifdef __clang__
#pragma clang diagnostic pop
#endif

template<class T, class InputIt, class ForwardIt>
constexpr auto alloc_uninitialized_copy(std::allocator<T>& /*allocator*/, InputIt const& first, InputIt const& last, ForwardIt d_first) {
	return adl_uninitialized_copy(first, last, d_first);
}

template<class Alloc, class InputIt, class ForwardIt, class = decltype(std::addressof(*std::declval<ForwardIt>())), class = std::enable_if_t<std::is_constructible_v<typename std::iterator_traits<ForwardIt>::value_type, typename std::iterator_traits<InputIt>::reference>>  // NOLINT(modernize-use-constraints) TODO(correaa)
		 >
#if (__cplusplus >= 202002L || (defined(_MSVC_LANG) && _MSVC_LANG >= 202002L))
constexpr
#endif
	auto alloc_uninitialized_copy(Alloc& alloc, InputIt first, InputIt last, ForwardIt d_first) {
	// ->std::decay_t<decltype(a.construct(std::addressof(*d_first), *first), d_first)> // problematic in clang-11 + gcc-9
	ForwardIt current  = d_first;
	using alloc_traits = std::allocator_traits<Alloc>;
	try {
		std::for_each(first, last, [&](auto const& elem) -> void {  // TODO(correaa) replace by adl_for_each
			alloc_traits::construct(alloc, std::addressof(*current), elem);
			++current;
		});
		return current;
	} catch(...) {
		std::for_each(d_first, current, [&](auto const& elem) -> void {
			std::allocator_traits<Alloc>::destroy(alloc, std::addressof(elem));
		});
		throw;
	}
}

template<class Alloc, class ForwardIt, class Size, class T>
auto alloc_uninitialized_fill_n(Alloc& alloc, ForwardIt first, Size n, T const& value)
	-> std::decay_t<decltype(std::allocator_traits<Alloc>::construct(alloc, std::addressof(*first), value), first)> {
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

namespace adl {

#ifdef BOOST_MULTI_ADL_HAS_THRUST
#if THRUST_VERSION < 300102
using ::thrust::distance;
#else
using ::cuda::std::distance;
#endif
#endif

class adl_distance_t {
 template<class... As> constexpr auto _(priority<1> /**/, As&&... args) const BOOST_MULTI_DECLRET(std::distance(std::forward<As>(args)...)) template<class... As> constexpr auto _(priority<2> /**/, As&&... args) const BOOST_MULTI_DECLRET(distance(std::forward<As>(args)...))
	 // template<class It1, class It2> constexpr auto _(priority<3>/**/, It1 it1, It2 it2     ) const BOOST_MULTI_DECLRET(it2 - it1)
	 template<class T, class... As>
	 constexpr auto _(priority<4> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRET(std::decay_t<T>::distance(std::forward<T>(arg), std::forward<As>(args)...)) template<class T, class... As> constexpr auto _(priority<5> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRET(std::forward<T>(arg).distance(std::forward<As>(args)...))

		 public : template<class... As>
				  constexpr auto operator()(As&&... args) const BOOST_MULTI_DECLRET(_(priority<5>{}, std::forward<As>(args)...))
};
}  // end namespace adl

inline constexpr adl::adl_distance_t adl_distance;

class adl_begin_t {
 template<class... As> constexpr auto _(priority<1> /**/, As&&... args) const BOOST_MULTI_DECLRET(std::begin(std::forward<As>(args)...))
	 //  template<class... As>          constexpr auto _(priority<2>/**/,          As&&... args) const BOOST_MULTI_DECLRET(                     begin(std::forward<As>(args)...))  // this is catching boost::range_iterator if Boost 1.53 is included
	 // #if defined(__NVCC__)  // this is no thrust::begin
	 //  template<class... As>          constexpr auto _(priority<2>/**/,          As&&... args) const BOOST_MULTI_DECLRET(::thrust::           begin(                      std::forward<As>(args)...))
	 // #endif
	 template<class T, class... As>
	 constexpr auto _(priority<3> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRET(std::decay_t<T>::begin(std::forward<T>(arg), std::forward<As>(args)...)) template<class T, class... As> constexpr auto _(priority<4> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRET(std::forward<T>(arg).begin(std::forward<As>(args)...))

		 public : template<class... As>
				  [[nodiscard]] constexpr auto operator()(As&&... args) const BOOST_MULTI_DECLRET(_(priority<4>{}, std::forward<As>(args)...))
};
inline constexpr adl_begin_t adl_begin;

class adl_end_t {
 template<class... As> constexpr auto _(priority<1> /**/, As&&... args) const BOOST_MULTI_DECLRET(std::end(std::forward<As>(args)...))
	 // template<class... As>          constexpr auto _(priority<2>/**/,          As&&... args) const BOOST_MULTI_DECLRET(                     end(std::forward<As>(args)...))
	 // #if defined(__NVCC__)  // there is no thrust::end
	 //  template<class... As>          constexpr auto _(priority<2>/**/,          As&&... args) const BOOST_MULTI_DECLRET(::thrust::           end(                      std::forward<As>(args)...))
	 // #endif
	 template<class T, class... As>
	 constexpr auto _(priority<3> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRET(std::decay_t<T>::end(std::forward<T>(arg), std::forward<As>(args)...)) template<class T, class... As> constexpr auto _(priority<4> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRET(std::forward<T>(arg).end(std::forward<As>(args)...))

		 public : template<class... As>
				  constexpr auto operator()(As&&... args) const BOOST_MULTI_DECLRET(_(priority<4>{}, std::forward<As>(args)...))
};
inline constexpr adl_end_t adl_end;

class adl_size_t {
 template<class... As> constexpr auto _(priority<1> /**/, As&&... args) const BOOST_MULTI_DECLRET(std::size(std::forward<As>(args)...)) template<class... As> constexpr auto _(priority<2> /**/, As&&... args) const BOOST_MULTI_DECLRET(size(std::forward<As>(args)...)) template<class T, class... As> constexpr auto _(priority<3> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRET(std::decay_t<T>::size(std::forward<T>(arg), std::forward<As>(args)...)) template<class T, class... As> constexpr auto _(priority<4> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRET(std::forward<T>(arg).size(std::forward<As>(args)...))

	 public : template<class... As>
			  [[nodiscard]] constexpr auto operator()(As&&... args) const BOOST_MULTI_DECLRET(_(priority<4>{}, std::forward<As>(args)...))
};
inline constexpr adl_size_t adl_size;

class adl_swap_ranges_t {
 template<class... As> constexpr auto _(priority<1> /**/, As&&... args) const BOOST_MULTI_DECLRET(std::swap_ranges(std::forward<As>(args)...)) template<class... As> constexpr auto _(priority<2> /**/, As&&... args) const BOOST_MULTI_DECLRET(swap_ranges(std::forward<As>(args)...)) template<class T, class... As> constexpr auto _(priority<3> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRET(std::decay_t<T>::swap_ranges(std::forward<T>(arg), std::forward<As>(args)...)) template<class T, class... As> constexpr auto _(priority<4> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRET(std::forward<T>(arg).swap_ranges(std::forward<As>(args)...))

	 public : template<class... As>
			  constexpr auto operator()(As&&... args) const BOOST_MULTI_DECLRET(_(priority<4>{}, std::forward<As>(args)...))
};
inline constexpr adl_swap_ranges_t adl_swap_ranges;

class adl_lexicographical_compare_t {
 template<class... As> constexpr auto _(priority<1> /**/, As&&... args) const BOOST_MULTI_DECLRET(std::lexicographical_compare(std::forward<As>(args)...)) template<class... As> constexpr auto _(priority<2> /**/, As&&... args) const BOOST_MULTI_DECLRET(lexicographical_compare(std::forward<As>(args)...)) template<class T, class... As> constexpr auto _(priority<3> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRET(std::decay_t<T>::lexicographical_compare(std::forward<T>(arg), std::forward<As>(args)...)) template<class T, class... As> constexpr auto _(priority<4> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRET(std::forward<T>(arg).lexicographical_compare(std::forward<As>(args)...))

	 public : template<class... As>
			  constexpr auto operator()(As&&... args) const BOOST_MULTI_DECLRET(_(priority<4>{}, std::forward<As>(args)...))
};
inline constexpr adl_lexicographical_compare_t adl_lexicographical_compare;

class adl_uninitialized_value_construct_n_t {
 template<class... As> static constexpr auto _(priority<1> /**/, As&&... args) BOOST_MULTI_DECLRET(std::uninitialized_value_construct_n(std::forward<As>(args)...))  // TODO(correaa) use boost alloc_X functions?
	 template<class... As>
	 static constexpr auto _(priority<2> /**/, As&&... args) BOOST_MULTI_DECLRET(uninitialized_value_construct_n(std::forward<As>(args)...)) template<class T, class... As> static constexpr auto _(priority<3> /**/, T&& arg, As&&... args) BOOST_MULTI_DECLRET(std::decay_t<T>::uninitialized_value_construct_n(std::forward<T>(arg), std::forward<As>(args)...)) template<class T, class... As> static constexpr auto _(priority<4> /**/, T&& arg, As&&... args) BOOST_MULTI_DECLRET(std::forward<T>(arg).uninitialized_value_construct_n(std::forward<As>(args)...))

		 public : template<class... As>
				  constexpr auto operator()(As&&... args) const BOOST_MULTI_DECLRET(_(priority<4>(), std::forward<As>(args)...))
};
inline constexpr adl_uninitialized_value_construct_n_t adl_uninitialized_value_construct_n;

class adl_alloc_uninitialized_value_construct_n_t {
 template<class Alloc, class... As> constexpr static auto _(priority<1> /**/, Alloc&& /*alloc*/, As&&... args) BOOST_MULTI_DECLRET(adl_uninitialized_value_construct_n(std::forward<As>(args)...))  // NOLINT(cppcoreguidelines-missing-std-forward)
																																																	//  template<class... As>              constexpr static auto _(priority<2>/**/,                    As&&... args) BOOST_MULTI_DECLRET(              xtd::  alloc_uninitialized_value_construct_n(std::forward<As>(args)...))  // TODO(correaa) use boost alloc_X functions?
	 template<class... As>
	 constexpr static auto _(priority<3> /**/, As&&... args) BOOST_MULTI_DECLRET(alloc_uninitialized_value_construct_n(std::forward<As>(args)...)) template<class T, class... As> constexpr static auto _(priority<4> /**/, T&& arg, As&&... args) BOOST_MULTI_DECLRET(std::decay_t<T>::alloc_uninitialized_value_construct_n(std::forward<T>(arg), std::forward<As>(args)...)) template<class T, class... As> constexpr static auto _(priority<5> /**/, T&& arg, As&&... args) BOOST_MULTI_DECLRET(std::forward<T>(arg).alloc_uninitialized_value_construct_n(std::forward<As>(args)...))

		 public : template<class... As>
				  constexpr auto operator()(As&&... args) const BOOST_MULTI_DECLRET(_(priority<5>(), std::forward<As>(args)...))
};
inline constexpr adl_alloc_uninitialized_value_construct_n_t adl_alloc_uninitialized_value_construct_n;

class adl_uninitialized_default_construct_n_t {
	template<class... As> constexpr static auto _(priority<1> /**/, As&&... args) { return std::uninitialized_default_construct_n(std::forward<As>(args)...); }
	// #if defined(__NVCC__)
	// template<class... As>          constexpr auto _(priority<2>/**/,          As&&... args) const BOOST_MULTI_DECLRET(             thrust::uninitialized_default_construct_n(                      std::forward<As>(args)...))
	// #endif
 template<class... As> constexpr static auto _(priority<3> /**/, As&&... args) BOOST_MULTI_DECLRET(uninitialized_default_construct_n(std::forward<As>(args)...)) template<class T, class... As> constexpr static auto _(priority<4> /**/, T&& arg, As&&... args) BOOST_MULTI_DECLRET(std::decay_t<T>::uninitialized_default_construct_n(std::forward<T>(arg), std::forward<As>(args)...)) template<class T, class... As> constexpr static auto _(priority<5> /**/, T&& arg, As&&... args) BOOST_MULTI_DECLRET(std::forward<T>(arg).uninitialized_default_construct_n(std::forward<As>(args)...))

	 public : template<class... As>
			  constexpr auto operator()(As&&... args) const BOOST_MULTI_DECLRET(_(priority<5>(), std::forward<As>(args)...))
};
inline constexpr adl_uninitialized_default_construct_n_t adl_uninitialized_default_construct_n;

class adl_alloc_uninitialized_default_construct_n_t {
 template<class Alloc, class... As> constexpr static auto _(priority<1> /**/, Alloc&& /*unused*/, As&&... args) BOOST_MULTI_JUSTRET(adl_uninitialized_default_construct_n(std::forward<As>(args)...))  // NOLINT(cppcoreguidelines-missing-std-forward)
	 template<class... As>
	 constexpr static auto _(priority<2> /**/, As&&... args) BOOST_MULTI_DECLRET(xtd::alloc_uninitialized_default_construct_n(std::forward<As>(args)...))  // TODO(correaa) use boost alloc_X functions?
#ifdef BOOST_MULTI_ADL_HAS_THRUST
#if defined(THRUST_VERSION) && (THRUST_VERSION < 200700)  // 200800)
	 // boost::multi::detail::what_value_t<THRUST_VERSION> a;
	 template<class Alloc, class It, class Size>
	 constexpr static auto _(priority<3> /**/, Alloc&& alloc, It first, Size n) BOOST_MULTI_DECLRET((thrust::detail::default_construct_range(std::forward<Alloc>(alloc), first, n)))
#else
	 // boost::multi::detail::what_value_t<THRUST_VERSION> b;
	 template<class Alloc, class It, class Size>
	 constexpr static auto _(priority<3> /**/, Alloc&& alloc, It first, Size n) BOOST_MULTI_DECLRET((thrust::detail::value_initialize_range(std::forward<Alloc>(alloc), first, n)))
#endif
#endif
		 template<class... As>
		 constexpr static auto _(priority<4> /**/, As&&... args) BOOST_MULTI_DECLRET(alloc_uninitialized_default_construct_n(std::forward<As>(args)...)) template<class T, class... As> constexpr static auto _(priority<5> /**/, T&& arg, As&&... args) BOOST_MULTI_DECLRET(std::decay_t<T>::alloc_uninitialized_default_construct_n(std::forward<T>(arg), std::forward<As>(args)...)) template<class T, class... As> constexpr static auto _(priority<6> /**/, T&& arg, As&&... args) BOOST_MULTI_DECLRET(std::forward<T>(arg).alloc_uninitialized_default_construct_n(std::forward<As>(args)...))

			 public : template<class... As>
					  constexpr auto operator()(As&&... args) const BOOST_MULTI_DECLRET(_(priority<6>(), std::forward<As>(args)...))
};
inline constexpr adl_alloc_uninitialized_default_construct_n_t adl_alloc_uninitialized_default_construct_n;

class adl_destroy_n_t {
 template<class... As> constexpr static auto _(priority<1> /**/, As&&... args) BOOST_MULTI_DECLRET(multi::destroy_n(std::forward<As>(args)...)) template<class... As> constexpr static auto _(priority<2> /**/, As&&... args) BOOST_MULTI_DECLRET(destroy_n(std::forward<As>(args)...)) template<class T, class... As> constexpr static auto _(priority<3> /**/, T&& arg, As&&... args) BOOST_MULTI_DECLRET(std::decay_t<T>::destroy_n(std::forward<T>(arg), std::forward<As>(args)...)) template<class T, class... As> constexpr static auto _(priority<4> /**/, T&& arg, As&&... args) BOOST_MULTI_DECLRET(std::forward<T>(arg).destroy_n(std::forward<As>(args)...))

	 public : template<class... As>
			  constexpr auto operator()(As&&... args) const BOOST_MULTI_DECLRET(_(priority<4>(), std::forward<As>(args)...))
};
inline constexpr adl_destroy_n_t adl_destroy_n;

class adl_alloc_destroy_n_t {
	template<class Alloc, class... As> constexpr auto _(priority<1> /**/, Alloc&& /*unused*/, As&&... args) const BOOST_MULTI_DECLRET(adl_destroy_n(std::forward<As>(args)...))  // NOLINT(cppcoreguidelines-missing-std-forward)
	#ifdef BOOST_MULTI_ADL_HAS_THRUST  // CLANG_FORMAT_OFF
	template<class Alloc, class It, class Size>
	constexpr auto _(priority<2> /**/, Alloc& alloc, It first, Size n) const BOOST_MULTI_DECLRET((thrust::detail::destroy_range(alloc, first, first + n)))  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
	#endif  // DO_NOT_FORMAT
	template<class... As>
	constexpr auto _(priority<3> /**/, As&&... args) const BOOST_MULTI_DECLRET(multi::alloc_destroy_n(std::forward<As>(args)...))  // TODO(correaa) use boost alloc_X functions?
	template<class... As>
	constexpr auto _(priority<4> /**/, As&&... args) const BOOST_MULTI_DECLRET(alloc_destroy_n(std::forward<As>(args)...)) template<class T, class... As> constexpr auto _(priority<5> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRET(std::decay_t<T>::alloc_destroy_n(std::forward<T>(arg), std::forward<As>(args)...)) template<class T, class... As> constexpr auto _(priority<6> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRET(std::forward<T>(arg).alloc_destroy_n(std::forward<As>(args)...))

 public:
 	template<class... As>
	constexpr auto operator()(As&&... args) const BOOST_MULTI_DECLRET(_(priority<6>(), std::forward<As>(args)...))
};
inline constexpr adl_alloc_destroy_n_t adl_alloc_destroy_n;

class adl_alloc_uninitialized_copy_t {
 template<class Alloc, class... As> constexpr auto _(priority<1> /**/, Alloc&& /*ll*/, As&&... args) const BOOST_MULTI_DECLRET(adl_uninitialized_copy(std::forward<As>(args)...))  // NOLINT(cppcoreguidelines-missing-std-forward)
	 template<class Alloc, class... As>
	 constexpr auto _(priority<2> /**/, Alloc&& alloc, As&&... args) const BOOST_MULTI_DECLRET(xtd::alloc_uninitialized_copy(std::forward<Alloc>(alloc), std::forward<As>(args)...)) template<class Alloc, class... As> constexpr auto _(priority<3> /**/, Alloc&& alloc, As&&... args) const BOOST_MULTI_DECLRET(alloc_uninitialized_copy(std::forward<Alloc>(alloc), std::forward<As>(args)...)) template<class Alloc, class... As> constexpr auto _(priority<4> /**/, Alloc&& alloc, As&&... args) const BOOST_MULTI_DECLRET(std::decay_t<Alloc>::alloc_uninitialized_copy(std::forward<Alloc>(alloc), std::forward<As>(args)...)) template<class Alloc, class... As> constexpr auto _(priority<5> /**/, Alloc&& alloc, As&&... args) const BOOST_MULTI_DECLRET(std::forward<Alloc>(alloc).alloc_uninitialized_copy(std::forward<As>(args)...))

		 public : template<class... As>
				  constexpr auto operator()(As&&... args) const BOOST_MULTI_DECLRET(_(priority<5>{}, std::forward<As>(args)...))
};
inline constexpr adl_alloc_uninitialized_copy_t adl_alloc_uninitialized_copy;

class adl_alloc_uninitialized_copy_n_t {
 template<class Alloc, class... As> constexpr auto _(priority<1> /**/, Alloc&& /*alloc*/, As&&... args) const BOOST_MULTI_DECLRET(adl_uninitialized_copy_n(std::forward<As>(args)...))  // NOLINT(cppcoreguidelines-missing-std-forward)
	 template<class... As>
	 constexpr auto _(priority<2> /**/, As&&... args) const BOOST_MULTI_DECLRET(alloc_uninitialized_copy_n(std::forward<As>(args)...))
	 //  template<class... As>              constexpr auto _(priority<3>/**/,                    As&&... args) const BOOST_MULTI_DECLRET(                xtd::alloc_uninitialized_copy_n(std::forward<As>(args)...))
	 // #if defined(__NVCC__)
	 //  there is no thrust alloc uninitialized copy
	 // #endif
	 template<class T, class... As>
	 constexpr auto _(priority<5> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRET(std::decay_t<T>::alloc_uninitialized_copy_n(std::forward<T>(arg), std::forward<As>(args)...)) template<class T, class... As> constexpr auto _(priority<6> /**/, T&& arg, As&&... args) const BOOST_MULTI_DECLRET(std::forward<T>(arg).alloc_uninitialized_copy_n(std::forward<As>(args)...))

		 public : template<class... As>
				  constexpr auto operator()(As&&... args) const { return _(priority<6>(), std::forward<As>(args)...); }
};
inline constexpr adl_alloc_uninitialized_copy_n_t adl_alloc_uninitialized_copy_n;

class alloc_uninitialized_move_n_t {
	// TODO(correaa) : fallback to no alloc version
	template<class... As> static constexpr auto _(priority<1> /**/, As&&... args) { return xtd::alloc_uninitialized_move_n(std::forward<As>(args)...); }
	template<class... As> static constexpr auto _(priority<2> /**/, As&&... args) BOOST_MULTI_DECLRET(alloc_uninitialized_move_n(std::forward<As>(args)...)) template<class T, class... As> static constexpr auto _(priority<3> /**/, T&& arg, As&&... args) BOOST_MULTI_DECLRET(std::forward<T>(arg).alloc_uninitialized_move_n(std::forward<As>(args)...))

	 public : template<class... As>
			  constexpr auto operator()(As&&... args) const { return _(priority<3>(), std::forward<As>(args)...); }
};
inline constexpr alloc_uninitialized_move_n_t adl_alloc_uninitialized_move_n;

class uninitialized_fill_n_t {
	template<class... As>          static constexpr auto _(priority<1>/**/,          As&&... args) BOOST_MULTI_DECLRET((              std::   uninitialized_fill_n(std::forward<As>(args)...)))
	template<class... As>          static constexpr auto _(priority<2>/**/,          As&&... args) BOOST_MULTI_DECLRET((                      uninitialized_fill_n(std::forward<As>(args)...)))
#ifdef BOOST_MULTI_ADL_HAS_THRUST
	template<class... As>          static constexpr auto _(priority<3>/**/,          As&&... args) BOOST_MULTI_DECLRET(( thrust_involved(args...),::thrust::uninitialized_fill_n(std::forward<As>(args)...)))
#endif
	template<class T, class... As> static constexpr auto _(priority<4>/**/, T&& arg, As&&... args) BOOST_MULTI_DECLRET(( std::forward<T>(arg).uninitialized_fill_n(std::forward<As>(args)...)))

			 public : template<class T1, class... As>
					  constexpr auto operator()(T1&& arg, As&&... args) const BOOST_MULTI_DECLRET(_(priority<4>(), std::forward<T1>(arg), std::forward<As>(args)...))
};
inline constexpr uninitialized_fill_n_t adl_uninitialized_fill_n;

class alloc_uninitialized_fill_n_t {
 template<class... As> constexpr auto _(priority<1> /**/, As&&... args) const BOOST_MULTI_DECLRET(xtd::alloc_uninitialized_fill_n(std::forward<As>(args)...)) template<class Alloc, class... As> constexpr auto _(priority<2> /**/, Alloc&& /*alloc*/, As&&... args) const BOOST_MULTI_DECLRET(adl_uninitialized_fill_n(std::forward<As>(args)...))  // NOLINT(cppcoreguidelines-missing-std-forward)
	 template<class... As>
	 constexpr auto _(priority<3> /**/, As&&... args) const BOOST_MULTI_DECLRET(alloc_uninitialized_fill_n(std::forward<As>(args)...)) template<class Alloc, class... As> constexpr auto _(priority<4> /**/, Alloc&& alloc, As&&... args) const BOOST_MULTI_DECLRET(std::forward<Alloc>(alloc).alloc_uninitialized_fill_n(std::forward<As>(args)...))

		 public : template<class T1, class... As>
				  constexpr auto operator()(T1&& arg, As&&... args) const BOOST_MULTI_DECLRET(_(priority<4>{}, std::forward<T1>(arg), std::forward<As>(args)...))
};
inline constexpr alloc_uninitialized_fill_n_t adl_alloc_uninitialized_fill_n;

// template<dimensionality_type N>
// struct recursive {
//  template<class Alloc, class InputIt, class ForwardIt>
//  static constexpr auto alloc_uninitialized_copy(Alloc& alloc, InputIt first, InputIt last, ForwardIt dest){
//      using std::begin; using std::end;
//      while(first!=last) {  // NOLINT(altera-unroll-loops) TODO(correaa) consider using an algorithm
//          recursive<N-1>::alloc_uninitialized_copy(alloc, begin(*first), end(*first), begin(*dest));
//          ++first;
//          ++dest;
//      }
//      return dest;
//  }
// };

// template<> struct recursive<1> {
//  template<class Alloc, class InputIt, class ForwardIt>
//  static auto alloc_uninitialized_copy(Alloc& alloc, InputIt first, InputIt last, ForwardIt dest){
//      return adl_alloc_uninitialized_copy(alloc, first, last, dest);
//  }
// };

}  // end namespace boost::multi

#undef BOOST_MULTI_DECLRET
#undef BOOST_MULTI_JUSTRET

#endif  // BOOST_MULTI_DETAIL_ADL_HPP
