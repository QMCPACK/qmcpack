// Copyright 2021-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_THRUST_HPP
#define BOOST_MULTI_ADAPTORS_THRUST_HPP
// #pragma once

#include <thrust/detail/raw_reference_cast.h>

#include "boost/multi/array.hpp"

namespace thrust {

template<class T, boost::multi::dimensionality_type D, class ElementPtr, class Layout>
__host__ __device__
boost::multi::subarray<T, D, ElementPtr, Layout> const&
raw_reference_cast(boost::multi::subarray<T, D, ElementPtr, Layout> const& ref) { return ref; }

template<class T, boost::multi::dimensionality_type D, class ElementPtr, class Layout>
__host__ __device__
boost::multi::subarray<T, D, ElementPtr, Layout>&
raw_reference_cast(boost::multi::subarray<T, D, ElementPtr, Layout>& ref) { return ref; }

template<class T, boost::multi::dimensionality_type D, class ElementPtr, class Layout>
__host__ __device__
boost::multi::subarray<T, D, ElementPtr, Layout>&
raw_reference_cast(boost::multi::subarray<T, D, ElementPtr, Layout>&& ref) { return ref; }

template<class T, boost::multi::dimensionality_type D, class ElementPtr, class Layout>
__host__ __device__
boost::multi::const_subarray<T, D, ElementPtr, Layout> const&
raw_reference_cast(boost::multi::const_subarray<T, D, ElementPtr, Layout> const& ref) { return ref; }

}  // namespace thrust


#include <exception>

#include <thrust/device_allocator.h>
#include <thrust/universal_allocator.h>

#if !defined(MULTI_USE_HIP)
#include <thrust/system/cuda/memory.h>  // for ::thrust::cuda::allocator
#else
#include <thrust/system/hip/memory.h>  // for ::thrust::hip::allocator
#endif
// #include <thrust/system/cuda/memory.h>  // ::thrust::cuda::allocator

// #include <thrust/detail/type_traits/pointer_traits.h>

// #include <utility>  // std::copy

#include <thrust/detail/pointer.h>      // for pointer
#include <thrust/mr/allocator.h>        // for allocator (ptr only), stateless_resource_allocator
#include <thrust/mr/memory_resource.h>  // for memory_resource

#if !defined(MULTI_USE_HIP)
#include <cuda_runtime_api.h>                            // for cudaGetDevice, cudaMemPrefetchAsync, cudaPointerGetAttributes
#include <driver_types.h>                                // for cudaErrorInvalidValue, cudaPointerAttributes, cudaSuccess, cudaErrorInvalidDevice, cudaMemoryTypeManaged
#include <thrust/system/cuda/detail/execution_policy.h>  // for tag
#include <thrust/system/cuda/memory_resource.h>          // for universal_memory_resource
#include <thrust/system/cuda/pointer.h>                  // for universal_pointer

#if THRUST_VERSION >= 300200  // CCCL 2
#include <thrust/detail/type_traits/pointer_traits.h>  // needed for specialization for Thrust 3 (CUDA 13.2)
#endif

#else
// #include <thrust/system/hip/detail/execution_policy.h>        // for tag
// #include <thrust/system/hip/memory_resource.h>                // for universal_memory_resource
// #include <thrust/system/hip/pointer.h>                        // for universal_pointer

// #include <hip_runtime_api.h>                                  // for cudaGetDevice, cudaMemPrefetchAsync, cudaPointerGetAttributes
#endif

#include "boost/multi/adaptors/thrust/fix_pointer_traits.hpp"

#include <cassert>
#include <iterator>  // for iterator_traits
#include <memory>    // for allocator_traits, allocator, pointer_traits
// #include <thrust/iterator/detail/iterator_traits.inl>          // for iterator_system
#include <type_traits>  // for decay_t

// // begin of nvcc trhust 11.5 workaround : https://github.com/NVIDIA/thrust/issues/1629
// namespace thrust {

// template<typename Element, typename Tag, typename Reference, typename Derived> class pointer;
// template<class T> struct pointer_traits;

// }  // end namespace thrust

// namespace std {

// template<class... As> struct pointer_traits<thrust::pointer<As...>>
// : thrust::detail::pointer_traits<thrust::pointer<As...>> {
//  template<class T>
//  using rebind = typename thrust::detail::pointer_traits<thrust::pointer<As...>>::template rebind<T>::other;
// };

// }  // end namespace std
// // end of nvcc thrust 11.5 workaround

#if !defined(MULTI_USE_HIP)
#define HICUP cuda
#define HICUP_(NAME) cuda##NAME
#else
#define HICUP hip
#define HICUP_(NAME) hip##NAME
#endif

namespace boost::multi {
template<class Alloc> struct allocator_traits;
}

namespace boost::multi {

template<class T>
struct pointer_traits<::thrust::pointer<T, ::thrust::HICUP::tag, T&>> : std::pointer_traits<::thrust::pointer<T, ::thrust::HICUP::tag, T&>> {
	using default_allocator_type = ::thrust::universal_allocator<std::decay_t<T>>;
};

}  // end namespace boost::multi

namespace boost::multi {

template<class TT>
struct allocator_traits<::thrust::mr::stateless_resource_allocator<TT, ::thrust::system::HICUP::universal_memory_resource>>
: std::allocator_traits<::thrust::mr::stateless_resource_allocator<TT, ::thrust::system::HICUP::universal_memory_resource>> {
 private:
	using Alloc = ::thrust::mr::stateless_resource_allocator<TT, ::thrust::system::HICUP::universal_memory_resource>;
	using base  = std::allocator_traits<Alloc>;

 public:
	using typename base::const_void_pointer;
	using typename base::pointer;
	using typename base::size_type;

	using base::allocate;
	[[nodiscard]] static constexpr auto allocate(Alloc& alloc, size_type n, const_void_pointer hint) -> pointer {
		auto ret = allocator_traits::allocate(alloc, n);
		if(!hint) {
			prefetch_to_device_(ret, n * sizeof(TT), get_current_device_());
			return ret;
		}
		prefetch_to_device_(ret, n * sizeof(TT), get_device_(hint));
		return ret;
	}

 private:
	using device_index = int;
	static auto get_current_device_() -> device_index {
		int device;  // NOLINT(cppcoreguidelines-init-variables) delayed init
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wswitch-enum"
#endif
		switch(HICUP_(GetDevice)(&device)) {
		case HICUP_(Success): break;
		case HICUP_(ErrorInvalidValue): assert(0);  // NOLINT(bugprone-branch-clone)
		// ... 125 cases not handled
		default: assert(0);
		}
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif
		return device;
	}
	static void prefetch_to_device_(const_void_pointer ptr, size_type byte_count, device_index dev) {
#if(CUDART_VERSION < 13000)  // CudaMemPrefetchAsync changes its interface on version 13 TODO(correaa) update API call
		switch(HICUP_(MemPrefetchAsync)(raw_pointer_cast(ptr), byte_count, dev)) {
		case HICUP_(Success): break;
		case HICUP_(ErrorInvalidValue): assert(0); break;   // NOLINT(bugprone-branch-clone)
		case HICUP_(ErrorInvalidDevice): assert(0); break;  // NOLINT(bugprone-branch-clone)
		default: assert(0);
		}
#endif
	}

	static auto get_device_(const_void_pointer ptr) -> device_index {
#if defined(__HIPCC__)
		hipPointerAttribute_t attr{};
#else  // #if defined(__NVCC__)
		cudaPointerAttributes attr{};
#endif
		switch(HICUP_(PointerGetAttributes)(&attr, raw_pointer_cast(ptr))) {
		case HICUP_(Success): break;
		case HICUP_(ErrorInvalidDevice): assert(0); break;  // NOLINT(bugprone-branch-clone)
		case HICUP_(ErrorInvalidValue): assert(0); break;   // NOLINT(bugprone-branch-clone)
		default: assert(0);                                 // 71 enumeration values not handled in switch: 'hipErrorOutOfMemory', 'hipErrorNotInitialized', 'hipErrorDeinitialized'...
		}
		assert(attr.type == HICUP_(MemoryTypeManaged));
		return attr.device;
	}
};

namespace thrust {

// template<class T>
// auto operator+(::thrust::device_reference<T const> const& ref) {
// 	return static_cast<T>(ref);
// }

template<std::size_t N> struct priority : std::conditional_t<N == 0, std::true_type, priority<N-1>> {};

template<class FF>
class result_helper{
	template<class F> static constexpr auto _(priority<0>/**/, F const& fun) -> void;
	template<class F> static constexpr auto _(priority<1>/**/, F const& fun) -> decltype(fun());
	template<class F> static constexpr auto _(priority<2>/**/, F const& fun) -> decltype(fun(multi::index{}));
	template<class F> static constexpr auto _(priority<3>/**/, F const& fun) -> decltype(fun(multi::index{}, multi::index{}));
	template<class F> static constexpr auto _(priority<4>/**/, F const& fun) -> decltype(fun(multi::index{}, multi::index{}, multi::index{}));

 public:
	using type = decltype(_(priority<4>{}, std::declval<FF>()));
	result_helper() = default;
	void dummy() const {}
};

template<class F, class R = void>
struct device_function_t : F {
	using F::operator();

	using result_type = R;
};

template<class R = void, class F>
auto device_function(F&& other) {
	return device_function_t<F, R>{std::forward<F>(other)};
}

#define BOOST_MULTI_CAPTURE(...) [__VA_ARGS__]

#define BOOST_MULTI_DEVICE_LAMBDA_LEGACY(Cap, ...) ([=] () { \
			[[maybe_unused]] auto host_dummy = (Cap __host__ __VA_ARGS__); \
			using result_type = typename multi::thrust::result_helper<decltype(host_dummy)>::type; \
			return multi::thrust::device_function<result_type>(Cap __device__ __VA_ARGS__ ); \
		}()) \

#define BOOST_MULTI_DEVICE_LAMBDA(...) () { \
			[[maybe_unused]] auto host_dummy = ([=] __host__ __VA_ARGS__); \
			using result_type = typename multi::thrust::result_helper<decltype(host_dummy)>::type; \
			return multi::thrust::device_function<result_type>([=] __device__ __VA_ARGS__ ); \
		}() \

template<dimensionality_type D, class Proj>
class device_restriction_iterator {
	typename extents_t<D>::iterator it_;
	Proj                               proj_;

	device_restriction_iterator(typename extents_t<D>::iterator it, Proj proj) : it_{it}, proj_{proj} {}

	template<dimensionality_type, class>
	friend class device_restriction;

 public:
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wuninitialized"
#endif
	[[deprecated("to enable default construction")]]
	device_restriction_iterator() : proj_{proj_} {}  // hack to enable default construction
#ifdef __clang__
#pragma clang diagnostic pop
#endif
	// device_restriction_iterator() = default;  // cppcheck-suppress uninitMemberVar ; partially formed
	// constexpr iterator() {}  // = default;  // NOLINT(hicpp-use-equals-default,modernize-use-equals-default) TODO(correaa) investigate workaround

	device_restriction_iterator(device_restriction_iterator const& other) = default;
	device_restriction_iterator(device_restriction_iterator&&)            = default;

	auto operator=(device_restriction_iterator&&) -> device_restriction_iterator&      = default;
	auto operator=(device_restriction_iterator const&) -> device_restriction_iterator& = default;

	~device_restriction_iterator() = default;

	using system = typename multi::detail::function_system<Proj>::type;

	using difference_type = std::ptrdiff_t;
	using value_type      = int;
	//	std::conditional_t<(D != 1), restriction<D - 1, bind_front_t<Proj>>, decltype(apply_(std::declval<Proj>(), std::declval<typename extents_t<D>::element>()))>

	using pointer = void;

	using reference = int;  // std::conditional_t<(D != 1), restriction<D - 1, bind_front_t<Proj>>, decltype(apply_(std::declval<Proj&>(), std::declval<typename extents_t<D>::element>()))>;

	using iterator_category = std::random_access_iterator_tag;

	constexpr auto operator++() -> auto& {
		++it_;
		return *this;
	}
	constexpr auto operator--() -> auto& {
		--it_;
		return *this;
	}

	constexpr auto operator+=(difference_type dd) -> auto& {
		it_ += dd;
		return *this;
	}
	constexpr auto operator-=(difference_type dd) -> auto& {
		it_ -= dd;
		return *this;
	}

	constexpr auto operator++(int) {
		device_restriction_iterator ret{*this};
		++(*this);
		return ret;
	}
	constexpr auto operator--(int) {
		device_restriction_iterator ret{*this};
		--(*this);
		return ret;
	}

	friend constexpr auto operator-(device_restriction_iterator const& self, device_restriction_iterator const& other) { return self.it_ - other.it_; }
	template<class = void>
	friend __host__ __device__ constexpr auto operator+(device_restriction_iterator const& self, difference_type n) {
		device_restriction_iterator ret{self};
		return ret += n;
	}
	friend constexpr auto operator-(device_restriction_iterator const& self, difference_type n) {
		device_restriction_iterator ret{self};
		return ret -= n;
	}

	friend constexpr auto operator+(difference_type n, device_restriction_iterator const& self) { return self + n; }

	friend constexpr auto operator==(device_restriction_iterator const& self, device_restriction_iterator const& other) noexcept -> bool { return self.it_ == other.it_; }
	friend constexpr auto operator!=(device_restriction_iterator const& self, device_restriction_iterator const& other) noexcept -> bool { return self.it_ != other.it_; }

	friend auto operator<=(device_restriction_iterator const& self, device_restriction_iterator const& other) noexcept -> bool { return self.it_ <= other.it_; }
	friend auto operator<(device_restriction_iterator const& self, device_restriction_iterator const& other) noexcept -> bool { return self.it_ < other.it_; }
	friend auto operator>(device_restriction_iterator const& self, device_restriction_iterator const& other) noexcept -> bool { return self.it_ > other.it_; }
	friend auto operator>=(device_restriction_iterator const& self, device_restriction_iterator const& other) noexcept -> bool { return self.it_ >= other.it_; }

	__host__ __device__ constexpr auto operator*() -> int {
		// decltype(auto) {
		// if constexpr(D != 1) {
		// 	using std::get;
		// 	// auto ll = [idx = get<0>(*it_), proj = proj_](auto... rest) { return proj(idx, rest...); };
		// 	return device_restriction<D - 1, bind_front_t<Proj>>(extents_t<D - 1>((*it_).tail()), bind_front_t<Proj>{get<0>(*it_), *Pproj_});
		// } else {
		using std::get;
		return proj_(get<0>(*it_));
	}

	__host__ __device__ auto operator[](difference_type dd) const -> decltype(auto) { return *((*this) + dd); }  // TODO(correaa) use ra_iterator_facade
};

template<dimensionality_type D, class Proj>
class device_restriction {  //: restriction<D, Proj, int> {
	multi::extents_t<D> exts_;
	Proj                   proj_;

 public:
	device_restriction(multi::extents_t<D> exts, Proj proj) : exts_{exts}, proj_{proj} {}
	using iterator = device_restriction_iterator<D, Proj>;

	auto begin() const -> iterator {
		return iterator{exts_.begin(), proj_};
	}
	auto end() const -> iterator {
		return iterator{exts_.end(), proj_};
	}
};

#ifdef __cpp_deduction_guides
template<dimensionality_type D, typename Fun>
device_restriction(multi::extents_t<D>, Fun) -> device_restriction<D, Fun>;

template<typename Fun> device_restriction(extents_t<0>, Fun) -> device_restriction<0, Fun>;
template<typename Fun> device_restriction(extents_t<1>, Fun) -> device_restriction<1, Fun>;
template<typename Fun> device_restriction(extents_t<2>, Fun) -> device_restriction<2, Fun>;
template<typename Fun> device_restriction(extents_t<3>, Fun) -> device_restriction<3, Fun>;
#endif

template<dimensionality_type D, typename F>
auto device_restricted(F&& fun, extents_t<D> const& ext) {  // nvc++ has 'restrict' reserved
	return device_restriction<D, F>(ext, std::forward<F>(fun));
}


}  // namespace thrust

}  // end namespace boost::multi

// this is important for algorithms to dispatch to the right thrust executor
namespace thrust {

// template<class It> struct iterator_system;  // not needed in cuda 12.0, doesn't work on cuda 12.5

template<class T, ::boost::multi::dimensionality_type D, class Pointer, bool IsConst, bool IsMove, typename Stride>
struct iterator_system<::boost::multi::detail::array_iterator<T, D, Pointer, IsConst, IsMove, Stride>> {
	using type = typename ::thrust::iterator_system<typename ::boost::multi::detail::array_iterator<T, D, Pointer, IsConst, IsMove, Stride>::element_ptr>::type;
};

template<typename Pointer, class LayoutType>
struct iterator_system<::boost::multi::elements_iterator_t<Pointer, LayoutType>> {  // TODO(correaa) might need changes for IsConst templating
	using type = typename ::thrust::iterator_system<typename ::boost::multi::elements_iterator_t<Pointer, LayoutType>::pointer>::type;
};

template<class T, class UF, class Ptr, class Ref>
struct iterator_system<::boost::multi::transform_ptr<T, UF, Ptr, Ref>> {  // TODO(correaa) might need changes for IsConst templating
	using type = typename ::thrust::iterator_system<Ptr>::type;
};

template<::boost::multi::dimensionality_type D, typename Proj>
struct iterator_system<::boost::multi::thrust::device_restriction_iterator<D, Proj>> {  // TODO(correaa) might need changes for IsConst templating
	using type = typename ::thrust::iterator_system<::thrust::device_ptr<void>>::type;
};

// namespace detail {
// template<class T1, class T2, class LO>
// struct pointer_traits<
//  boost::multi::basic_array_ptr<
//      boost::multi::subarray<T1, 1L, thrust::pointer<T2, thrust::cuda_cub::tag, thrust::tagged_reference<T2, thrust::cuda_cub::tag>, thrust::use_default>, LO>,
//      LO
//  >
// >
// {
//  using Ptr = boost::multi::basic_array_ptr<
//      boost::multi::subarray<T1, 1L, thrust::pointer<T2, thrust::cuda_cub::tag, thrust::tagged_reference<T2, thrust::cuda_cub::tag>, thrust::use_default>, LO>,
//      LO
//  >;
//  using pointer = Ptr;
//  using reference = thrust::tagged_reference<T2, thrust::cuda_cub::tag>;
//   typedef typename pointer_element<Ptr>::type    element_type;
//   typedef typename pointer_difference<Ptr>::type difference_type;

//   template<typename U>
//     struct rebind
//   {
//     typedef typename rebind_pointer<Ptr,U>::type other;
//   };

// //  __host__ __device__
// //   inline static pointer pointer_to(typename pointer_traits_detail::pointer_to_param<element_type>::type r)
// //   {
// //     // XXX this is supposed to be pointer::pointer_to(&r); (i.e., call a static member function of pointer called pointer_to)
// //     //     assume that pointer has a constructor from raw pointer instead

// //     return pointer(&r);
// //   }

//   // thrust additions follow
//   //typedef typename pointer_raw_pointer<Ptr>::type raw_pointer;
//  using raw_pointer = boost::multi::basic_array_ptr<
//      boost::multi::subarray<T1, 1L, T2*, LO>,
//      LO
//  >;

//   __host__ __device__
//   inline static raw_pointer get(pointer ptr)
//   {
//  return reinterpret_cast<raw_pointer&>(ptr); //     return ptr.get();
//   }
// };
// }

}  // end namespace thrust

#if THRUST_VERSION >= 300200  // this is needed by CCCL 2
namespace thrust::detail {

template<typename T, ::boost::multi::dimensionality_type D, typename ElementPtr, class Layout, bool IsConst>
struct pointer_element<::boost::multi::detail::subarray_ptr<T, D, ElementPtr, Layout, IsConst>> {
	using type = std::conditional_t<D == 1,
		std::conditional_t<IsConst,
			std::add_const_t<typename thrust::detail::pointer_element<ElementPtr>::type>,
			typename thrust::detail::pointer_element<ElementPtr>::type
		>,
		std::conditional_t<IsConst,
			::boost::multi::const_subarray<T, D, ElementPtr, Layout>,
			::boost::multi::subarray<T, D, ElementPtr, Layout>
		>
		// void
	>;
};

}  // end namespace thrust::detail
#endif

namespace boost::multi::thrust {

// defines multi::thrust::device_array
// defines multi::thrust::host_array

template<typename T, multi::dimensionality_type D, class Alloc = ::thrust::device_allocator<T>> using device_array       = multi::array<T, D, Alloc>;
template<typename T, multi::dimensionality_type D, class Alloc = ::thrust::universal_allocator<T>> using universal_array = multi::array<T, D, Alloc>;
template<typename T, multi::dimensionality_type D, class Alloc = std::allocator<T>> using host_array                     = multi::array<T, D, Alloc>;

// defines multi::thrust::device::array
// defines multi::thrust::host  ::array
namespace device {
template<class T, multi::dimensionality_type D> using array = device_array<T, D>;
}  // end namespace device
namespace universal {
template<class T, multi::dimensionality_type D> using array = universal_array<T, D>;
}  // end namespace universal
namespace host {
template<class T, multi::dimensionality_type D> using array = host_array<T, D>;
}  // end namespace host

// defines multi::thrust::cuda::array
// defines multi::thrust::cuda::managed::array
namespace cuda {
template<class T, multi::dimensionality_type D> using array = multi::array<T, D, ::thrust::HICUP::allocator<T>>;

// namespace managed {
//  template<class T, multi::dimensionality_type D> using array = multi::array<T, D, boost::multi::thrust::hip::managed::allocator<T>>;
// }  // end namespace managed
}  // end namespace cuda

namespace mr {
template<class T, multi::dimensionality_type D, class MR> using array = array<T, D, ::thrust::mr::allocator<T, MR>>;
}
namespace pmr {
template<class T, multi::dimensionality_type D, class Pointer> using array = mr::array<T, D, ::thrust::mr::memory_resource<Pointer>>;
template<class T, multi::dimensionality_type D> using universal_array      = pmr::array<T, D, ::thrust::universal_ptr<void>>;
}  // end namespace pmr

namespace cuda {

template<class T, multi::dimensionality_type D> using universal_array = multi::array<T, D, ::thrust::HICUP::universal_allocator<T>>;

namespace universal {
template<class T, multi::dimensionality_type D> using array = multi::thrust::cuda::universal_array<T, D>;
}  // end namespace universal

namespace pmr {
template<class T, multi::dimensionality_type D> using universal_array = ::boost::multi::thrust::pmr::array<T, D, ::thrust::HICUP::universal_pointer<void>>;
}  // end namespace pmr
}  // end namespace cuda

}  // end namespace boost::multi::thrust

namespace boost::multi {

template<class Q, class R>
constexpr auto default_allocator_of(::thrust::pointer<Q, ::thrust::HICUP::tag, Q&> /*unused*/) {
	return ::thrust::HICUP::universal_allocator<typename std::iterator_traits<::thrust::pointer<Q, ::thrust::HICUP::tag, Q&>>::value_type>{};
}

}  // end namespace boost::multi

#undef HICUP
#undef HICUP_

#endif  // BOOST_MULTI_ADAPTORS_THRUST_HPP
