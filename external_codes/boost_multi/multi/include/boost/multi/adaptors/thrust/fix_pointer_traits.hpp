// Copyright 2021-2026 Alfredo A. Correa

#ifndef BOOST_MULTI_ADAPTORS_THRUST_FIX_POINTER_TRAITS_HPP
#define BOOST_MULTI_ADAPTORS_THRUST_FIX_POINTER_TRAITS_HPP
#pragma once

#include <thrust/detail/type_traits/pointer_traits.h>
#include <thrust/memory.h>

#include <memory>  // for std::pointer_traits

// #if(__CUDACC_VER_MAJOR__ * 10000 + __CUDACC_VER_MINOR__ * 100 + __CUDACC_VER_BUILD__ < 120500)
// begin of nvcc thrust 11.5 workaround : https://github.com/NVIDIA/thrust/issues/1629
namespace thrust {

// template<typename Element, typename Tag, typename Reference, typename Derived> class pointer;
// template<class T> struct pointer_traits;

}  // end namespace thrust

namespace boost::multi::detail {
// Primary: for references, form rebound reference via a helper
template<class Ref, class U> struct rebind_reference {
	using type = Ref;
};

// Helper: forming U& safely (void short-circuits to void)
template<class U> struct make_lvalue_ref {
	using type = U&;
};
template<> struct make_lvalue_ref<void> {
	using type = void;
};
template<> struct make_lvalue_ref<void const> {
	using type = void const;
};
template<> struct make_lvalue_ref<void volatile> {
	using type = void volatile;
};
template<> struct make_lvalue_ref<void const volatile> {
	using type = void const volatile;
};

template<class T, class U>
struct rebind_reference<T&, U> {
	using type = typename make_lvalue_ref<U>::type;
};

template<class Tagged, class U> struct make_tagged_ref;

template<class T, class Tag, class U>
struct make_tagged_ref<::thrust::tagged_reference<T, Tag>, U> {
	using type = ::thrust::tagged_reference<U, Tag>;
};

template<class T, class Tag>
struct make_tagged_ref<::thrust::tagged_reference<T, Tag>, void> {
	using type = void;
};

template<class T, class Tag>
struct make_tagged_ref<::thrust::tagged_reference<T, Tag>, void const> {
	using type =
		void const;
};
template<class T, class Tag>
struct make_tagged_ref<::thrust::tagged_reference<T, Tag>, void volatile> {
	using type = void volatile;
};
template<class T, class Tag>
struct make_tagged_ref<::thrust::tagged_reference<T, Tag>, void const volatile> {
	using type = void const volatile;
};

template<class T, class Tag, class U>
struct rebind_reference<::thrust::tagged_reference<T, Tag>, U> {
	using type = typename make_tagged_ref<::thrust::tagged_reference<T, Tag>, U>::type;
};
}  // namespace boost::multi::detail

template<class T, class Tag, class Ref, class Derived>
struct std::pointer_traits<::thrust::pointer<T, Tag, Ref, Derived>> {  // NOLINT(cert-dcl58-cpp) to specialize pointer traits
	using element_type    = T;
	using pointer         = ::thrust::pointer<T, Tag, Ref, Derived>;
	using difference_type = typename ::thrust::pointer<T, Tag, Ref, Derived>::difference_type;

	template<class U>
	using rebind = ::thrust::pointer<
		U, Tag,
		typename boost::multi::detail::rebind_reference<Ref, U>::type,
		Derived>;
};

// template<class... As>
// struct std::pointer_traits<::thrust::pointer<As...>>  // NOLINT(cert-dcl58-cpp) normal way to specialize pointer_traits
// : ::thrust::detail::pointer_traits<thrust::pointer<As...>> {
// 	template<class T>
// 	using rebind = typename ::thrust::detail::pointer_traits<::thrust::pointer<As...>>::template rebind<T>::other;
// };

// template<class T, class Tag, class Derived>
// struct std::pointer_traits<
// 	::thrust::pointer<T, Tag, ::thrust::tagged_reference<T, Tag>, Derived>
// > : ::thrust::detail::pointer_traits<
// 		::thrust::pointer<T, Tag, ::thrust::tagged_reference<T, Tag>, Derived>
// 	>
// {
// 	template<class U>
// 	using rebind = ::thrust::pointer<
// 		U, Tag,
// 		::thrust::tagged_reference<U, Tag>,  // keep Reference in sync with Element
// 		Derived
// 	>;
// };

// end of nvcc thrust 11.5 workaround
// #endif

#endif
