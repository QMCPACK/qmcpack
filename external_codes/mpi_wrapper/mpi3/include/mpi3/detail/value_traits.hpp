// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2017-2023 Alfredo A. Correa

#ifndef BOOST_MPI3_DETAIL_VALUE_TRAITS_HPP
#define BOOST_MPI3_DETAIL_VALUE_TRAITS_HPP

//#include "../../mpi3/detail/package_archive.hpp"
//#include "../../mpi3/communicator.hpp"

#include <mpi3/detail/iterator.hpp>

#include <mpi3/type.hpp>

#include<iterator>
#include<type_traits>

namespace boost {
namespace mpi3 {
namespace detail {

/*
template<class T> 
struct is_memcopyable : std::integral_constant<bool, std::is_trivially_copyable<T>{}>{};

template<class T, class... Ts> 
struct is_memcopyable<std::tuple<T, Ts...>> : 
    std::integral_constant<bool, 
        is_memcopyable<T>{} and is_memcopyable<std::tuple<Ts...>>{}
    >
{};

template<class T1, class T2> 
struct is_memcopyable<std::pair<T1, T2>> : 
    std::integral_constant<bool, 
        is_memcopyable<T1>{} and is_memcopyable<T2>{}
    >
{};
*/
/*
template<
	class T, 
	typename = decltype(std::declval<mpi3::detail::package_oarchive&>() << std::declval<T>()),
	typename = decltype(std::declval<mpi3::detail::package_iarchive&>() >> std::declval<T>())
>
std::true_type  is_serializable_aux(T);
std::false_type is_serializable_aux(...);
*/

template<class T>
struct is_serializable : decltype(is_serializable_aux(std::declval<T>())){};

struct value_unspecified_tag{};
struct nonmemcopyable_serializable_tag : value_unspecified_tag{
	using base = value_unspecified_tag;
};
struct memcopyable_tag : value_unspecified_tag{
	using base = value_unspecified_tag;
};
struct basic_tag : memcopyable_tag{using base = memcopyable_tag;};

//template<class V, typename = std::enable_if_t<is_memcopyable<V>{} and not is_basic<V>{}>>
//memcopyable_tag value_category_aux(V&&);
template<class V, std::enable_if_t<is_basic<V>{} or has_datatype<V>{}, int> =0>
basic_tag value_category_aux(V const&);

//template<class V, std::enable_if_t<has_datatype<V>{}, int> =0>
//basic_tag value_category_aux(V const&);

template<class V, std::size_t N>
value_unspecified_tag value_category_aux(V[N]);  // NOLINT(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays) : legacy compatibility
value_unspecified_tag value_category_aux(... );

template<class V>
struct value_category{using type = decltype(value_category_aux(std::declval<V>()));};  // NOLINT(cppcoreguidelines-pro-type-vararg,hicpp-vararg)

template<class V>
using value_category_t = typename value_category<V>::type;

}  // end namespace detail
}  // end namespace mpi3
}  // end namespace boost

#endif
