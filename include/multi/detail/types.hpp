// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2021 Alfredo A. Correa

#ifndef MULTI_DETAIL_TYPES_HPP
#define MULTI_DETAIL_TYPES_HPP

#include "index_range.hpp"

#include<array>
#include<cassert>
#include<cstddef>
#include<tuple>        // for make_tuple
#include<type_traits>  // for make_signed_t
#include<utility>      // for forward

namespace boost {
namespace multi {

using size_t = std::make_signed_t<std::size_t>;
using size_type = std::make_signed_t<std::size_t>;

using index               = std::make_signed_t<size_type>;
using difference_type     = std::make_signed_t<index>;
using index_range         = range<index>;
using index_extension     = extension_t<index>;
using dimensionality_type = index;

using iextension = index_extension;
using irange     = index_range;

namespace detail {

template<typename, typename>
struct append_to_type_seq{};

template<typename T, typename... Ts, template<typename...> class TT>
struct append_to_type_seq<T, TT<Ts...> > {
    using type = TT<Ts..., T>;
};

template<typename T, dimensionality_type N, template<typename...> class TT = std::tuple>
struct repeat {
    using type = typename
        append_to_type_seq<
            T,
            typename repeat<T, N-1, TT>::type
        >::type;
};

template<typename T, template<typename...> class TT>
struct repeat<T, 0, TT> {
	using type = TT<>;
};

template<class T, std::size_t N>
constexpr auto array_size_impl(const std::array<T, N>&)
    -> std::integral_constant<std::size_t, N>;

template<class... T>
constexpr auto array_size_impl(const std::tuple<T...>&)
    -> std::integral_constant<std::size_t, std::tuple_size<std::tuple<T...>>{}>;

template<class Array>
using array_size = decltype(array_size_impl(std::declval<const Array&>()));

template<class Array>
constexpr auto static_size() -> std::decay_t<decltype(array_size<Array>::value)> {
	return array_size<Array>::value;
}
template<class Array>
constexpr auto static_size(Array const& /*unused*/) -> decltype(static_size<Array>()) {
	return static_size<Array>();
}

// TODO(correaa) consolidate with tuple_tail defined somewhere else
template<class Tuple>
constexpr auto head(Tuple&& t)
->decltype(std::get<0>(std::forward<Tuple>(t))) {
	return std::get<0>(std::forward<Tuple>(t)); }

template<typename Tuple, std::size_t... Ns>
constexpr auto tail_impl(std::index_sequence<Ns...> /*unused*/, Tuple&& t) {
	return std::make_tuple(std::get<Ns+1U>(std::forward<Tuple>(t))...);
}

template<class Tuple>
constexpr auto tail(Tuple const& t) {
	return tail_impl(std::make_index_sequence<(static_size<Tuple>())-1>(), t);
}

}  // end namespace detail

template<typename T, dimensionality_type D>
struct initializer_list{
	using type = std::initializer_list<typename initializer_list<T, D-1>::type>;
};
template<typename T>
struct initializer_list<T, 1>{using type = std::initializer_list<T>;};

template<typename T, dimensionality_type D>
using initializer_list_t = typename initializer_list<T, D>::type;

template<dimensionality_type D> using index_extensions = typename detail::repeat<index_extension, D>::type;

template<dimensionality_type D, class Tuple>
constexpr auto contains(index_extensions<D> const& ie, Tuple const& tp) {
//  using detail::head;
//  using detail::tail;
	return contains(head(ie), head(tp)) and contains(tail(ie), tail(tp));
}

}  // end namespace multi
}  // end namespace boost
#endif
