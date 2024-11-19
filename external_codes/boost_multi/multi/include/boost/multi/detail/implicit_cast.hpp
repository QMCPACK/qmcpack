// Copyright 2023-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_DETAIL_IMPLICIT_CAST_HPP
#define BOOST_MULTI_DETAIL_IMPLICIT_CAST_HPP
#pragma once

#include <type_traits>
#include <utility>

namespace boost::multi::detail {

template<class From, class To> constexpr bool is_implicitly_convertible_v = std::is_convertible_v<From, To>;
template<class From, class To> constexpr bool is_explicitly_convertible_v = std::is_constructible_v<To, From>;

template<class To, class From, std::enable_if_t<std::is_convertible_v<From, To>, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa)
constexpr auto implicit_cast(From&& ref) -> To {return static_cast<To>(std::forward<From&&>(ref));}

template<class To, class From, std::enable_if_t<std::is_constructible_v<To, From> &&  ! std::is_convertible_v<From, To>, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa)
constexpr auto explicit_cast(From&& ref) -> To {return static_cast<To>(std::forward<From&&>(ref));}

}  // end namespace boost::multi::detail
#endif
