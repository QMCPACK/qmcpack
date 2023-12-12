// Copyright 2023 Alfredo A. Correa

#ifndef MULTI_DETAIL_IMPLICIT_CAST_HPP
#define MULTI_DETAIL_IMPLICIT_CAST_HPP
#pragma once

namespace boost::multi::detail {

template<class From, class To> constexpr bool is_implicitly_convertible_v = std::is_convertible_v<From, To>;
template<class From, class To> constexpr bool is_explicitly_convertible_v = std::is_constructible_v<To, From>;

template<class To, class From, std::enable_if_t<std::is_convertible<From, To>::value, int> =0>  // ::value (not _v) needed by intel's icpc 19
constexpr auto implicit_cast(From&& r) -> To {return static_cast<To>(r);}  // NOLINT(readability-identifier-length) std naming

template<class To, class From, std::enable_if_t<std::is_constructible<To, From>::value and not std::is_convertible<From, To>::value, int> =0>  // ::value (not _v) needed by intel's icpc 19
constexpr auto explicit_cast(From&& r) -> To {return static_cast<To>(r);}  // NOLINT(readability-identifier-length) std naming

}  // end namespace boost::multi::detail
#endif
