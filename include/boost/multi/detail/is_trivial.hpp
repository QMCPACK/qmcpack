// Copyright 2022-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_DETAIL_IS_TRIVIAL_HPP
#define BOOST_MULTI_DETAIL_IS_TRIVIAL_HPP

#include <type_traits>

namespace boost {  // NOLINT(modernize-concat-nested-namespaces)
namespace multi {

template<class T> struct is_trivially_default_constructible : std::is_trivially_default_constructible<T> {};
template<class T> struct is_trivial : std::is_trivial<T> {};


}  // end namespace multi
}  // end namespace boost

#endif  // BOOST_MULTI_DETAIL_IS_TRIVIAL_HPP
