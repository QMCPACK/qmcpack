// Copyright 2018-2026 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_DETAIL_TYPES_HPP
#define BOOST_MULTI_DETAIL_TYPES_HPP
// #pragma once

#include <cstddef>      // for std::size_t
#include <type_traits>  // for make_signed_t

namespace boost::multi {

using ssize_t    = std::make_signed_t<std::size_t>;
using usize_t = std::make_unsigned_t<multi::ssize_t>;

using size_type [[deprecated("use boost::multi::ssize_t (for the default multi size type (signed) or std::size_t for the default STL size type (unsiged)")]] = std::make_signed_t<std::size_t>;
using size_t [[deprecated("use boost::multi::ssize_t (for the default multi size type (signed) or std::size_t for the default STL size type (unsiged)")]] = std::make_signed_t<std::size_t>;

using index           = std::make_signed_t<ssize_t>;
using difference_type = std::make_signed_t<index>;

using dimensionality_t    = index;
using dimensionality_type = dimensionality_t;

}  // end namespace boost::multi
#endif  // BOOST_MULTI_DETAIL_TYPES_HPP
