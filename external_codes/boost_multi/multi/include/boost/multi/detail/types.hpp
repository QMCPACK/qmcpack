// Copyright 2018-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_DETAIL_TYPES_HPP
#define BOOST_MULTI_DETAIL_TYPES_HPP

#include<cstddef>      // for std::size_t
#include<type_traits>  // for make_signed_t

namespace boost::multi {

using size_t    = std::make_signed_t<std::size_t>;
using size_type = std::make_signed_t<std::size_t>;

using index               = std::make_signed_t<size_type>;
using difference_type     = std::make_signed_t<index>;

using dimensionality_type = index;

}  // end namespace boost::multi
#endif  // BOOST_MULTI_DETAIL_TYPES_HPP
