// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2018-2021 Alfredo A. Correa

#ifndef MULTI_DETAIL_TYPES_HPP
#define MULTI_DETAIL_TYPES_HPP

// #include "index_range.hpp"

#include<cstddef>      // for std::size_t
// #include<tuple>        // for make_tuple
#include<type_traits>  // for make_signed_t
#include<utility>      // for forward

namespace boost::multi {

using size_t    = std::make_signed_t<std::size_t>;
using size_type = std::make_signed_t<std::size_t>;

using index               = std::make_signed_t<size_type>;
using difference_type     = std::make_signed_t<index>;
using dimensionality_type = index;

}  // end namespace boost::multi
#endif
