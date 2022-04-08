// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2022 Alfredo A. Correa

#pragma once

#include<type_traits>

namespace boost {  // NOLINT(modernize-concat-nested-namespaces)
namespace multi {

template<class T> struct is_trivially_default_constructible : std::is_trivially_default_constructible<T> {};

}  // end namespace multi
}  // end namespace boost
