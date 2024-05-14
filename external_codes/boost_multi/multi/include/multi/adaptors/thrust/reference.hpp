// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2022-2023 Alfredo A. Correa

#pragma once

#include<thrust/detail/reference.h>  // hipthrust needs this

namespace thrust {

    // TODO(correaa) consider restrict this for universal memory only
    template<class T1, class Tag1, class T2, class Tag2>
    HD constexpr auto operator*(
        thrust::tagged_reference<T1, Tag1> const& r1, 
        thrust::tagged_reference<T2, Tag2> const& r2 
    )
    ->decltype(thrust::raw_reference_cast(r1) * thrust::raw_reference_cast(r2)) {
        return thrust::raw_reference_cast(r1) * thrust::raw_reference_cast(r2); }

    template<class T1, class Tag1, class T2>
    HD constexpr auto operator*(
        thrust::tagged_reference<T1, Tag1> const& r1, 
        T2 const& r2
    )
    ->decltype(thrust::raw_reference_cast(r1) * r2) {
        return thrust::raw_reference_cast(r1) * r2; }

    template<class T1, class T2, class Tag2>
    HD constexpr auto operator*(
        T1 const& r1,
        thrust::tagged_reference<T2, Tag2> const& r2
    )
    ->decltype(r2 * thrust::raw_reference_cast(r1)) {
        return r2 * thrust::raw_reference_cast(r1); }

}
