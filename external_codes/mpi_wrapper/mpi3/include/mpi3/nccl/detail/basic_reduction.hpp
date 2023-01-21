// -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
// Copyright 2022 Alfredo A. Correa

#ifndef MPI3_NCCL_DETAIL_BASIC_REDUCTION_HPP_
#define MPI3_NCCL_DETAIL_BASIC_REDUCTION_HPP_

#include <nccl.h>

namespace boost {
namespace mpi3 {
namespace nccl {

template<class T = void> using  plus = std::plus<T>;
template<class T = void> using  multiplies = std::multiplies<>;
template<class T = void> struct min {};
template<class T = void> struct max {};
template<class T = void> struct average {};

namespace detail {

template<class T> struct basic_reduction_t;

template<class Op> struct basic_reduction_t;

template<> struct basic_reduction_t<std::plus<>      > : std::integral_constant<ncclRedOp_t, ncclSum > {};
template<> struct basic_reduction_t<std::multiplies<>> : std::integral_constant<ncclRedOp_t, ncclProd> {};
template<> struct basic_reduction_t<nccl::min<>      > : std::integral_constant<ncclRedOp_t, ncclMin > {};
template<> struct basic_reduction_t<nccl::max<>      > : std::integral_constant<ncclRedOp_t, ncclMax > {};
template<> struct basic_reduction_t<nccl::average<>  > : std::integral_constant<ncclRedOp_t, ncclAvg > {};

template<class T> auto basic_reduction = basic_reduction_t<T>::value;

}

}}}

#endif
