//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_CUBLAS_LU_HPP
#define QMCPLUSPLUS_CUBLAS_LU_HPP

#include <complex>
#include <type_traits>
#include <cuda_runtime_api.h>
#include <cuComplex.h>

namespace qmcplusplus
{
namespace cuBLAS_LU
{
template<class...>
struct disjunction : std::false_type
{};
template<class B1>
struct disjunction<B1> : B1
{};
template<class B1, class... Bn>
struct disjunction<B1, Bn...> : std::conditional_t<bool(B1::value), B1, disjunction<Bn...>>
{};

template<typename V1, typename V2, typename T>
struct OnTypesEqual : std::integral_constant<bool, std::is_same<V1, V2>::value>
{
  using type = T;
};

template<typename T>
struct default_type : std::true_type
{
  using type = T;
};

template<typename T>
using TypesMapper = typename disjunction<OnTypesEqual<T, std::complex<double>, cuDoubleComplex>,
                                         OnTypesEqual<T, std::complex<float>, cuComplex>,
                                         OnTypesEqual<T, std::complex<double>*, cuDoubleComplex*>,
                                         OnTypesEqual<T, std::complex<float>*, cuComplex*>,
                                         OnTypesEqual<T, const std::complex<double>*, const cuDoubleComplex*>,
                                         OnTypesEqual<T, const std::complex<float>*, const cuComplex*>,
                                         default_type<void>>::type;

// There must be a way to make this a template function. But it isn't easy.
#define CUDATYPECAST(var) reinterpret_cast<TypesMapper<decltype(var)>>(var)

void computeInverseAndDetLog_batched(cublasHandle_t& h_cublas,
                                     cudaStream_t& hstream,
                                     const int n,
                                     const int lda,
                                     double** Ms,
                                     double** Cs,
                                     double* LU_diags,
                                     int* pivots,
                                     int* infos,
                                     std::complex<double>* log_dets,
                                     const int batch_size);

template<typename T>
void computeGetrf_batched(cublasHandle_t& h_cublas,
                          const int n,
                          const int lda,
                          T* Ms[],
                          int* pivots,
                          int* infos,
                          const int batch_size);

template<typename T>
void computeLogDet_batched(cudaStream_t& hstream,
                           const int n,
                           const T* LU_diags,
                           const int* pivots,
                           std::complex<double>* logdets,
                           const int batch_size);

void computeLUDiag_batched(cudaStream_t& hstream,
                           const int n,
                           const int lda,
                           double* Ms[],
                           double* LU_diags,
                           const int batch_size);


void computeGetri_batched(cublasHandle_t& h_cublas,
                          const int n,
                          const int lda,
                          double* Ms[],
                          double* Cs[],
                          int* pivots,
                          int* infos,
                          const int batch_size);

void peekinvM_batched(cudaStream_t& hstream,
                      double* M[],
                      double* invM[],
                      int* pivots,
                      int* infos,
                      const int batch_size);


} // namespace cuBLAS_LU
} // namespace qmcplusplus
#endif
