//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_SPARSE_CUDA_GPU_PTR_HPP
#define AFQMC_SPARSE_CUDA_GPU_PTR_HPP

#include <type_traits>
#include <cassert>
#include <vector>
#include "AFQMC/Memory/custom_pointers.hpp"
#include "AFQMC/Numerics/detail/CUDA/cusparse_wrapper.hpp"
#include "AFQMC/Numerics/detail/CUDA/cublas_wrapper.hpp"
#include <cassert>
#include <complex>

#include "multi/array.hpp"

namespace device
{
extern boost::multi::array<std::complex<double>, 1, device::device_allocator<std::complex<double>>>* cusparse_buffer;


template<typename T, typename Q>
void csrmv(const char transa,
           const int M,
           const int K,
           const T alpha,
           const char* matdescra,
           device_pointer<T> const A,
           device_pointer<int> const indx,
           device_pointer<int> const pntrb,
           device_pointer<int> const pntre,
           device_pointer<Q> const x,
           const T beta,
           device_pointer<T> y)
{
  static_assert(std::is_same<typename std::decay<Q>::type, T>::value, "Wrong dispatch.\n");
  // somehow need to check if the matrix is compact!
  int pb, pe;
  arch::memcopy(std::addressof(pb), to_address(pntrb), sizeof(int), arch::memcopyD2H, "sparse_cuda_gpu_ptr::csrmv");
  arch::memcopy(std::addressof(pe), to_address(pntre + (M - 1)), sizeof(int), arch::memcopyD2H, "sparse_cuda_gpu_ptr::csrmv");
  int nnz = pe - pb;
  if (CUSPARSE_STATUS_SUCCESS !=
      cusparse::cusparse_csrmv(*A.handles.cusparse_handle, transa, M, K, nnz, alpha,
                               qmc_cuda::afqmc_cusparse_matrix_descr, to_address(A), to_address(pntrb),
                               to_address(indx), to_address(x), beta, to_address(y)))
    throw std::runtime_error("Error: cusparse_csrmv returned error code.");
}

template<typename T, typename Q>
void csrmm(const char transa,
           const int M,
           const int N,
           const int K,
           const T alpha,
           const char* matdescra,
           device_pointer<T> A,
           device_pointer<int> indx,
           device_pointer<int> pntrb,
           device_pointer<int> pntre,
           device_pointer<Q> B,
           const int ldb,
           const T beta,
           device_pointer<T> C,
           const int ldc)
{
  static_assert(std::is_same<typename std::decay<Q>::type, T>::value, "Wrong dispatch.\n");
  // somehow need to check if the matrix is compact!
  int pb, pe;
  arch::memcopy(std::addressof(pb), to_address(pntrb), sizeof(int), arch::memcopyD2H, "lapack_sparse_gpu_ptr::csrmm");
  arch::memcopy(std::addressof(pe), to_address(pntre + (M - 1)), sizeof(int), arch::memcopyD2H, "lapack_sparse_gpu_ptr::csrmm");
  int nnz = pe - pb;
  if (transa == 'N')
  {
    /*
      // CSR_A * B = C  -->  (Fortran)  B^T * CSC_(CSR_A) = C^T
      if(CUSPARSE_STATUS_SUCCESS != cusparse::cusparse_gemmi(*A.handles.cusparse_handle,
            N,M,K,nnz,alpha,to_address(B),ldb,
            to_address(A),to_address(pntrb),to_address(indx),
            beta,to_address(C),ldc))
        throw std::runtime_error("Error: cusparse_csrmm(gemmi) returned error code.");
// */
    // replace this with call to gemmi!!!
    // /*
    char transb('T');
    // setup work space for column major matrix C
    if (cusparse_buffer == nullptr)
    {
      cusparse_buffer = new boost::multi::array<
          std::complex<double>, 1,
          device::device_allocator<std::complex<double>>>(typename boost::multi::layout_t<1u>::extensions_type{M * M},
                                                          device::device_allocator<std::complex<double>>{});
    }
    else if (cusparse_buffer->num_elements() < M * N)
      cusparse_buffer->reextent(typename boost::multi::layout_t<1u>::extensions_type{M * N});
    device_pointer<T> C_(cusparse_buffer->origin().pointer_cast<T>());

    // if beta != 0, transpose C into C_
    if (std::abs(beta) > 1e-12)
      if (CUBLAS_STATUS_SUCCESS !=
          cublas::cublas_geam(*A.handles.cublas_handle, transb, transa, M, N, T(1), to_address(C), ldc, T(0),
                              to_address(C_), M, to_address(C_), M))
        throw std::runtime_error("Error: cublas_geam returned error code.");

    // call csrmm2 on C_
    if (CUSPARSE_STATUS_SUCCESS !=
        cusparse::cusparse_csrmm2(*A.handles.cusparse_handle, transa, transb, M, N, K, nnz, alpha,
                                  qmc_cuda::afqmc_cusparse_matrix_descr, to_address(A), to_address(pntrb),
                                  to_address(indx), to_address(B), ldb, beta, to_address(C_), M))
      throw std::runtime_error("Error: cusparse_csrmm returned error code.");

    // transpose work matrix to row major on result C
    if (CUBLAS_STATUS_SUCCESS !=
        cublas::cublas_geam(*A.handles.cublas_handle, transb, transa, N, M, T(1), to_address(C_), M, T(0),
                            to_address(C), ldc, to_address(C), ldc))
      throw std::runtime_error("Error: cublas_geam returned error code.");
    // */
  }
  else
  {
    char transT('T');
    char transN('N');
    // setup work space for column major matrix B,C
    if (cusparse_buffer == nullptr)
    {
      cusparse_buffer = new boost::multi::array<
          std::complex<double>, 1,
          device::device_allocator<std::complex<double>>>(typename boost::multi::layout_t<1u>::extensions_type{(M + K) *
                                                                                                               N},
                                                          device::device_allocator<std::complex<double>>{});
    }
    else if (cusparse_buffer->num_elements() < (M + K) * N)
      cusparse_buffer->reextent(typename boost::multi::layout_t<1u>::extensions_type{(M + K) * N});
    // A is MxK
    // B should be MxN
    device_pointer<T> B_(cusparse_buffer->origin().pointer_cast<T>());
    // C should be KxN
    device_pointer<T> C_(B_ + M * N);

    // if beta != 0, transpose C into C_
    if (std::abs(beta) > 1e-12)
      if (CUBLAS_STATUS_SUCCESS !=
          cublas::cublas_geam(*A.handles.cublas_handle, transT, transN, K, N, T(1), to_address(C), ldc, T(0),
                              to_address(C_), K, to_address(C_), K))
        throw std::runtime_error("Error: cublas_geam returned error code. C");

    if (CUBLAS_STATUS_SUCCESS !=
        cublas::cublas_geam(*A.handles.cublas_handle, transT, transN, M, N, T(1), to_address(B), ldb, T(0),
                            to_address(B_), M, to_address(B_), M))
      throw std::runtime_error("Error: cublas_geam returned error code. B");

    // call csrmm2 on C_
    if (CUSPARSE_STATUS_SUCCESS !=
        cusparse::cusparse_csrmm(*A.handles.cusparse_handle, transa, M, N, K, nnz, alpha,
                                 qmc_cuda::afqmc_cusparse_matrix_descr, to_address(A), to_address(pntrb),
                                 to_address(indx), to_address(B_), M, beta, to_address(C_), K))
      throw std::runtime_error("Error: cusparse_csrmm returned error code.");

    // transpose work matrix to row major on result C
    if (CUBLAS_STATUS_SUCCESS !=
        cublas::cublas_geam(*A.handles.cublas_handle, transT, transN, N, K, T(1), to_address(C_), K, T(0),
                            to_address(C), ldc, to_address(C), ldc))
      throw std::runtime_error("Error: cublas_geam returned error code. C_");
  }
}

}; // namespace device


#endif
