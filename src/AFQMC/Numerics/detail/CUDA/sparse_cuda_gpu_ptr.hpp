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
#include <complex>
#include <cuda_runtime.h>
#include "cusparse.h"
#include "AFQMC/Memory/CUDA/cuda_utilities.h"
#include "AFQMC/Memory/custom_pointers.hpp"
#include "AFQMC/Memory/buffer_managers.h"
#include "AFQMC/Numerics/detail/CUDA/cublas_wrapper.hpp"
#if CUSPARSE_VER_MAJOR < 11
#include "AFQMC/Numerics/detail/CUDA/cusparse_wrapper_deprecated.hpp"
#endif
#if CUSPARSE_VERSION >= 11400
#define CUSPARSE_CSRMV_ALG1 CUSPARSE_SPMV_CSR_ALG1
#endif

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
  using qmc_cuda::cusparse_check;
  static_assert(std::is_same<typename std::decay<Q>::type, T>::value, "Wrong dispatch.\n");
  // somehow need to check if the matrix is compact!
  int pb, pe;
  arch::memcopy(std::addressof(pb), to_address(pntrb), sizeof(int), arch::memcopyD2H, "sparse_cuda_gpu_ptr::csrmv");
  arch::memcopy(std::addressof(pe), to_address(pntre + (M - 1)), sizeof(int), arch::memcopyD2H,
                "sparse_cuda_gpu_ptr::csrmv");
  int nnz = pe - pb;
#if CUSPARSE_VER_MAJOR < 11
  cusparse_check(cusparse::cusparse_csrmv(*A.handles.cusparse_handle, transa, M, K, nnz, alpha,
                                          qmc_cuda::afqmc_cusparse_matrix_descr, to_address(A), to_address(pntrb),
                                          to_address(indx), to_address(x), beta, to_address(y)),
                 "Error: cusparse_csrmv returned error code.");
#else
  cusparseSpMatDescr_t matA;
  cusparseDnVecDescr_t vecx, vecy;
  size_t bufferSize = 0;

  cusparse_check(cusparseCreateCsr(&matA, M, K, nnz, (void*)to_address(pntrb), (void*)to_address(indx),
                                   (void*)to_address(A), CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                   CUSPARSE_INDEX_BASE_ZERO, qmc_cuda::cusparse_data_type<T>()),
                 "csrmv::cusparseCreateCsr");
  cusparse_check(cusparseCreateDnVec(&vecx, (transa == 'N' ? K : M), (void*)to_address(x),
                                     qmc_cuda::cusparse_data_type<T>()),
                 "csrmv::cusparseCreateDnMat");
  cusparse_check(cusparseCreateDnVec(&vecy, (transa == 'N' ? M : K), (void*)to_address(y),
                                     qmc_cuda::cusparse_data_type<T>()),
                 "csrmv::cusparseCreateDnMat");
  cusparse_check(cusparseSpMV_bufferSize(*A.handles.cusparse_handle, qmc_cuda::cusparseOperation(transa), &alpha, matA,
                                         vecx, &beta, vecy, qmc_cuda::cusparse_data_type<T>(), CUSPARSE_CSRMV_ALG1,
                                         &bufferSize),
                 "csrmv::cusparseSpMV_bufferSize");

  if (bufferSize > 0)
  {
    using qmcplusplus::afqmc::DeviceBufferManager;
    using pointer_t = typename DeviceBufferManager::template allocator_t<T>::pointer;
    DeviceBufferManager buffer_manager;
    auto alloc(buffer_manager.get_generator().template get_allocator<T>());
    auto ptr = alloc.allocate(bufferSize);
    cusparse_check(cusparseSpMV(*A.handles.cusparse_handle, qmc_cuda::cusparseOperation(transa), &alpha, matA, vecx,
                                &beta, vecy, qmc_cuda::cusparse_data_type<T>(), CUSPARSE_CSRMV_ALG1,
                                (void*)to_address(ptr)),
                   "csrmv::cusparseSpMV");
    alloc.deallocate(ptr, bufferSize);
  }
  else
  {
    void* dBuffer = NULL;
    cusparse_check(cusparseSpMV(*A.handles.cusparse_handle, qmc_cuda::cusparseOperation(transa), &alpha, matA, vecx,
                                &beta, vecy, qmc_cuda::cusparse_data_type<T>(), CUSPARSE_CSRMV_ALG1, dBuffer),
                   "csrmv::cusparseSpMV");
  }

  cusparse_check(cusparseDestroySpMat(matA), "csrmv::destroyA");
  cusparse_check(cusparseDestroyDnVec(vecx), "csrmv::destroyX");
  cusparse_check(cusparseDestroyDnVec(vecy), "csrmv::destroyY");
#endif
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
  using qmc_cuda::cusparse_check;
  static_assert(std::is_same<typename std::decay<Q>::type, T>::value, "Wrong dispatch.\n");
  // somehow need to check if the matrix is compact!
  int pb, pe;
  arch::memcopy(std::addressof(pb), to_address(pntrb), sizeof(int), arch::memcopyD2H, "lapack_sparse_gpu_ptr::csrmm");
  arch::memcopy(std::addressof(pe), to_address(pntre + (M - 1)), sizeof(int), arch::memcopyD2H,
                "lapack_sparse_gpu_ptr::csrmm");
  int nnz = pe - pb;
#if CUSPARSE_VER_MAJOR < 11
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
#else
  cusparseSpMatDescr_t matA;
  cusparseDnMatDescr_t matB, matC;
  size_t bufferSize = 0;
  size_t M_         = ((transa == 'N') ? M : K);
  size_t K_         = ((transa == 'N') ? K : M);

  cusparse_check(cusparseCreateCsr(&matA, M, K, nnz, (void*)to_address(pntrb), (void*)to_address(indx),
                                   (void*)to_address(A), CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                   CUSPARSE_INDEX_BASE_ZERO, qmc_cuda::cusparse_data_type<T>()),
                 "csrmm::cusparseCreateCsr");
  cusparse_check(cusparseCreateDnMat(&matB, K_, N, ldb, (void*)to_address(B), qmc_cuda::cusparse_data_type<T>(),
                                     CUSPARSE_ORDER_ROW),
                 "csrmm::cusparseCreateDnMat");
  cusparse_check(cusparseCreateDnMat(&matC, M_, N, ldc, (void*)to_address(C), qmc_cuda::cusparse_data_type<T>(),
                                     CUSPARSE_ORDER_ROW),
                 "csrmm::cusparseCreateDnMat");
  cusparse_check(cusparseSpMM_bufferSize(*A.handles.cusparse_handle, qmc_cuda::cusparseOperation(transa),
                                         CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, matB, &beta, matC,
                                         qmc_cuda::cusparse_data_type<T>(), CUSPARSE_SPMM_CSR_ALG2, &bufferSize),
                 "csrmm::cusparseSpMM_bufferSize");

  if (bufferSize > 0)
  {
    using qmcplusplus::afqmc::DeviceBufferManager;
    using pointer_t = typename DeviceBufferManager::template allocator_t<T>::pointer;
    DeviceBufferManager buffer_manager;
    auto alloc(buffer_manager.get_generator().template get_allocator<T>());
    auto ptr = alloc.allocate(bufferSize);
    cusparse_check(cusparseSpMM(*A.handles.cusparse_handle, qmc_cuda::cusparseOperation(transa),
                                CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, matB, &beta, matC,
                                qmc_cuda::cusparse_data_type<T>(), CUSPARSE_SPMM_CSR_ALG2, (void*)to_address(ptr)),
                   "csrmm::cusparseSpMM");
    alloc.deallocate(ptr, bufferSize);
  }
  else
  {
    void* dBuffer = NULL;
    cusparse_check(cusparseSpMM(*A.handles.cusparse_handle, qmc_cuda::cusparseOperation(transa),
                                CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, matA, matB, &beta, matC,
                                qmc_cuda::cusparse_data_type<T>(), CUSPARSE_SPMM_CSR_ALG2, dBuffer),
                   "csrmm::cusparseSpMM");
  }

  cusparse_check(cusparseDestroySpMat(matA), "csrmm::destroyA");
  cusparse_check(cusparseDestroyDnMat(matB), "csrmm::destroyB");
  cusparse_check(cusparseDestroyDnMat(matC), "csrmm::destroyC");
#endif
}

}; // namespace device


#endif
