//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_DELAYED_UPDATE_CUDA_H
#define QMCPLUSPLUS_DELAYED_UPDATE_CUDA_H

#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "CUDA/CUDAruntime.hpp"
#include "CUDA/CUDAallocator.hpp"
#include "CUDA/cuBLAS.hpp"
#include "QMCWaveFunctions/detail/CUDA/delayed_update_helper.h"
#include "PrefetchedRange.h"
#if defined(QMC_CUDA2HIP)
#include "rocSolverInverter.hpp"
#else
#include "cuSolverInverter.hpp"
#endif

namespace qmcplusplus
{
/** implements delayed update on NVIDIA GPU using cuBLAS and cusolverDN
 * @tparam T base precision for most computation
 * @tparam T_FP high precision for matrix inversion, T_FP >= T
 */
template<typename T, typename T_FP>
class DelayedUpdateCUDA
{
  // Data staged during for delayed acceptRows
  Matrix<T, CUDAHostAllocator<T>> U;
  Matrix<T, CUDAHostAllocator<T>> Binv;
  Matrix<T> V;
  //Matrix<T> tempMat; // for debugging only
  Matrix<T, CUDAAllocator<T>> temp_gpu;
  /// GPU copy of U, V, Binv, Ainv
  Matrix<T, CUDAAllocator<T>> U_gpu;
  Matrix<T, CUDAAllocator<T>> V_gpu;
  Matrix<T, CUDAAllocator<T>> Binv_gpu;
  Matrix<T, CUDAAllocator<T>> Ainv_gpu;
  // auxiliary arrays for B
  Vector<T> p;
  Vector<int, CUDAHostAllocator<int>> delay_list;
  Vector<int, CUDAAllocator<int>> delay_list_gpu;
  /// current number of delays, increase one for each acceptance, reset to 0 after updating Ainv
  int delay_count;

#if defined(QMC_CUDA2HIP)
  rocSolverInverter<T_FP> rocsolver_inverter;
#else
  cuSolverInverter<T_FP> cusolver_inverter;
#endif

  // the range of prefetched_Ainv_rows
  PrefetchedRange prefetched_range;
  // Ainv prefetch buffer
  Matrix<T, CUDAHostAllocator<T>> Ainv_buffer;

  // CUDA specific variables
  cublasHandle_t h_cublas;
  cudaStream_t hstream;

  /// reset delay count to 0
  inline void clearDelayCount()
  {
    delay_count = 0;
    prefetched_range.clear();
  }

public:
  /// default constructor
  DelayedUpdateCUDA() : delay_count(0)
  {
    cudaErrorCheck(cudaStreamCreate(&hstream), "cudaStreamCreate failed!");
    cublasErrorCheck(cublasCreate(&h_cublas), "cublasCreate failed!");
    cublasErrorCheck(cublasSetStream(h_cublas, hstream), "cublasSetStream failed!");
  }

  ~DelayedUpdateCUDA()
  {
    cublasErrorCheck(cublasDestroy(h_cublas), "cublasDestroy failed!");
    cudaErrorCheck(cudaStreamDestroy(hstream), "cudaStreamDestroy failed!");
  }

  /** resize the internal storage
   * @param norb number of electrons/orbitals
   * @param delay, maximum delay 0<delay<=norb
   */
  inline void resize(int norb, int delay)
  {
    //tempMat.resize(norb, delay);
    V.resize(delay, norb);
    U.resize(delay, norb);
    p.resize(delay);
    Binv.resize(delay, delay);
    // prefetch 8% more rows corresponding to roughly 96% acceptance ratio
    Ainv_buffer.resize(std::min(static_cast<int>(delay * 1.08), norb), norb);

    temp_gpu.resize(norb, delay);
    delay_list.resize(delay);
    U_gpu.resize(delay, norb);
    V_gpu.resize(delay, norb);
    Binv_gpu.resize(delay, delay);
    delay_list_gpu.resize(delay);
    Ainv_gpu.resize(norb, norb);
  }

  /** compute the inverse of the transpose of matrix A and its determinant value in log
   * @tparam TREAL real type
   */
  template<typename TREAL>
  void invert_transpose(const Matrix<T>& logdetT, Matrix<T>& Ainv, std::complex<TREAL>& log_value)
  {
    clearDelayCount();
#if defined(QMC_CUDA2HIP)
    rocsolver_inverter.invert_transpose(logdetT, Ainv, Ainv_gpu, log_value);
#else
    cusolver_inverter.invert_transpose(logdetT, Ainv, Ainv_gpu, log_value);
#endif
  }

  /** initialize internal objects when Ainv is refreshed
   * @param Ainv inverse matrix
   */
  inline void initializeInv(const Matrix<T>& Ainv)
  {
    cudaErrorCheck(cudaMemcpyAsync(Ainv_gpu.data(), Ainv.data(), Ainv.size() * sizeof(T), cudaMemcpyHostToDevice,
                                   hstream),
                   "cudaMemcpyAsync failed!");
    clearDelayCount();
    // H2D transfer must be synchronized regardless of host memory being pinned or not.
    cudaErrorCheck(cudaStreamSynchronize(hstream), "cudaStreamSynchronize failed!");
  }

  inline int getDelayCount() const { return delay_count; }

  /** compute the row of up-to-date Ainv
   * @param Ainv inverse matrix
   * @param rowchanged the row id corresponding to the proposed electron
   */
  template<typename VVT>
  inline void getInvRow(const Matrix<T>& Ainv, int rowchanged, VVT& invRow)
  {
    if (!prefetched_range.checkRange(rowchanged))
    {
      int last_row = std::min(rowchanged + Ainv_buffer.rows(), Ainv.rows());
      cudaErrorCheck(cudaMemcpyAsync(Ainv_buffer.data(), Ainv_gpu[rowchanged],
                                     invRow.size() * (last_row - rowchanged) * sizeof(T), cudaMemcpyDeviceToHost,
                                     hstream),
                     "cudaMemcpyAsync failed!");
      prefetched_range.setRange(rowchanged, last_row);
      cudaErrorCheck(cudaStreamSynchronize(hstream), "cudaStreamSynchronize failed!");
    }
    // save AinvRow to new_AinvRow
    std::copy_n(Ainv_buffer[prefetched_range.getOffset(rowchanged)], invRow.size(), invRow.data());
    if (delay_count > 0)
    {
      constexpr T cone(1);
      constexpr T czero(0);
      const int norb     = Ainv.rows();
      const int lda_Binv = Binv.cols();
      // multiply V (NxK) Binv(KxK) U(KxN) AinvRow right to the left
      BLAS::gemv('T', norb, delay_count, cone, U.data(), norb, invRow.data(), 1, czero, p.data(), 1);
      BLAS::gemv('N', delay_count, delay_count, -cone, Binv.data(), lda_Binv, p.data(), 1, czero, Binv[delay_count], 1);
      BLAS::gemv('N', norb, delay_count, cone, V.data(), norb, Binv[delay_count], 1, cone, invRow.data(), 1);
    }
  }

  /** accept a move with the update delayed
   * @param Ainv inverse matrix
   * @param rowchanged the row id corresponding to the proposed electron
   * @param psiV new orbital values
   *
   * Before delay_count reaches the maximum delay, only Binv is updated with a recursive algorithm
   */
  template<typename VVT, typename RATIOT>
  inline void acceptRow(Matrix<T>& Ainv, int rowchanged, const VVT& psiV, const RATIOT ratio_new)
  {
    // update Binv from delay_count to delay_count+1
    constexpr T cone(1);
    constexpr T czero(0);
    const int norb     = Ainv.rows();
    const int lda_Binv = Binv.cols();
    std::copy_n(Ainv_buffer[prefetched_range.getOffset(rowchanged)], norb, V[delay_count]);
    std::copy_n(psiV.data(), norb, U[delay_count]);
    delay_list[delay_count] = rowchanged;
    // the new Binv is [[X Y] [Z sigma]]
    BLAS::gemv('T', norb, delay_count + 1, -cone, V.data(), norb, psiV.data(), 1, czero, p.data(), 1);
    // sigma
    const T sigma                  = static_cast<T>(RATIOT(1) / ratio_new);
    Binv[delay_count][delay_count] = sigma;
    // Y
    BLAS::gemv('T', delay_count, delay_count, sigma, Binv.data(), lda_Binv, p.data(), 1, czero,
               Binv.data() + delay_count, lda_Binv);
    // X
    BLAS::ger(delay_count, delay_count, cone, Binv[delay_count], 1, Binv.data() + delay_count, lda_Binv, Binv.data(),
              lda_Binv);
    // Z
    for (int i = 0; i < delay_count; i++)
      Binv[delay_count][i] *= sigma;
    delay_count++;
    // update Ainv when maximal delay is reached
    if (delay_count == lda_Binv)
      updateInvMat(Ainv, false);
  }

  /** update the full Ainv and reset delay_count
   * @param Ainv inverse matrix
   */
  inline void updateInvMat(Matrix<T>& Ainv, bool transfer_to_host = true)
  {
    // update the inverse matrix
    if (delay_count > 0)
    {
      constexpr T cone(1);
      constexpr T czero(0);
      constexpr T cminusone(-1);
      const int norb     = Ainv.rows();
      const int lda_Binv = Binv.cols();
      cudaErrorCheck(cudaMemcpyAsync(U_gpu.data(), U.data(), norb * delay_count * sizeof(T), cudaMemcpyHostToDevice,
                                     hstream),
                     "cudaMemcpyAsync failed!");
      cublasErrorCheck(cuBLAS::gemm(h_cublas, CUBLAS_OP_T, CUBLAS_OP_N, delay_count, norb, norb, &cone, U_gpu.data(),
                                    norb, Ainv_gpu.data(), norb, &czero, temp_gpu.data(), lda_Binv),
                       "cuBLAS::gemm failed!");
      cudaErrorCheck(cudaMemcpyAsync(delay_list_gpu.data(), delay_list.data(), delay_count * sizeof(int),
                                     cudaMemcpyHostToDevice, hstream),
                     "cudaMemcpyAsync failed!");
      applyW_stageV_cuda(delay_list_gpu.data(), delay_count, temp_gpu.data(), norb, temp_gpu.cols(), V_gpu.data(),
                         Ainv_gpu.data(), hstream);
      cudaErrorCheck(cudaMemcpyAsync(Binv_gpu.data(), Binv.data(), lda_Binv * delay_count * sizeof(T),
                                     cudaMemcpyHostToDevice, hstream),
                     "cudaMemcpyAsync failed!");
      cublasErrorCheck(cuBLAS::gemm(h_cublas, CUBLAS_OP_N, CUBLAS_OP_N, norb, delay_count, delay_count, &cone,
                                    V_gpu.data(), norb, Binv_gpu.data(), lda_Binv, &czero, U_gpu.data(), norb),
                       "cuBLAS::gemm failed!");
      cublasErrorCheck(cuBLAS::gemm(h_cublas, CUBLAS_OP_N, CUBLAS_OP_N, norb, norb, delay_count, &cminusone,
                                    U_gpu.data(), norb, temp_gpu.data(), lda_Binv, &cone, Ainv_gpu.data(), norb),
                       "cuBLAS::gemm failed!");
      clearDelayCount();
    }

    // transfer Ainv_gpu to Ainv and wait till completion
    if (transfer_to_host)
    {
      cudaErrorCheck(cudaMemcpyAsync(Ainv.data(), Ainv_gpu.data(), Ainv.size() * sizeof(T), cudaMemcpyDeviceToHost,
                                     hstream),
                     "cudaMemcpyAsync failed!");
      // no need to wait because : For transfers from device memory to pageable host memory, the function will return only once the copy has completed.
      //cudaErrorCheck(cudaStreamSynchronize(hstream), "cudaStreamSynchronize failed!");
    }
  }
};
} // namespace qmcplusplus

#endif // QMCPLUSPLUS_DELAYED_UPDATE_CUDA_H
