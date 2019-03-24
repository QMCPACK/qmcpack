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

#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/OhmmsMatrix.h>
#include "CUDA/CUDAallocator.hpp"
#include <Numerics/CUDA/cuBLAS.hpp>
#include "QMCWaveFunctions/Fermion/delayed_update_helper.h"
#include <cuda_runtime_api.h>
#include "CUDA/cudaError.h"

namespace qmcplusplus
{

/// helper class for the prefetched range of a vector
class Range
{
  // [first, last) rows of Ainv
  int first, last;

public:
  Range() : first(0), last(0){};
  void setRange(int first_in, int last_in)
  {
    first = first_in;
    last  = last_in;
  }
  inline void clear() { first = last = 0; };
  inline int getOffset(int index) const
  {
    if (!checkRange(index))
      throw std::runtime_error("index not in range \n");
    return index - first;
  }
  inline bool checkRange(int index) const { return (index >= first) && (index < last); };
};

template<typename T, typename T_FP>
class DelayedUpdateCUDA
{
  // Data staged during for delayed acceptRows
  Matrix<T, CUDAHostAllocator<T>> U, Binv;
  Matrix<T> V;
  //Matrix<T> tempMat; // for debugging only
  Matrix<T, CUDAAllocator<T>, MemorySpace::CUDA> U_gpu, V_gpu, Binv_gpu, temp_gpu, Ainv_gpu;
  // auxiliary arrays for B
  Vector<T> p;
  Vector<int, CUDAHostAllocator<int>> delay_list;
  Vector<int, CUDAAllocator<int>, MemorySpace::CUDA> delay_list_gpu;
  int delay_count;

  // the range of prefetched_Ainv_rows
  Range prefetched_range;
  // Ainv prefetch buffer
  Matrix<T, CUDAHostAllocator<T>> Ainv_buffer;

  // CUDA specific variables
  cublasHandle_t handle;
  cudaStream_t hstream;

  inline void waitStream() { cudaErrorCheck(cudaStreamSynchronize(hstream), "cudaStreamSynchronize failed!"); }

public:
  /// default constructor
  DelayedUpdateCUDA() : delay_count(0)
  {
    cudaErrorCheck(cudaStreamCreate(&hstream), "cudaStreamCreate failed!");
    cublasErrorCheck(cublasCreate(&handle), "cublasCreate failed!");
    cublasErrorCheck(cublasSetStream(handle, hstream), "cublasSetStream failed!");
  }

  ~DelayedUpdateCUDA()
  {
    cublasErrorCheck(cublasDestroy(handle), "cublasDestroy failed!");
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
    Ainv_gpu.resize(norb, norb);
    delay_list_gpu.resize(delay);
  }

  /** initialize internal objects when Ainv is refreshed
   * @param Ainv inverse matrix
   */
  inline void initializeInv(const Matrix<T>& Ainv)
  {
    cudaErrorCheck(cudaMemcpyAsync(Ainv_gpu.data(), Ainv.data(), Ainv.size() * sizeof(T), cudaMemcpyHostToDevice,
                                   hstream),
                   "cudaMemcpyAsync failed!");
    // safe mechanism
    delay_count = 0;
    prefetched_range.clear();
  }

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
      waitStream();
    }
    // save AinvRow to new_AinvRow
    std::copy_n(Ainv_buffer[prefetched_range.getOffset(rowchanged)], invRow.size(), invRow.data());
    if (delay_count > 0)
    {
      const T cone(1);
      const T czero(0);
      const int norb     = Ainv.rows();
      const int lda_Binv = Binv.cols();
      // multiply V (NxK) Binv(KxK) U(KxN) AinvRow right to the left
      BLAS::gemv('T', norb, delay_count, cone, U.data(), norb, invRow.data(), 1, czero, p.data(), 1);
      BLAS::gemv('N', delay_count, delay_count, cone, Binv.data(), lda_Binv, p.data(), 1, czero, Binv[delay_count], 1);
      BLAS::gemv('N', norb, delay_count, -cone, V.data(), norb, Binv[delay_count], 1, cone, invRow.data(), 1);
    }
  }

  /** accept a move with the update delayed
   * @param Ainv inverse matrix
   * @param rowchanged the row id corresponding to the proposed electron
   * @param psiV new orbital values
   *
   * Before delay_count reaches the maximum delay, only Binv is updated with a recursive algorithm
   */
  template<typename VVT>
  inline void acceptRow(Matrix<T>& Ainv, int rowchanged, const VVT& psiV)
  {
    // update Binv from delay_count to delay_count+1
    const T cminusone(-1);
    const T czero(0);
    const int norb     = Ainv.rows();
    const int lda_Binv = Binv.cols();
    std::copy_n(Ainv_buffer[prefetched_range.getOffset(rowchanged)], norb, V[delay_count]);
    std::copy_n(psiV.data(), norb, U[delay_count]);
    delay_list[delay_count] = rowchanged;
    // the new Binv is [[X Y] [Z x]]
    BLAS::gemv('T', norb, delay_count + 1, cminusone, V.data(), norb, psiV.data(), 1, czero, p.data(), 1);
    // x
    T y = -p[delay_count];
    for (int i = 0; i < delay_count; i++)
      y += Binv[delay_count][i] * p[i];
    Binv[delay_count][delay_count] = y = T(1) / y;
    // Y
    BLAS::gemv('T', delay_count, delay_count, y, Binv.data(), lda_Binv, p.data(), 1, czero, Binv.data() + delay_count,
               lda_Binv);
    // X
    BLAS::ger(delay_count, delay_count, cminusone, Binv[delay_count], 1, Binv.data() + delay_count, lda_Binv,
              Binv.data(), lda_Binv);
    // Z
    for (int i = 0; i < delay_count; i++)
      Binv[delay_count][i] *= -y;
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
      const T cone(1);
      const T czero(0);
      const int norb     = Ainv.rows();
      const int lda_Binv = Binv.cols();
      const T cminusone(-1);
      cudaErrorCheck(cudaMemcpyAsync(U_gpu.data(), U.data(), norb * delay_count * sizeof(T), cudaMemcpyHostToDevice,
                                     hstream),
                     "cudaMemcpyAsync failed!");
      cuBLAS::gemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, delay_count, norb, norb, &cone, U_gpu.data(), norb,
                   Ainv_gpu.data(), norb, &czero, temp_gpu.data(), lda_Binv);
      cudaErrorCheck(cudaMemcpyAsync(delay_list_gpu.data(), delay_list.data(), delay_count * sizeof(int),
                                     cudaMemcpyHostToDevice, hstream),
                     "cudaMemcpyAsync failed!");
      applyW_stageV_cuda(delay_list_gpu.data(), delay_count, temp_gpu.data(), norb, temp_gpu.cols(), V_gpu.data(),
                         Ainv_gpu.data(), hstream);
      cudaErrorCheck(cudaMemcpyAsync(Binv_gpu.data(), Binv.data(), lda_Binv * delay_count * sizeof(T),
                                     cudaMemcpyHostToDevice, hstream),
                     "cudaMemcpyAsync failed!");
      cuBLAS::gemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, norb, delay_count, delay_count, &cone, V_gpu.data(), norb,
                   Binv_gpu.data(), lda_Binv, &czero, U_gpu.data(), norb);
      cuBLAS::gemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, norb, norb, delay_count, &cminusone, U_gpu.data(), norb,
                   temp_gpu.data(), lda_Binv, &cone, Ainv_gpu.data(), norb);
      delay_count = 0;
      // Ainv is invalid, reset range
      prefetched_range.clear();
    }

    // transfer Ainv_gpu to Ainv and wait till completion
    if (transfer_to_host)
    {
      cudaErrorCheck(cudaMemcpyAsync(Ainv.data(), Ainv_gpu.data(), Ainv.size() * sizeof(T), cudaMemcpyDeviceToHost,
                                     hstream),
                     "cudaMemcpyAsync failed!");
      // no need to wait because : For transfers from device memory to pageable host memory, the function will return only once the copy has completed.
      //waitStream();
    }
  }
};
} // namespace qmcplusplus

#endif // QMCPLUSPLUS_DELAYED_UPDATE_CUDA_H
