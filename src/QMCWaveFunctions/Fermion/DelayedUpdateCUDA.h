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
#include "CUDA/CUDAallocator.hpp"
#include "CUDA/cuBLAS.hpp"
#include "CUDA/cusolver.hpp"
#include "QMCWaveFunctions/detail/CUDA/delayed_update_helper.h"
#include "QMCWaveFunctions/Fermion/DiracMatrix.h"
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

/** implements delayed update on NVIDIA GPU using cuBLAS and cusolverDN
 * @tparam T base precision for most computation
 * @tparam T_FP high precision for matrix inversion, T_FP >= T
 */
template<typename T, typename T_FP>
class DelayedUpdateCUDA
{
  /// define real type
  using real_type = typename scalar_traits<T>::real_type;
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
  /// scratch memory for cusolverDN
  Matrix<T_FP, CUDAAllocator<T_FP>> Mat1_gpu;
  /** scratch memory for cusolverDN and shared storage with Ainv_gpu.
   * The full precision build has Mat2_gpu = Ainv_gpu
   * The mixed precision build has the first half of Mat2_gpu containing Ainv_gpu
   */
  Matrix<T_FP, CUDAAllocator<T_FP>> Mat2_gpu;
  /// pivot array + info
  Vector<int, CUDAHostAllocator<int>> ipiv;
  Vector<int, CUDAAllocator<int>> ipiv_gpu;
  /// workspace
  Vector<T_FP, CUDAAllocator<T_FP>> work_gpu;
  /// diagonal terms of LU matrix
  Vector<T_FP, CUDAHostAllocator<T_FP>> LU_diag;
  Vector<T_FP, CUDAAllocator<T_FP>> LU_diag_gpu;

  // the range of prefetched_Ainv_rows
  Range prefetched_range;
  // Ainv prefetch buffer
  Matrix<T, CUDAHostAllocator<T>> Ainv_buffer;

  // CUDA specific variables
  cublasHandle_t h_cublas;
  cusolverDnHandle_t h_cusolver;
  cudaStream_t hstream;

  inline void waitStream() { cudaErrorCheck(cudaStreamSynchronize(hstream), "cudaStreamSynchronize failed!"); }

public:
  /// default constructor
  DelayedUpdateCUDA() : delay_count(0)
  {
    cudaErrorCheck(cudaStreamCreate(&hstream), "cudaStreamCreate failed!");
    cublasErrorCheck(cublasCreate(&h_cublas), "cublasCreate failed!");
    cublasErrorCheck(cublasSetStream(h_cublas, hstream), "cublasSetStream failed!");
    cusolverErrorCheck(cusolverDnCreate(&h_cusolver), "cusolverCreate failed!");
    cusolverErrorCheck(cusolverDnSetStream(h_cusolver, hstream), "cusolverSetStream failed!");
  }

  ~DelayedUpdateCUDA()
  {
    cusolverErrorCheck(cusolverDnDestroy(h_cusolver), "cusolverDestroy failed!");
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
    Mat1_gpu.resize(norb, norb);
    Mat2_gpu.resize(norb, norb);
    LU_diag.resize(norb);
    LU_diag_gpu.resize(norb);
    // Ainv_gpu reuses all or part of the memory of Mat2_gpu
    Ainv_gpu.attachReference(reinterpret_cast<T*>(Mat2_gpu.data()), norb, norb);
    // prepare cusolver auxiliary arrays
    ipiv.resize(norb + 1);
    ipiv_gpu.resize(norb + 1);
    int lwork;
    cusolverErrorCheck(cusolver::getrf_bufferSize(h_cusolver, norb, norb, Mat2_gpu.data(), norb, &lwork),
                       "cusolver::getrf_bufferSize failed!");
    work_gpu.resize(lwork);
  }

  /** compute the inverse of the transpose of matrix A and its determinant value in log
   * when T_FP and T are the same
   * @tparam TREAL real type
   */
  template<typename TMAT, typename TREAL, typename = std::enable_if_t<std::is_same<TMAT, T>::value>>
  std::enable_if_t<std::is_same<TMAT, T_FP>::value>
  invert_transpose(const Matrix<TMAT>& logdetT, Matrix<TMAT>& Ainv, std::complex<TREAL>& log_value)
  {
    // safe mechanism
    delay_count = 0;
    int norb    = logdetT.rows();
    cudaErrorCheck(cudaMemcpyAsync(Mat1_gpu.data(), logdetT.data(), logdetT.size() * sizeof(T), cudaMemcpyHostToDevice,
                                   hstream),
                   "cudaMemcpyAsync failed!");
    cusolverErrorCheck(cusolver::getrf(h_cusolver, norb, norb, Mat1_gpu.data(), norb, work_gpu.data(),
                                       ipiv_gpu.data() + 1, ipiv_gpu.data()),
                       "cusolver::getrf failed!");
    cudaErrorCheck(cudaMemcpyAsync(ipiv.data(), ipiv_gpu.data(), ipiv_gpu.size() * sizeof(int), cudaMemcpyDeviceToHost,
                                   hstream),
                   "cudaMemcpyAsync failed!");
    extract_matrix_diagonal_cuda(norb, Mat1_gpu.data(), norb, LU_diag_gpu.data(), hstream);
    cudaErrorCheck(cudaMemcpyAsync(LU_diag.data(), LU_diag_gpu.data(), LU_diag.size() * sizeof(T_FP),
                                   cudaMemcpyDeviceToHost, hstream),
                   "cudaMemcpyAsync failed!");
    // check LU success
    waitStream();
    if (ipiv[0] != 0)
    {
      std::ostringstream err;
      err << "cusolver::getrf calculation failed with devInfo = " << ipiv[0] << std::endl;
      std::cerr << err.str();
      throw std::runtime_error(err.str());
    }
    make_identity_matrix_cuda(norb, Mat2_gpu.data(), norb, hstream);
    cusolverErrorCheck(cusolver::getrs(h_cusolver, CUBLAS_OP_T, norb, norb, Mat1_gpu.data(), norb, ipiv_gpu.data() + 1,
                                       Mat2_gpu.data(), norb, ipiv_gpu.data()),
                       "cusolver::getrs failed!");
    cudaErrorCheck(cudaMemcpyAsync(ipiv.data(), ipiv_gpu.data(), sizeof(int), cudaMemcpyDeviceToHost, hstream),
                   "cudaMemcpyAsync failed!");
    computeLogDet(LU_diag.data(), norb, ipiv.data() + 1, log_value);
    cudaErrorCheck(cudaMemcpyAsync(Ainv.data(), Ainv_gpu.data(), Ainv.size() * sizeof(T), cudaMemcpyDeviceToHost,
                                   hstream),
                   "cudaMemcpyAsync failed!");
    // no need to wait because : For transfers from device memory to pageable host memory, the function will return only once the copy has completed.
    //waitStream();
    if (ipiv[0] != 0)
    {
      std::ostringstream err;
      err << "cusolver::getrs calculation failed with devInfo = " << ipiv[0] << std::endl;
      std::cerr << err.str();
      throw std::runtime_error(err.str());
    }
  }

  /** compute the inverse of the transpose of matrix A and its determinant value in log
   * when T_FP and T are the same
   * @tparam TREAL real type
   */
  template<typename TMAT, typename TREAL, typename = std::enable_if_t<std::is_same<TMAT, T>::value>>
  std::enable_if_t<!std::is_same<TMAT, T_FP>::value>
  invert_transpose(const Matrix<TMAT>& logdetT, Matrix<TMAT>& Ainv, std::complex<TREAL>& log_value)
  {
    // safe mechanism
    delay_count = 0;
    int norb    = logdetT.rows();
    cudaErrorCheck(cudaMemcpyAsync(Mat1_gpu.data(), logdetT.data(), logdetT.size() * sizeof(T), cudaMemcpyHostToDevice,
                                   hstream),
                   "cudaMemcpyAsync failed!");
    copy_matrix_cuda(norb, norb, (T*)Mat1_gpu.data(), norb, Mat2_gpu.data(), norb, hstream);
    cusolverErrorCheck(cusolver::getrf(h_cusolver, norb, norb, Mat2_gpu.data(), norb, work_gpu.data(),
                                       ipiv_gpu.data() + 1, ipiv_gpu.data()),
                       "cusolver::getrf failed!");
    cudaErrorCheck(cudaMemcpyAsync(ipiv.data(), ipiv_gpu.data(), ipiv_gpu.size() * sizeof(int), cudaMemcpyDeviceToHost,
                                   hstream),
                   "cudaMemcpyAsync failed!");
    extract_matrix_diagonal_cuda(norb, Mat2_gpu.data(), norb, LU_diag_gpu.data(), hstream);
    cudaErrorCheck(cudaMemcpyAsync(LU_diag.data(), LU_diag_gpu.data(), LU_diag.size() * sizeof(T_FP),
                                   cudaMemcpyDeviceToHost, hstream),
                   "cudaMemcpyAsync failed!");
    // check LU success
    waitStream();
    if (ipiv[0] != 0)
    {
      std::ostringstream err;
      err << "cusolver::getrf calculation failed with devInfo = " << ipiv[0] << std::endl;
      std::cerr << err.str();
      throw std::runtime_error(err.str());
    }
    make_identity_matrix_cuda(norb, Mat1_gpu.data(), norb, hstream);
    cusolverErrorCheck(cusolver::getrs(h_cusolver, CUBLAS_OP_T, norb, norb, Mat2_gpu.data(), norb, ipiv_gpu.data() + 1,
                                       Mat1_gpu.data(), norb, ipiv_gpu.data()),
                       "cusolver::getrs failed!");
    copy_matrix_cuda(norb, norb, Mat1_gpu.data(), norb, (T*)Mat2_gpu.data(), norb, hstream);
    cudaErrorCheck(cudaMemcpyAsync(ipiv.data(), ipiv_gpu.data(), sizeof(int), cudaMemcpyDeviceToHost, hstream),
                   "cudaMemcpyAsync failed!");
    computeLogDet(LU_diag.data(), norb, ipiv.data() + 1, log_value);
    cudaErrorCheck(cudaMemcpyAsync(Ainv.data(), Ainv_gpu.data(), Ainv.size() * sizeof(T), cudaMemcpyDeviceToHost,
                                   hstream),
                   "cudaMemcpyAsync failed!");
    // no need to wait because : For transfers from device memory to pageable host memory, the function will return only once the copy has completed.
    //waitStream();
    if (ipiv[0] != 0)
    {
      std::ostringstream err;
      err << "cusolver::getrs calculation failed with devInfo = " << ipiv[0] << std::endl;
      std::cerr << err.str();
      throw std::runtime_error(err.str());
    }
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
