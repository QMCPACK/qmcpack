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

#ifndef QMCPLUSPLUS_DELAYED_UPDATE_SYCL_H
#define QMCPLUSPLUS_DELAYED_UPDATE_SYCL_H

#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "SYCL/SYCLruntime.hpp"
#include "SYCL/SYCLallocator.hpp"
#include "QMCWaveFunctions/detail/SYCL/sycl_determinant_helper.hpp"
#include "DiracMatrix.h"
#include "oneapi/mkl/blas.hpp"
//#define SYCL_BLOCKING

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
class DelayedUpdateSYCL
{
  // Data staged during for delayed acceptRows
  Matrix<T> U;
  Matrix<T> Binv;
  Matrix<T> V;
  //Matrix<T> tempMat; // for debugging only
  Matrix<T, SYCLAllocator<T>> temp_gpu;
  /// GPU copy of U, V, Binv, Ainv
  Matrix<T, SYCLAllocator<T>> U_gpu;
  Matrix<T, SYCLAllocator<T>> V_gpu;
  Matrix<T, SYCLAllocator<T>> Binv_gpu;
  Matrix<T, SYCLAllocator<T>> Ainv_gpu;
  // auxiliary arrays for B
  Vector<T> p;
  // using host allocator
  Vector<int, SYCLHostAllocator<int>> delay_list;
  /// current number of delays, increase one for each acceptance, reset to 0 after updating Ainv
  int delay_count;

  DiracMatrix<T_FP> host_inverter_;

  // the range of prefetched_Ainv_rows
  Range prefetched_range;
  // Ainv prefetch buffer
  Matrix<T> Ainv_buffer;

  sycl::queue *m_queue=nullptr;
  sycl::event live_event;

  /// reset delay count to 0
  inline void clearDelayCount()
  {
    delay_count = 0;
    prefetched_range.clear();
  }

public:
  /// default constructor
  DelayedUpdateSYCL() : delay_count(0)
  {
    m_queue=get_default_queue();
  }

  ~DelayedUpdateSYCL()
  {
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
    //delay_list_gpu.resize(delay);
    Ainv_gpu.resize(norb, norb);

  }

  /** compute the inverse of the transpose of matrix A and its determinant value in log
   * @tparam TREAL real type
   */
  template<typename TREAL>
  void invert_transpose(const Matrix<T>& logdetT, Matrix<T>& Ainv, std::complex<TREAL>& log_value)
  {
    host_inverter_.invert_transpose(logdetT, Ainv, log_value);
    initializeInv(Ainv);
  }

  /** initialize internal objects when Ainv is refreshed
   * @param Ainv inverse matrix
   */
  inline void initializeInv(const Matrix<T>& Ainv)
  {
#ifdef SYCL_BLOCKING
    m_queue->memcpy(Ainv_gpu.data(), Ainv.data(), Ainv.size() * sizeof(T)).wait();
#else
    live_event = m_queue->memcpy(Ainv_gpu.data(), Ainv.data(), Ainv.size() * sizeof(T));
#endif
    clearDelayCount();
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
#ifdef SYCL_BLOCKING
      m_queue->memcpy(Ainv_buffer.data(), Ainv_gpu[rowchanged],
                      invRow.size() * (last_row - rowchanged) * sizeof(T)).wait();
#else
      live_event = m_queue->memcpy(Ainv_buffer.data(), Ainv_gpu[rowchanged],
                                   invRow.size() * (last_row - rowchanged) * sizeof(T), live_event);
#endif
      prefetched_range.setRange(rowchanged, last_row);
    }

#ifndef SYCL_BLOCKING
    sycl::event::wait({live_event});
#endif

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
    {
      updateInvMat(Ainv, false);
    }
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
      constexpr auto trans    = oneapi::mkl::transpose::trans;
      constexpr auto nontrans = oneapi::mkl::transpose::nontrans;

      const int norb     = Ainv.rows();
      const int lda_Binv = Binv.cols();

#if 0
      m_queue->memcpy(U_gpu.data(), U.data(), norb * delay_count * sizeof(T)).wait();

      oneapi::mkl::blas::gemm(*m_queue, trans, nontrans, 
                              delay_count, norb, norb, cone, U_gpu.data(),
                              norb, Ainv_gpu.data(), norb, czero, temp_gpu.data(), 
                              lda_Binv).wait();

      m_queue->memcpy(Binv_gpu.data(), Binv.data(), lda_Binv * delay_count * sizeof(T)).wait();

      applyW_stageV_sycl(*m_queue, delay_list.data(), delay_count, 
                         temp_gpu.data(), norb, temp_gpu.cols(), V_gpu.data(), Ainv_gpu.data()).wait();

      oneapi::mkl::blas::gemm(*m_queue, nontrans, nontrans, norb, delay_count, delay_count, cone,
                                   V_gpu.data(), norb, Binv_gpu.data(), lda_Binv, czero, U_gpu.data(), norb).wait();
      oneapi::mkl::blas::gemm(*m_queue, nontrans, nontrans, norb, norb, delay_count, -cone,
                              U_gpu.data(), norb, temp_gpu.data(), lda_Binv, cone, Ainv_gpu.data(), norb).wait();

#endif

#if 1
      auto u_ = m_queue->memcpy(U_gpu.data(), U.data(), norb * delay_count * sizeof(T));
      auto b_ = m_queue->memcpy(Binv_gpu.data(), Binv.data(), lda_Binv * delay_count * sizeof(T));

      oneapi::mkl::blas::gemm(*m_queue, trans, nontrans, 
                              delay_count, norb, norb, cone, U_gpu.data(),
                              norb, Ainv_gpu.data(), norb, czero, temp_gpu.data(), 
                              lda_Binv, {u_} ).wait();

      u_ = applyW_stageV_sycl(*m_queue, delay_list.data(), delay_count, 
                              temp_gpu.data(), norb, temp_gpu.cols(), V_gpu.data(), Ainv_gpu.data());

      b_ = oneapi::mkl::blas::gemm(*m_queue, nontrans, nontrans, norb, delay_count, delay_count, cone,
                                   V_gpu.data(), norb, Binv_gpu.data(), lda_Binv, czero, U_gpu.data(), norb, {b_, u_} );

#ifdef SYCL_BLOCKING
      oneapi::mkl::blas::gemm(*m_queue, nontrans, nontrans, norb, norb, delay_count, -cone,
                              U_gpu.data(), norb, temp_gpu.data(), lda_Binv, cone, Ainv_gpu.data(), norb, {b_} ).wait();
#else
      live_event = oneapi::mkl::blas::gemm(*m_queue, nontrans, nontrans, norb, norb, delay_count, -cone,
                                           U_gpu.data(), norb, temp_gpu.data(), lda_Binv, cone, Ainv_gpu.data(), norb, {b_} );
#endif

#endif

      clearDelayCount();
    }

    // transfer Ainv_gpu to Ainv and wait till completion
    if (transfer_to_host)
    {
#ifdef SYCL_BLOCKING
      m_queue->memcpy(Ainv.data(), Ainv_gpu.data(), Ainv.size() * sizeof(T)).wait();
#else
      m_queue->memcpy(Ainv.data(), Ainv_gpu.data(), Ainv.size() * sizeof(T), {live_event}).wait();
#endif
    }
  }
};
} // namespace qmcplusplus

#endif // QMCPLUSPLUS_DELAYED_UPDATE_CUDA_H
