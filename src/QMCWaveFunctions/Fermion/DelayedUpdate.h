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

#ifndef QMCPLUSPLUS_DELAYED_UPDATE_H
#define QMCPLUSPLUS_DELAYED_UPDATE_H

#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "CPU/BLAS.hpp"
#include "CPU/BlasThreadingEnv.h"
#include "DiracMatrix.h"
#include "Concurrency/OpenMP.h"

namespace qmcplusplus
{
/** implements delayed update on CPU using BLAS
 * @tparam T base precision for most computation
 * @tparam T_FP high precision for matrix inversion, T_FP >= T
 */
template<typename T, typename T_FP>
class DelayedUpdate
{
  /// orbital values of delayed electrons
  Matrix<T> U;
  /// rows of Ainv corresponding to delayed electrons
  Matrix<T> V;
  /// Matrix inverse of B, at maximum KxK
  Matrix<T> Binv;
  /// scratch space, used during inverse update
  Matrix<T> tempMat;
  /// temporal scratch space used by SM-1
  Vector<T> temp;
  /// new column of B
  Vector<T> p;
  /// list of delayed electrons
  std::vector<int> delay_list;
  /// current number of delays, increase one for each acceptance, reset to 0 after updating Ainv
  int delay_count;
  /// matrix inversion engine
  DiracMatrix<T_FP> detEng;

public:
  /// default constructor
  DelayedUpdate() : delay_count(0) {}

  /** resize the internal storage
   * @param norb number of electrons/orbitals
   * @param delay, maximum delay 0<delay<=norb
   */
  inline void resize(int norb, int delay)
  {
    V.resize(delay, norb);
    U.resize(delay, norb);
    p.resize(delay);
    temp.resize(norb);
    tempMat.resize(norb, delay);
    Binv.resize(delay, delay);
    delay_list.resize(delay);
  }

  /** compute the inverse of the transpose of matrix A
   * @param logdetT orbital value matrix
   * @param Ainv inverse matrix
   */
  template<typename TREAL>
  inline void invert_transpose(const Matrix<T>& logdetT, Matrix<T>& Ainv, std::complex<TREAL>& log_value)
  {
    detEng.invert_transpose(logdetT, Ainv, log_value);
    // safe mechanism
    delay_count = 0;
  }

  /** initialize internal objects when Ainv is refreshed
   * @param Ainv inverse matrix
   */
  inline void initializeInv(const Matrix<T>& Ainv)
  {
    // safe mechanism
    delay_count = 0;
  }

  inline int getDelayCount() const { return delay_count; }

  /** compute the row of up-to-date Ainv
   * @param Ainv inverse matrix
   * @param rowchanged the row id corresponding to the proposed electron
   */
  template<typename VVT>
  inline void getInvRow(const Matrix<T>& Ainv, int rowchanged, VVT& invRow)
  {
    if (delay_count == 0)
    {
      // Ainv is fresh, directly access Ainv
      std::copy_n(Ainv[rowchanged], invRow.size(), invRow.data());
      return;
    }
    constexpr T cone(1);
    constexpr T czero(0);
    const int norb     = Ainv.rows();
    const int lda_Binv = Binv.cols();
    // save Ainv[rowchanged] to invRow
    std::copy_n(Ainv[rowchanged], norb, invRow.data());
    // multiply V (NxK) Binv(KxK) U(KxN) invRow right to the left
    BLAS::gemv('T', norb, delay_count, cone, U.data(), norb, invRow.data(), 1, czero, p.data(), 1);
    BLAS::gemv('N', delay_count, delay_count, -cone, Binv.data(), lda_Binv, p.data(), 1, czero, Binv[delay_count], 1);
    BLAS::gemv('N', norb, delay_count, cone, V.data(), norb, Binv[delay_count], 1, cone, invRow.data(), 1);
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
    constexpr T cone(1);
    constexpr T czero(0);
    const int norb     = Ainv.rows();
    const int lda_Binv = Binv.cols();
    std::copy_n(Ainv[rowchanged], norb, V[delay_count]);
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
      updateInvMat(Ainv);
  }

  /** update the full Ainv and reset delay_count
   * @param Ainv inverse matrix
   */
  inline void updateInvMat(Matrix<T>& Ainv)
  {
    if (delay_count == 0)
      return;
    // update the inverse matrix
    constexpr T cone(1);
    constexpr T czero(0);
    const int norb = Ainv.rows();
    if (delay_count == 1)
    {
      // this is a special case invoking the Fahy's variant of Sherman-Morrison update.
      // Only use the first norb elements of tempMat as a temporal array
      BLAS::gemv('T', norb, norb, cone, Ainv.data(), norb, U[0], 1, czero, temp.data(), 1);
      temp[delay_list[0]] -= cone;
      BLAS::ger(norb, norb, -Binv[0][0], V[0], 1, temp.data(), 1, Ainv.data(), norb);
    }
    else
    {
      const int lda_Binv = Binv.cols();
      // number of threads at the next level, forced to 1 if the problem is small.
      const int num_threads = (norb < 256 ? 1 : getNextLevelNumThreads());
      if (num_threads == 1 || BlasThreadingEnv::NestedThreadingSupported())
      {
        // threading depends on BLAS
        BlasThreadingEnv knob(num_threads);
        BLAS::gemm('T', 'N', delay_count, norb, norb, cone, U.data(), norb, Ainv.data(), norb, czero, tempMat.data(),
                   lda_Binv);
        for (int i = 0; i < delay_count; i++)
          tempMat(delay_list[i], i) -= cone;
        BLAS::gemm('N', 'N', norb, delay_count, delay_count, cone, V.data(), norb, Binv.data(), lda_Binv, czero,
                   U.data(), norb);
        BLAS::gemm('N', 'N', norb, norb, delay_count, -cone, U.data(), norb, tempMat.data(), lda_Binv, cone,
                   Ainv.data(), norb);
      }
      else
      {
        // manually threaded version of the above GEMM calls
#pragma omp parallel
        {
          const int block_size = getAlignedSize<T>((norb + num_threads - 1) / num_threads);
          int num_block        = (norb + block_size - 1) / block_size;
#pragma omp for
          for (int ix = 0; ix < num_block; ix++)
          {
            int x_offset = ix * block_size;
            BLAS::gemm('T', 'N', delay_count, std::min(norb - x_offset, block_size), norb, cone, U.data(), norb,
                       Ainv[x_offset], norb, czero, tempMat[x_offset], lda_Binv);
          }
#pragma omp master
          for (int i = 0; i < delay_count; i++)
            tempMat(delay_list[i], i) -= cone;
#pragma omp for
          for (int iy = 0; iy < num_block; iy++)
          {
            int y_offset = iy * block_size;
            BLAS::gemm('N', 'N', std::min(norb - y_offset, block_size), delay_count, delay_count, cone,
                       V.data() + y_offset, norb, Binv.data(), lda_Binv, czero, U.data() + y_offset, norb);
          }
#pragma omp for collapse(2) nowait
          for (int iy = 0; iy < num_block; iy++)
            for (int ix = 0; ix < num_block; ix++)
            {
              int x_offset = ix * block_size;
              int y_offset = iy * block_size;
              BLAS::gemm('N', 'N', std::min(norb - y_offset, block_size), std::min(norb - x_offset, block_size),
                         delay_count, -cone, U.data() + y_offset, norb, tempMat[x_offset], lda_Binv, cone,
                         Ainv[x_offset] + y_offset, norb);
            }
        }
      }
    }
    delay_count = 0;
  }
};
} // namespace qmcplusplus

#endif // QMCPLUSPLUS_DELAYED_UPDATE_H
