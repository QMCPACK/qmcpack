//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_DIRAC_MATRIX_COMPUTE_CUDA_H
#define QMCPLUSPLUS_DIRAC_MATRIX_COMPUTE_CUDA_H

#include "OhmmsPETE/OhmmsMatrix.h"
#include "OMPTarget/OMPallocator.hpp"
#include "Platforms/PinnedAllocator.h"
#include "Platforms/CUDA/CUDALinearAlgebraHandles.h"
#include "type_traits/scalar_traits.h"
#include "Message/OpenMP.h"
#include "CPU/SIMD/simd.hpp"
#include "ResourceCollection.h"

namespace qmcplusplus
{
/** helper class to compute matrix inversion and the log value of determinant
 *  of a batch of DiracMatrixes.
 * @tparam T_FP the datatype used in the actual computation of matrix inversion
 *  
 *  There is one per crowd not one per MatrixUpdateEngine.
 *  this prevents needing to deal with the resources.
 */
template<typename T_FP>
class DiracMatrixComputeCUDA : public Resource
{
  // Why not just use QMCTraits::FullPrecRealType?
  using FullPrecReal = typename scalar_traits<T_FP>::real_type;

  template<typename T>
  using OffloadPinnedAllocator = OMPallocator<T, PinnedAlignedAllocator<T>>;

  template<typename T>
  using OffloadPinnedMatrix = Matrix<T, OffloadPinnedAllocator<T>>;

  template<typename T>
  using OffloadPinnedVector = Vector<T, OffloadPinnedAllocator<T>>;

  //Unlike for the DirecMatrix.h this is going to be contiguous vectors for each walker of size lda.
  OffloadPinnedVector<int> m_pivots_;
  // Contiguous Matrices for each walker, n^2 elements
  OffloadPinnedVector<T_FP> psiM_fp_;
  OffloadPinnedVector<T_FP> LU_diags_fp_;
  OffloadPinnedVector<T_FP> logdets_fp_;
  OffloadPinnedVector<int> pivots_;
  OffloadPinnedVector<int> infos_;

  Vector<char, OffloadPinnedAllocator<char>> psiM_ptrs_;
  /// reset internal work space
  // inline void reset(RefVector < T_FP * invMat_ptr, const int lda)
  // {
  //   m_pivot.resize(lda);

  //   T_FP tmp;
  //   real_type_fp lw;
  //   Xgetri(lda, invMat_ptr, lda, m_pivots_.data(), &tmp, Lwork);
  //   convert(tmp, lw);
  //   Lwork = static_cast<int>(lw);
  //   m_work.resize(Lwork);
  //   LU_diag.resize(lda);
  // }

  /** Calculates the actual inv and log determinant on accelerator
   */
  template<typename TREAL>
  inline void computeInvertAndLog(cudaStream_t hstream,
                                  OffloadPinnedVector<TREAL>& psi_Ms,
                                  const int n,
                                  const int lda,
                                  OffloadPinnedVector<std::complex<TREAL>>& log_dets)
  {
    int nw = log_dets.size();
    psiM_ptrs_.resize(sizeof(TREAL*) * nw);
    Vector<const T_FP*> ptr_buffer(reinterpret_cast<const TREAL**>(psiM_ptrs_.data()), nw);
    for (int iw = 0; iw < nw; ++iw)
      ptr_buffer[iw] = psi_Ms.device_data() + iw * n * n;
    pivots_.resize(n * nw);
    infos_.resize(nw);
    LU_diags_fp_.resize(n * nw);
    cudaErrorCheck(cudaMemcpyAsync(psiM_ptrs_.device_data(), psiM_ptrs_.data(), psiM_ptrs_.size(),
                                   cudaMemcpyHostToDevice, hstream),
                   "cudaMemcpyAsync psiM_ptrs_ failed!");
    computeInverseAndDetLog_batched(hstream, n, lda, psiM_ptrs_.device_data(), LU_diags_fp_.device_data(),
                                    pivots_.device_data(), infos_.device_data(), log_dets.device_data(), nw);
  }

public:
  template<typename TMAT, typename TREAL>
  inline void mw_invertTranspose(CUDALinearAlgebraHandles& cuda_handles,
                                 const RefVector<OffloadPinnedMatrix<TMAT>>& a_mats,
                                  RefVector<OffloadPinnedMatrix<TMAT>>& inv_a_mats,
                                  OffloadPinnedVector<std::complex<TREAL>>& log_dets)
  {
    const int n   = inv_a_mats[0].get().rows();
    const int lda = inv_a_mats[0].get().cols();
    size_t nsqr{n * n};
    psiM_fp_.resize(n * lda * a_mats.size());
    for (int iw = 0; iw < a_mats.size(); ++iw)
      simd::transpose(a_mats[iw].get().data(), n, a_mats[iw].get().cols(), psiM_fp_.data() + nsqr * iw, n, lda);
    cudaErrorCheck(cudaMemcpyAsync(psiM_fp_.device_data(), psiM_fp_.data(), psiM_fp_.size(), cudaMemcpyHostToDevice,
                                   cuda_handles.hstream),
                   "cudaMemcpyAsync failed copying DiracMatrixBatch::psiM_fp to device");
    computeInvertAndLog(cuda_handles.h_cublas, psiM_fp_, n, lda, log_dets);
    cudaErrorCheck(cudaMemcpyAsync(psiM_fp_.data(), psiM_fp_.device_data(), psiM_fp_.size(), cudaMemcpyDeviceToHost,
                                   cuda_handles.hstream),
                   "cudaMemcpyAsync failed copying back DiracMatrixBatch::psiM_fp from device");
    for (int iw = 0; iw < a_mats.size(); ++iw)
    {
      Matrix<TMAT> data_ref_matrix;
      data_ref_matrix.attachReference(psiM_fp_.data() + nsqr * iw, n, n);
      // Use ohmms matrix to do element wise assignment with possible narrowing conversion.
      inv_a_mats[iw].get() = data_ref_matrix;
    }
  }
};
} // namespace qmcplusplus

#endif //QMCPLUSPLUS_DIRAC_MATRIX_COMPUTE_CUDA_H
