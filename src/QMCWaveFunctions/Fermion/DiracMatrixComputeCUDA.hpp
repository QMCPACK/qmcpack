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
#include "Platforms/CUDA/cuBLAS.hpp"
#include "Platforms/CUDA/cuBLAS_LU.hpp"
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
  OffloadPinnedVector<T_FP> invM_fp_;
  OffloadPinnedVector<T_FP> LU_diags_fp_;
  OffloadPinnedVector<T_FP> logdets_fp_;
  OffloadPinnedVector<int> pivots_;
  OffloadPinnedVector<int> infos_;

  Vector<char, OffloadPinnedAllocator<char>> psiM_ptrs_;
  Vector<char, OffloadPinnedAllocator<char>> invM_ptrs_;

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
   *
   *  \param[in]      h_cublas    cublas handle, hstream handle is retrieved from it.			
   *  \param[in,out]  psi_Ms      matrices flattened into single pinned vector, returned with LU matrices.
   *  \param[out]     inv_Ms      matrices flattened into single pinned vector.				
   *  \param[in]      n           matrices rank.								
   *  \param[in]      lda         leading dimension of each matrix					
   *  \param[out]     log_values  log determinant value for each matrix, batch_size = log_values.size()
   *
   *  It's essential that the hstream be the same used by h_cublas due to the async transfer and
   *  kernel dependencies.
   */
  template<typename TREAL>
  inline void computeInvertAndLog(cublasHandle_t h_cublas,
                                  OffloadPinnedVector<TREAL>& psi_Ms,
                                  OffloadPinnedVector<TREAL>& inv_Ms,
                                  const int n,
                                  const int lda,
                                  OffloadPinnedVector<std::complex<TREAL>>& log_values)
  {
    // This is probably dodgy
    int nw = log_values.size();
    psiM_ptrs_.resize(sizeof(TREAL*) * nw);
    invM_ptrs_.resize(sizeof(TREAL*) * nw);
    Vector<const T_FP*> M_ptr_buffer(reinterpret_cast<const TREAL**>(psiM_ptrs_.data()), nw);
    Vector<const T_FP*> invM_ptr_buffer(reinterpret_cast<const TREAL**>(invM_ptrs_.data()), nw);
    for (int iw = 0; iw < nw; ++iw)
    {
      M_ptr_buffer[iw]    = psi_Ms.device_data() + iw * n * lda;
      invM_ptr_buffer[iw] = inv_Ms.device_data() + iw * n * lda;
    }
    pivots_.resize(n * nw);
    infos_.resize(nw);
    LU_diags_fp_.resize(n * nw);
    cudaStream_t hstream;
    cublasErrorCheck(cublasGetStream(h_cublas, &hstream), "cublasGetStream failed!");
    cudaErrorCheck(cudaMemcpyAsync(psi_Ms.device_data(), psi_Ms.data(), psi_Ms.size() * sizeof(T_FP),
                                   cudaMemcpyHostToDevice, hstream),
                   "cudaMemcpyAsync failed copying DiracMatrixBatch::psiM_fp to device");
    cudaErrorCheck(cudaMemcpyAsync(psiM_ptrs_.device_data(), psiM_ptrs_.data(), psiM_ptrs_.size(),
                                   cudaMemcpyHostToDevice, hstream),
                   "cudaMemcpyAsync psiM_ptrs_ failed!");
    cudaErrorCheck(cudaMemcpyAsync(invM_ptrs_.device_data(), invM_ptrs_.data(), invM_ptrs_.size(),
                                   cudaMemcpyHostToDevice, hstream),
                   "cudaMemcpyAsync invM_ptrs_ failed!");
    TREAL** psiM_mw_ptr = reinterpret_cast<TREAL**>(psiM_ptrs_.device_data());
    TREAL** invM_mw_ptr = reinterpret_cast<TREAL**>(invM_ptrs_.device_data());

    cuBLAS_LU::computeInverseAndDetLog_batched(h_cublas, hstream, n, lda, psiM_mw_ptr, invM_mw_ptr,
                                               LU_diags_fp_.device_data(), pivots_.device_data(), infos_.device_data(),
                                               log_values.device_data(), nw);
#if NDEBUG
    // This is very useful to see whether the data after all kernels and cublas calls are run is wrong on the device or due to copy.
    // cuBLAS_LU::peekinvM_batched(hstream, psiM_mw_ptr, invM_mw_ptr, pivots_.device_data(), infos_.device_data(),
    //                             log_values.device_data(), nw);
#endif
    cudaErrorCheck(cudaMemcpyAsync(invM_ptrs_.data(), invM_ptrs_.device_data(), invM_ptrs_.size(),
                                   cudaMemcpyDeviceToHost, hstream),
                   "cudaMemcpyAsync invM_ptrs_ failed!");
    cudaErrorCheck(cudaMemcpyAsync(inv_Ms.data(), inv_Ms.device_data(), inv_Ms.size() * sizeof(TREAL),
                                   cudaMemcpyDeviceToHost, hstream),
                   "cudaMemcpyAsync failed copying back DiracMatrixBatch::invM_fp from device");
    cudaErrorCheck(cudaMemcpyAsync(log_values.data(), log_values.device_data(),
                                   log_values.size() * sizeof(std::complex<TREAL>), cudaMemcpyDeviceToHost, hstream),
                   "cudaMemcpyAsync log_values failed!");
    cudaErrorCheck(cudaStreamSynchronize(hstream), "cudaStreamSynchronize failed!");
  }

public:
  DiracMatrixComputeCUDA() : Resource("DiracMatrixComputeCUDA") {}

  Resource* makeClone() const override { return new DiracMatrixComputeCUDA(*this); }

  /** Given a_mat returns inverted amit and log determinant of a_matches.
   *  \param [in] a_mat a matrix input
   *  \param [out] inv_a_mat inverted matrix
   *  \param [out] log determinant is in logvalues[0]
   *
   *  I consider this single call to be semi depricated so the log determinant values
   *  vector is used to match the primary batched interface to the accelerated routings.
   */
  template<typename TMAT, typename TREAL>
  void invert_transpose(CUDALinearAlgebraHandles& cuda_handles,
                        const OffloadPinnedMatrix<TMAT>& a_mat,
                        OffloadPinnedMatrix<TMAT>& inv_a_mat,
                        OffloadPinnedVector<std::complex<TREAL>>& log_values)
  {
    const int n   = a_mat.rows();
    const int lda = a_mat.cols();
    psiM_fp_.resize(n * lda);
    invM_fp_.resize(n * lda);
    std::fill(log_values.begin(), log_values.end(), std::complex<TREAL>{0.0,0.0});
    cudaErrorCheck(cudaMemcpyAsync(log_values.device_data(), log_values.data(), log_values.size() * sizeof(std::complex<TREAL>),
                                   cudaMemcpyHostToDevice, cuda_handles.hstream),
                   "cudaMemcpyAsync failed copying DiracMatrixBatch::log_values to device");
    simd::transpose(a_mat.data(), n, lda, psiM_fp_.data(), n, lda);
    cudaErrorCheck(cudaMemcpyAsync(psiM_fp_.device_data(), psiM_fp_.data(), psiM_fp_.size() * sizeof(T_FP),
                                   cudaMemcpyHostToDevice, cuda_handles.hstream),
                   "cudaMemcpyAsync failed copying DiracMatrixBatch::psiM_fp to device");
    computeInvertAndLog(cuda_handles.h_cublas, psiM_fp_, invM_fp_, n, lda, log_values);
    OffloadPinnedMatrix<T_FP> data_ref_matrix;
    data_ref_matrix.attachReference(invM_fp_.data(), n, n);
    // Use ohmms matrix to do element wise assignment with possible narrowing conversion.
    inv_a_mat = data_ref_matrix;
    cudaErrorCheck(cudaMemcpyAsync(inv_a_mat.device_data(), inv_a_mat.data(), inv_a_mat.size()*sizeof(TMAT), cudaMemcpyHostToDevice, cuda_handles.hstream), "cudaMemcpyAsync of inv_a_mat to device failed!");
  }

  /** as it stands there is no point in a_mats and inv_a_mats being Pinned since they aren't 
   *  used directly.
   */
  template<typename TMAT, typename TREAL>
  inline void mw_invertTranspose(CUDALinearAlgebraHandles& cuda_handles,
                                 const RefVector<const OffloadPinnedMatrix<TMAT>>& a_mats,
                                 RefVector<OffloadPinnedMatrix<TMAT>>& inv_a_mats,
                                 OffloadPinnedVector<std::complex<TREAL>>& log_values)
  {
    int nw        = a_mats.size();
    const int n   = inv_a_mats[0].get().rows();
    const int lda = inv_a_mats[0].get().cols();
    size_t nsqr   = n * n;
    psiM_fp_.resize(n * lda * nw);
    invM_fp_.resize(n * lda * nw);
    std::fill(log_values.begin(), log_values.end(), std::complex<TREAL>{0.0,0.0});
    cudaErrorCheck(cudaMemcpyAsync(log_values.device_data(), log_values.data(), log_values.size() * sizeof(std::complex<TREAL>),
                                   cudaMemcpyHostToDevice, cuda_handles.hstream),
                                  "cudaMemcpyAsync failed copying DiracMatrixBatch::log_values to device");
    for (int iw = 0; iw < nw; ++iw)
      simd::transpose(a_mats[iw].get().data(), n, a_mats[iw].get().cols(), psiM_fp_.data() + nsqr * iw, n, lda);
    computeInvertAndLog(cuda_handles.h_cublas, psiM_fp_, invM_fp_, n, lda, log_values);
    for (int iw = 0; iw < a_mats.size(); ++iw)
    {
      OffloadPinnedMatrix<T_FP> data_ref_matrix;
      data_ref_matrix.attachReference(invM_fp_.data() + nsqr * iw, n, n);
      // Use ohmms matrix to do element wise assignment with possible narrowing conversion.
      inv_a_mats[iw].get() = data_ref_matrix;
      cudaErrorCheck(cudaMemcpyAsync(inv_a_mats[iw].get().device_data(), inv_a_mats[iw].get().data(), inv_a_mats[iw].get().size() * sizeof(TMAT), cudaMemcpyHostToDevice, cuda_handles.hstream), "cudaMemcpyAsync of inv_a_mat to device failed!");
    }
  }
};
} // namespace qmcplusplus

#endif //QMCPLUSPLUS_DIRAC_MATRIX_COMPUTE_CUDA_H
