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

#include <type_traits>

#include "OhmmsPETE/OhmmsMatrix.h"
#include "DualAllocatorAliases.hpp"
#include "Platforms/CUDA/CUDALinearAlgebraHandles.h"
#include "Platforms/CUDA/cuBLAS.hpp"
#include "detail/CUDA/cuBLAS_LU.hpp"
#include "type_traits/complex_help.hpp"
#include "Concurrency/OpenMP.h"
#include "CPU/SIMD/simd.hpp"
#include "ResourceCollection.h"

namespace qmcplusplus
{
/** class defining a compute and memory resource to compute matrix inversion and 
 *  the log determinants of a batch of DiracMatrixes.
 *  Multiplicty is one per crowd not one per UpdateEngine
 *  It matches the multiplicity of the accelerator call 
 *  and batched resource requirement.
 *
 *  @tparam VALUE_FP the datatype used in the actual computation of matrix inversion
 *
 *  There are no per walker variables, resources specific to the per crowd
 *  compute object are owned here. The compute object itself is the resource
 *  to the per walker DiracDeterminantBatched.
 *  Resources used by this object but owned by the 
 *  surrounding scope are passed as arguments.
 */
template<typename VALUE_FP>
class DiracMatrixComputeCUDA : public Resource
{
  using FullPrecReal = RealAlias<VALUE_FP>;
  using LogValue     = std::complex<FullPrecReal>;

  template<typename T>
  using DualMatrix = Matrix<T, PinnedDualAllocator<T>>;

  template<typename T>
  using DualVector = Vector<T, PinnedDualAllocator<T>>;

  // Contiguous memory for fp precision Matrices for each walker, n^2 * nw elements
  DualVector<VALUE_FP> psiM_fp_;
  DualVector<VALUE_FP> invM_fp_;

  // working vectors
  DualVector<VALUE_FP> LU_diags_fp_;
  DualVector<int> pivots_;
  DualVector<int> infos_;

  //DualMatrix<T_FP> temp_mat_;

  /** Transfer buffer for device pointers to matrices.
   *  The element count is usually low and the transfer launch cost are more than the transfer themselves.
   *  For this reason, it is beneficial to fusing multiple lists of pointers.
   *  Right now this buffer packs nw psiM pointers and then packs nw invM pointers.
   *  Use only within a function scope and do not rely on previous value.
   */
  DualVector<VALUE_FP*> psiM_invM_ptrs_;

  // cuBLAS geam wants these.
  VALUE_FP host_one{1.0};
  VALUE_FP host_zero{0.0};

  /** Calculates the actual inv and log determinant on accelerator
   *
   *  \param[in]      h_cublas    cublas handle, hstream handle is retrieved from it.			
   *  \param[in,out]  a_mats      dual A matrices, they will be transposed on the device side as a side effect.
   *  \param[out]     inv_a_mats  dual invM matrices
   *  \param[in]      n           matrices rank.								
   *  \param[out]     log_values  log determinant value for each matrix, batch_size = log_values.size()
   *
   *  On Volta so far little seems to be achieved by having the mats continuous.
   *
   *  List of operations:
   *  1. matrix-by-matrix. Copy a_mat to inv_a_mat on host, transfer inv_a_mat to device, transpose inv_a_mat to a_mat on device.
   *  2. batched. LU and invert
   *  3. matrix-by-matrix. Transfer inv_a_mat to host
   *
   *  Pros and cons:
   *  1. \todo try to do like mw_computeInvertAndLog_stride, copy and transpose to psiM_fp_ and fuse transfer.
   *  3. \todo Remove Transfer inv_a_mat to host and let the upper level code handle it.
   */
  inline void mw_computeInvertAndLog(CUDALinearAlgebraHandles& cuda_handles,
                                     const RefVector<const DualMatrix<VALUE_FP>>& a_mats,
                                     const RefVector<DualMatrix<VALUE_FP>>& inv_a_mats,
                                     const int n,
                                     DualVector<LogValue>& log_values)
  {
    const int nw = a_mats.size();
    assert(a_mats.size() == inv_a_mats.size());

    psiM_invM_ptrs_.resize(nw * 2);
    const int lda           = a_mats[0].get().cols();
    const int ldinv         = inv_a_mats[0].get().cols();
    cudaStream_t hstream    = cuda_handles.hstream;
    cublasHandle_t h_cublas = cuda_handles.h_cublas;
    psiM_fp_.resize(n * ldinv * nw);

    for (int iw = 0; iw < nw; ++iw)
    {
      psiM_invM_ptrs_[iw]      = psiM_fp_.device_data() + iw * n * ldinv;
      psiM_invM_ptrs_[iw + nw] = inv_a_mats[iw].get().device_data();
      // Since inv_a_mat can have a different leading dimension from a_mat first we remap copy on the host
      simd::remapCopy(n, n, a_mats[iw].get().data(), lda, inv_a_mats[iw].get().data(), ldinv);
      // Then copy a_mat in inv_a_mats to the device
      cudaErrorCheck(cudaMemcpyAsync(inv_a_mats[iw].get().device_data(), inv_a_mats[iw].get().data(),
                                     inv_a_mats[iw].get().size() * sizeof(VALUE_FP), cudaMemcpyHostToDevice, hstream),
                     "cudaMemcpyAsync failed copying DiracMatrixBatch::psiM to device");
      // On the device Here we transpose to a_mat;
      cublasErrorCheck(cuBLAS::geam(h_cublas, CUBLAS_OP_T, CUBLAS_OP_N, n, n, &host_one,
                                    inv_a_mats[iw].get().device_data(), ldinv, &host_zero,
                                    a_mats[iw].get().device_data(), lda, psiM_invM_ptrs_[iw], ldinv),
                       "cuBLAS::geam failed.");
    }
    pivots_.resize(n * nw);
    infos_.resize(nw);
    LU_diags_fp_.resize(n * nw);
    cudaErrorCheck(cudaMemcpyAsync(psiM_invM_ptrs_.device_data(), psiM_invM_ptrs_.data(),
                                   psiM_invM_ptrs_.size() * sizeof(VALUE_FP*), cudaMemcpyHostToDevice, hstream),
                   "cudaMemcpyAsync psiM_invM_ptrs_ failed!");
    cuBLAS_LU::computeInverseAndDetLog_batched(h_cublas, hstream, n, ldinv, psiM_invM_ptrs_.device_data(),
                                               psiM_invM_ptrs_.device_data() + nw, LU_diags_fp_.device_data(),
                                               pivots_.device_data(), infos_.data(), infos_.device_data(),
                                               log_values.device_data(), nw);
    for (int iw = 0; iw < nw; ++iw)
    {
      cudaErrorCheck(cudaMemcpyAsync(inv_a_mats[iw].get().data(), inv_a_mats[iw].get().device_data(),
                                     inv_a_mats[iw].get().size() * sizeof(VALUE_FP), cudaMemcpyDeviceToHost, hstream),
                     "cudaMemcpyAsync failed copying DiracMatrixBatch::inv_psiM to host");
    }
    cudaErrorCheck(cudaMemcpyAsync(log_values.data(), log_values.device_data(), log_values.size() * sizeof(LogValue),
                                   cudaMemcpyDeviceToHost, hstream),
                   "cudaMemcpyAsync log_values failed!");
    cudaErrorCheck(cudaStreamSynchronize(hstream), "cudaStreamSynchronize failed!");
  }


  /** Calculates the actual inv and log determinant on accelerator with psiMs and invMs widened to full precision
   *  and copied into continuous vectors.
   *
   *  \param[in]      h_cublas    cublas handle, hstream handle is retrieved from it.			
   *  \param[in,out]  psi_Ms      matrices flattened into single pinned vector, returned with LU matrices.
   *  \param[out]     inv_Ms      matrices flattened into single pinned vector.				
   *  \param[in]      n           matrices rank.								
   *  \param[in]      lda         leading dimension of each matrix					
   *  \param[out]     log_values  log determinant value for each matrix, batch_size = log_values.size()
   *
   *  List of operations:
   *  1. batched. Transfer psi_Ms to device
   *  2. batched. LU and invert
   *  3. batched. Transfer inv_Ms to host
   *  \todo Remove 1 and 3. Handle transfer at upper level.
   */
  inline void mw_computeInvertAndLog_stride(CUDALinearAlgebraHandles& cuda_handles,
                                            DualVector<VALUE_FP>& psi_Ms,
                                            DualVector<VALUE_FP>& inv_Ms,
                                            const int n,
                                            const int lda,
                                            DualVector<LogValue>& log_values)
  {
    // This is probably dodgy
    const int nw = log_values.size();
    psiM_invM_ptrs_.resize(nw * 2);
    for (int iw = 0; iw < nw; ++iw)
    {
      psiM_invM_ptrs_[iw]      = psi_Ms.device_data() + iw * n * lda;
      psiM_invM_ptrs_[iw + nw] = inv_Ms.device_data() + iw * n * lda;
    }
    pivots_.resize(n * nw);
    infos_.resize(nw);
    LU_diags_fp_.resize(n * nw);

    cudaStream_t hstream    = cuda_handles.hstream;
    cublasHandle_t h_cublas = cuda_handles.h_cublas;
    cudaErrorCheck(cudaMemcpyAsync(psi_Ms.device_data(), psi_Ms.data(), psi_Ms.size() * sizeof(VALUE_FP),
                                   cudaMemcpyHostToDevice, hstream),
                   "cudaMemcpyAsync failed copying DiracMatrixBatch::psiM_fp to device");
    cudaErrorCheck(cudaMemcpyAsync(psiM_invM_ptrs_.device_data(), psiM_invM_ptrs_.data(),
                                   psiM_invM_ptrs_.size() * sizeof(VALUE_FP*), cudaMemcpyHostToDevice, hstream),
                   "cudaMemcpyAsync psiM_invM_ptrs_ failed!");
    cuBLAS_LU::computeInverseAndDetLog_batched(h_cublas, hstream, n, lda, psiM_invM_ptrs_.device_data(),
                                               psiM_invM_ptrs_.device_data() + nw, LU_diags_fp_.device_data(),
                                               pivots_.device_data(), infos_.data(), infos_.device_data(),
                                               log_values.device_data(), nw);
#if NDEBUG
    // This is very useful to see whether the data after all kernels and cublas calls are run is wrong on the device or due to copy.
    // cuBLAS_LU::peekinvM_batched(hstream, psiM_mw_ptr, invM_mw_ptr, pivots_.device_data(), infos_.device_data(),
    //                             log_values.device_data(), nw);
#endif
    cudaErrorCheck(cudaMemcpyAsync(inv_Ms.data(), inv_Ms.device_data(), inv_Ms.size() * sizeof(VALUE_FP),
                                   cudaMemcpyDeviceToHost, hstream),
                   "cudaMemcpyAsync failed copying back DiracMatrixBatch::invM_fp from device");
    cudaErrorCheck(cudaMemcpyAsync(log_values.data(), log_values.device_data(), log_values.size() * sizeof(LogValue),
                                   cudaMemcpyDeviceToHost, hstream),
                   "cudaMemcpyAsync log_values failed!");
    cudaErrorCheck(cudaStreamSynchronize(hstream), "cudaStreamSynchronize failed!");
  }

public:
  DiracMatrixComputeCUDA() : Resource("DiracMatrixComputeCUDA") {}

  DiracMatrixComputeCUDA(const DiracMatrixComputeCUDA& other) : Resource(other.getName()) {}

  Resource* makeClone() const override { return new DiracMatrixComputeCUDA(*this); }

  /** Given a_mat returns inverted amit and log determinant of a_matches.
   *  \param [in] a_mat a matrix input
   *  \param [out] inv_a_mat inverted matrix
   *  \param [out] log determinant is in logvalues[0]
   *
   *  I consider this single call to be semi depricated so the log determinant values
   *  vector is used to match the primary batched interface to the accelerated routings.
   *  There is no optimization (yet) for TMAT same type as TREAL
   */
  template<typename TMAT>
  void invert_transpose(CUDALinearAlgebraHandles& cuda_handles,
                        DualMatrix<TMAT>& a_mat,
                        DualMatrix<TMAT>& inv_a_mat,
                        DualVector<LogValue>& log_values)
  {
    const int n   = a_mat.rows();
    const int lda = a_mat.cols();
    psiM_fp_.resize(n * lda);
    invM_fp_.resize(n * lda);
    std::fill(log_values.begin(), log_values.end(), LogValue{0.0, 0.0});
    // making sure we know the log_values are zero'd on the device.
    cudaErrorCheck(cudaMemcpyAsync(log_values.device_data(), log_values.data(), log_values.size() * sizeof(LogValue),
                                   cudaMemcpyHostToDevice, cuda_handles.hstream),
                   "cudaMemcpyAsync failed copying DiracMatrixBatch::log_values to device");
    simd::transpose(a_mat.data(), n, lda, psiM_fp_.data(), n, lda);
    cudaErrorCheck(cudaMemcpyAsync(psiM_fp_.device_data(), psiM_fp_.data(), psiM_fp_.size() * sizeof(VALUE_FP),
                                   cudaMemcpyHostToDevice, cuda_handles.hstream),
                   "cudaMemcpyAsync failed copying DiracMatrixBatch::psiM_fp to device");
    mw_computeInvertAndLog_stride(cuda_handles, psiM_fp_, invM_fp_, n, lda, log_values);
    DualMatrix<VALUE_FP> data_ref_matrix;

    data_ref_matrix.attachReference(invM_fp_.data(), n, n);

    // We can't use operator= with different lda, ldb which can happen so we use this assignment which is over the
    // smaller of the two's dimensions
    inv_a_mat.assignUpperLeft(data_ref_matrix);
    cudaErrorCheck(cudaMemcpyAsync(inv_a_mat.device_data(), inv_a_mat.data(), inv_a_mat.size() * sizeof(TMAT),
                                   cudaMemcpyHostToDevice, cuda_handles.hstream),
                   "cudaMemcpyAsync of inv_a_mat to device failed!");
  }

  /** Mixed precision specialization
   *  When TMAT is not full precision we need to still do the inversion and log
   *  at full precision. This is not yet optimized to transpose on the GPU
   *
   *  List of operations:
   *  1. matrix-by-matrix. Transpose a_mat to psiM_fp_ used on host
   *  2. batched. Call mw_computeInvertAndLog_stride, H2D, invert, D2H
   *  3. matrix-by-matrix. Copy invM_fp_ to inv_a_mat on host. Transfer inv_a_mat to device.
   *
   *  Pros and cons:
   *  1. transfer is batched but double the transfer size due to precision promotion
   *  3. \todo Copy invM_fp_ to inv_a_mat on device is desired. Transfer inv_a_mat to host should be handled by the upper level code.
   */
  template<typename TMAT>
  inline std::enable_if_t<!std::is_same<VALUE_FP, TMAT>::value> mw_invertTranspose(
      CUDALinearAlgebraHandles& cuda_handles,
      const RefVector<const DualMatrix<TMAT>>& a_mats,
      const RefVector<DualMatrix<TMAT>>& inv_a_mats,
      DualVector<LogValue>& log_values)
  {
    assert(log_values.size() == a_mats.size());
    const int nw  = a_mats.size();
    const int n   = a_mats[0].get().rows();
    const int lda = a_mats[0].get().cols();
    size_t nsqr   = n * n;
    psiM_fp_.resize(n * lda * nw);
    invM_fp_.resize(n * lda * nw);
    std::fill(log_values.begin(), log_values.end(), LogValue{0.0, 0.0});
    // making sure we know the log_values are zero'd on the device.
    cudaErrorCheck(cudaMemcpyAsync(log_values.device_data(), log_values.data(), log_values.size() * sizeof(LogValue),
                                   cudaMemcpyHostToDevice, cuda_handles.hstream),
                   "cudaMemcpyAsync failed copying DiracMatrixBatch::log_values to device");
    for (int iw = 0; iw < nw; ++iw)
      simd::transpose(a_mats[iw].get().data(), n, a_mats[iw].get().cols(), psiM_fp_.data() + nsqr * iw, n, lda);
    mw_computeInvertAndLog_stride(cuda_handles, psiM_fp_, invM_fp_, n, lda, log_values);
    for (int iw = 0; iw < a_mats.size(); ++iw)
    {
      DualMatrix<VALUE_FP> data_ref_matrix;
      data_ref_matrix.attachReference(invM_fp_.data() + nsqr * iw, n, lda);
      // We can't use operator= with different lda, ldb which can happen so we use this assignment which is over the
      // smaller of the two's dimensions
      inv_a_mats[iw].get().assignUpperLeft(data_ref_matrix);
      cudaErrorCheck(cudaMemcpyAsync(inv_a_mats[iw].get().device_data(), inv_a_mats[iw].get().data(),
                                     inv_a_mats[iw].get().size() * sizeof(TMAT), cudaMemcpyHostToDevice,
                                     cuda_handles.hstream),
                     "cudaMemcpyAsync of inv_a_mat to device failed!");
    }
  }

  /** Batched inversion and calculation of log determinants.
   *  When TMAT is full precision we can use the a_mat and inv_mat directly
   *  Side effect of this is after this call the device copy of  a_mats contains
   *  the LU factorization matrix.
   */
  template<typename TMAT>
  inline std::enable_if_t<std::is_same<VALUE_FP, TMAT>::value> mw_invertTranspose(
      CUDALinearAlgebraHandles& cuda_handles,
      const RefVector<const DualMatrix<TMAT>>& a_mats,
      const RefVector<DualMatrix<TMAT>>& inv_a_mats,
      DualVector<LogValue>& log_values)
  {
    assert(log_values.size() == a_mats.size());
    const int n = a_mats[0].get().rows();
    mw_computeInvertAndLog(cuda_handles, a_mats, inv_a_mats, n, log_values);
  }
};

} // namespace qmcplusplus

#endif //QMCPLUSPLUS_DIRAC_MATRIX_COMPUTE_CUDA_H
