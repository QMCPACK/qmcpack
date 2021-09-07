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
#include "Message/OpenMP.h"
#include "CPU/SIMD/simd.hpp"
#include "ResourceCollection.h"

namespace qmcplusplus
{
/** class defining a compute and memory resource to compute matrix inversion and 
 *  the log value of determinant of a batch of DiracMatrixes.
 *  Multiplicty is one per crowd not one per MatrixUpdateEngine.
 *  Matching the multiplicity of the call and resource requirement.
 *  At least for CUDA this is useful since memory movement
 *  and computation can be handled on the same stream which is idea for maximum async.
 *
 *  @tparam VALUE_FP the datatype used in the actual computation of matrix inversion
 *              currently we only support std::complex<double>
 */
template<typename VALUE_FP>
class DiracMatrixComputeCUDA : public Resource
{
  using FullPrecReal = RealAlias<VALUE_FP>;

  template<typename T>
  using DualMatrix = Matrix<T, PinnedDualAllocator<T>>;

  template<typename T>
  using DualVector = Vector<T, PinnedDualAllocator<T>>;

  cudaStream_t hstream_;

  // Contiguous memory for fp precision Matrices for each walker, n^2 * nw elements
  DualVector<VALUE_FP> psiM_fp_;
  DualVector<VALUE_FP> invM_fp_;

  // working vectors
  DualVector<VALUE_FP> LU_diags_fp_;
  DualVector<int> pivots_;
  DualVector<int> infos_;

  //DualMatrix<T_FP> temp_mat_;

  // For device pointers to matrices
  Vector<char, OffloadPinnedAllocator<char>> psiM_ptrs_;
  Vector<char, OffloadPinnedAllocator<char>> invM_ptrs_;

  // cuBLAS geam wants these.
  VALUE_FP host_one{1.0};
  VALUE_FP host_zero{0.0};

  /** Calculates the actual inv and log determinant on accelerator
   *
   *  \param[in]      h_cublas    cublas handle, hstream handle is retrieved from it.			
   *  \param[in,out]  a_mats      dual A matrices, they will be transposed on the device side as a side effect.
   *  \param[out]     inv_a_mats  dual invM matrices
   *  \param[in]      n           matrices rank.								
   *  \param[in]      lda         leading dimension of each matrix					
   *  \param[out]     log_values  log determinant value for each matrix, batch_size = log_values.size()
   *
   *  This is the prefered interface for calling computeInvertAndLog when the VALUE_FP and TMAT are equal.
   *  On Volta so far little seems to be achieved by having the mats continuous.
   */
  template<typename TMAT>
  inline void mw_computeInvertAndLog(CUDALinearAlgebraHandles& cuda_handles,
                                     RefVector<DualMatrix<TMAT>>& a_mats,
                                     RefVector<DualMatrix<TMAT>>& inv_a_mats,
                                     const int n,
                                     const int lda,
                                     DualVector<std::complex<FullPrecReal>>& log_values)
  {
    int nw = a_mats.size();
    assert(a_mats.size() == inv_a_mats.size());

    psiM_ptrs_.resize(sizeof(TMAT*) * nw);
    invM_ptrs_.resize(sizeof(TMAT*) * nw);
    //temp_mat_.resize(n,lda);
    int ldinv = inv_a_mats[0].get().cols();
    Vector<const VALUE_FP*> M_ptr_buffer(reinterpret_cast<const TMAT**>(psiM_ptrs_.data()), nw);
    Vector<const VALUE_FP*> invM_ptr_buffer(reinterpret_cast<const TMAT**>(invM_ptrs_.data()), nw);
    cudaStream_t hstream = cuda_handles.hstream;
    cublasHandle_t h_cublas = cuda_handles.h_cublas;

    for (int iw = 0; iw < nw; ++iw)
    {
      M_ptr_buffer[iw]    = a_mats[iw].get().device_data();
      invM_ptr_buffer[iw] = inv_a_mats[iw].get().device_data();
      // Since inv_a_mat can have a different leading dimension from a_mat first we remap copy on the host
      simd::remapCopy(n, n, a_mats[iw].get().data(), lda, inv_a_mats[iw].get().data(), ldinv);
      // Then copy a_mat in inv_a_mats to the device
      cudaErrorCheck(cudaMemcpyAsync((void*)(inv_a_mats[iw].get().device_data()), (void*)(inv_a_mats[iw].get().data()),
                                     inv_a_mats[iw].get().size() * sizeof(TMAT), cudaMemcpyHostToDevice, hstream),
                     "cudaMemcpyAsync failed copying DiracMatrixBatch::psiM to device");
      // On the device Here we transpose to a_mat;
      cublasErrorCheck(cuBLAS::geam(h_cublas, CUBLAS_OP_T, CUBLAS_OP_N, n, n, &host_one,
                                    inv_a_mats[iw].get().device_data(), ldinv, &host_zero,
                                    a_mats[iw].get().device_data(), lda, a_mats[iw].get().device_data(), lda),
                       "cuBLAS::geam failed.");
    }
    pivots_.resize(n * nw);
    infos_.resize(nw);
    LU_diags_fp_.resize(n * nw);
    cudaErrorCheck(cudaMemcpyAsync(psiM_ptrs_.device_data(), psiM_ptrs_.data(), psiM_ptrs_.size(),
                                   cudaMemcpyHostToDevice, hstream),
                   "cudaMemcpyAsync psiM_ptrs_ failed!");
    cudaErrorCheck(cudaMemcpyAsync(invM_ptrs_.device_data(), invM_ptrs_.data(), invM_ptrs_.size(),
                                   cudaMemcpyHostToDevice, hstream),
                   "cudaMemcpyAsync invM_ptrs_ failed!");
    TMAT** psiM_mw_ptr = reinterpret_cast<TMAT**>(psiM_ptrs_.device_data());
    TMAT** invM_mw_ptr = reinterpret_cast<TMAT**>(invM_ptrs_.device_data());
    cuBLAS_LU::computeInverseAndDetLog_batched(h_cublas, hstream, n, lda, psiM_mw_ptr, invM_mw_ptr,
                                               LU_diags_fp_.device_data(), pivots_.device_data(), infos_.data(),
                                               infos_.device_data(), log_values.device_data(), nw);
    for (int iw = 0; iw < nw; ++iw)
    {
      cudaErrorCheck(cudaMemcpyAsync(inv_a_mats[iw].get().data(), inv_a_mats[iw].get().device_data(),
                                     inv_a_mats[iw].get().size() * sizeof(TMAT), cudaMemcpyDeviceToHost, hstream),
                     "cudaMemcpyAsync failed copying DiracMatrixBatch::inv_psiM to host");
    }
    cudaErrorCheck(cudaMemcpyAsync(log_values.data(), log_values.device_data(),
                                   log_values.size() * sizeof(std::complex<FullPrecReal>), cudaMemcpyDeviceToHost,
                                   hstream),
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
   */
  inline void mw_computeInvertAndLog(CUDALinearAlgebraHandles& cuda_handles,
                                     DualVector<VALUE_FP>& psi_Ms,
                                     DualVector<VALUE_FP>& inv_Ms,
                                     const int n,
                                     const int lda,
                                     DualVector<std::complex<FullPrecReal>>& log_values)
  {
    // This is probably dodgy
    int nw = log_values.size();
    psiM_ptrs_.resize(sizeof(VALUE_FP*) * nw);
    invM_ptrs_.resize(sizeof(VALUE_FP*) * nw);
    Vector<const VALUE_FP*> M_ptr_buffer(reinterpret_cast<const VALUE_FP**>(psiM_ptrs_.data()), nw);
    Vector<const VALUE_FP*> invM_ptr_buffer(reinterpret_cast<const VALUE_FP**>(invM_ptrs_.data()), nw);
    for (int iw = 0; iw < nw; ++iw)
    {
      M_ptr_buffer[iw]    = psi_Ms.device_data() + iw * n * lda;
      invM_ptr_buffer[iw] = inv_Ms.device_data() + iw * n * lda;
    }
    pivots_.resize(n * nw);
    infos_.resize(nw);
    LU_diags_fp_.resize(n * nw);

    cudaStream_t hstream = cuda_handles.hstream;
    cublasHandle_t h_cublas = cuda_handles.h_cublas;
    cudaErrorCheck(cudaMemcpyAsync(psi_Ms.device_data(), psi_Ms.data(), psi_Ms.size() * sizeof(VALUE_FP),
                                   cudaMemcpyHostToDevice, hstream),
                   "cudaMemcpyAsync failed copying DiracMatrixBatch::psiM_fp to device");
    cudaErrorCheck(cudaMemcpyAsync(psiM_ptrs_.device_data(), psiM_ptrs_.data(), psiM_ptrs_.size(),
                                   cudaMemcpyHostToDevice, hstream),
                   "cudaMemcpyAsync psiM_ptrs_ failed!");
    cudaErrorCheck(cudaMemcpyAsync(invM_ptrs_.device_data(), invM_ptrs_.data(), invM_ptrs_.size(),
                                   cudaMemcpyHostToDevice, hstream),
                   "cudaMemcpyAsync invM_ptrs_ failed!");
    VALUE_FP** psiM_mw_ptr = reinterpret_cast<VALUE_FP**>(psiM_ptrs_.device_data());
    VALUE_FP** invM_mw_ptr = reinterpret_cast<VALUE_FP**>(invM_ptrs_.device_data());

    cuBLAS_LU::computeInverseAndDetLog_batched(h_cublas, hstream, n, lda, psiM_mw_ptr, invM_mw_ptr,
                                               LU_diags_fp_.device_data(), pivots_.device_data(), infos_.data(),
                                               infos_.device_data(), log_values.device_data(), nw);
#if NDEBUG
    // This is very useful to see whether the data after all kernels and cublas calls are run is wrong on the device or due to copy.
    // cuBLAS_LU::peekinvM_batched(hstream, psiM_mw_ptr, invM_mw_ptr, pivots_.device_data(), infos_.device_data(),
    //                             log_values.device_data(), nw);
#endif
    cudaErrorCheck(cudaMemcpyAsync(invM_ptrs_.data(), invM_ptrs_.device_data(), invM_ptrs_.size(),
                                   cudaMemcpyDeviceToHost, hstream),
                   "cudaMemcpyAsync invM_ptrs_ failed!");
    cudaErrorCheck(cudaMemcpyAsync(inv_Ms.data(), inv_Ms.device_data(), inv_Ms.size() * sizeof(VALUE_FP),
                                   cudaMemcpyDeviceToHost, hstream),
                   "cudaMemcpyAsync failed copying back DiracMatrixBatch::invM_fp from device");
    cudaErrorCheck(cudaMemcpyAsync(log_values.data(), log_values.device_data(),
                                   log_values.size() * sizeof(std::complex<FullPrecReal>), cudaMemcpyDeviceToHost,
                                   hstream),
                   "cudaMemcpyAsync log_values failed!");
    cudaErrorCheck(cudaStreamSynchronize(hstream), "cudaStreamSynchronize failed!");
  }

public:
  DiracMatrixComputeCUDA(cudaStream_t hstream) : Resource("DiracMatrixComputeCUDA"), hstream_(hstream) {}

  DiracMatrixComputeCUDA(const DiracMatrixComputeCUDA& other, cudaStream_t hstream)
      : Resource(other.getName()), hstream_(hstream)
  {}

  Resource* makeClone(cudaStream_t hstream) const { return new DiracMatrixComputeCUDA(*this, hstream); }
  Resource* makeClone() const override { return new DiracMatrixComputeCUDA(*this, this->hstream_); }

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
                        DualVector<std::complex<FullPrecReal>>& log_values)
  {
    const int n        = a_mat.rows();
    const int lda      = a_mat.cols();
    psiM_fp_.resize(n * lda);
    invM_fp_.resize(n * lda);
    std::fill(log_values.begin(), log_values.end(), std::complex<FullPrecReal>{0.0, 0.0});
    // making sure we know the log_values are zero'd on the device.
    cudaErrorCheck(cudaMemcpyAsync(log_values.device_data(), log_values.data(),
                                   log_values.size() * sizeof(std::complex<FullPrecReal>), cudaMemcpyHostToDevice,
                                   cuda_handles.hstream),
                   "cudaMemcpyAsync failed copying DiracMatrixBatch::log_values to device");
    simd::transpose(a_mat.data(), n, lda, psiM_fp_.data(), n, lda);
    cudaErrorCheck(cudaMemcpyAsync(psiM_fp_.device_data(), psiM_fp_.data(), psiM_fp_.size() * sizeof(VALUE_FP),
                                   cudaMemcpyHostToDevice, cuda_handles.hstream),
                   "cudaMemcpyAsync failed copying DiracMatrixBatch::psiM_fp to device");
    mw_computeInvertAndLog(cuda_handles, psiM_fp_, invM_fp_, n, lda, log_values);
    DualMatrix<VALUE_FP> data_ref_matrix;

    data_ref_matrix.attachReference(invM_fp_.data(), n, n);

    // We can't use operator= with different lda, ldb which can happen so we use this assignment which is over the
    // smaller of the two's dimensions
    inv_a_mat.assignUpperRight(data_ref_matrix);
    cudaErrorCheck(cudaMemcpyAsync(inv_a_mat.device_data(), inv_a_mat.data(), inv_a_mat.size() * sizeof(TMAT),
                                   cudaMemcpyHostToDevice, cuda_handles.hstream),
                   "cudaMemcpyAsync of inv_a_mat to device failed!");
  }

  /** Mixed precision specialization
   *  When TMAT is not full precision we need to still do the inversion and log
   *  at full precision. This is not yet optimized to transpose on the GPU
   */
  template<typename TMAT>
  inline std::enable_if_t<!std::is_same<VALUE_FP, TMAT>::value> mw_invertTranspose(
      CUDALinearAlgebraHandles& cuda_handles,
      RefVector<DualMatrix<TMAT>>& a_mats,
      RefVector<DualMatrix<TMAT>>& inv_a_mats,
      DualVector<std::complex<FullPrecReal>>& log_values,
      const std::vector<bool>& compute_mask)
  {
    assert(log_values.size() == a_mats.size());
    int nw             = a_mats.size();
    const int n        = a_mats[0].get().rows();
    const int lda      = a_mats[0].get().cols();
    size_t nsqr        = n * n;
    psiM_fp_.resize(n * lda * nw);
    invM_fp_.resize(n * lda * nw);
    std::fill(log_values.begin(), log_values.end(), std::complex<FullPrecReal>{0.0, 0.0});
    // making sure we know the log_values are zero'd on the device.
    cudaErrorCheck(cudaMemcpyAsync(log_values.device_data(), log_values.data(),
                                   log_values.size() * sizeof(std::complex<FullPrecReal>), cudaMemcpyHostToDevice,
                                   cuda_handles.hstream),
                   "cudaMemcpyAsync failed copying DiracMatrixBatch::log_values to device");
    for (int iw = 0; iw < nw; ++iw)
      simd::transpose(a_mats[iw].get().data(), n, a_mats[iw].get().cols(), psiM_fp_.data() + nsqr * iw, n, lda);
    mw_computeInvertAndLog(cuda_handles.h_cublas, psiM_fp_, invM_fp_, n, lda, log_values);
    for (int iw = 0; iw < a_mats.size(); ++iw)
    {
      DualMatrix<VALUE_FP> data_ref_matrix;
      data_ref_matrix.attachReference(invM_fp_.data() + nsqr * iw, n, lda);
      // We can't use operator= with different lda, ldb which can happen so we use this assignment which is over the
      // smaller of the two's dimensions
      inv_a_mats[iw].get().assignUpperRight(data_ref_matrix);
      cudaErrorCheck(cudaMemcpyAsync(inv_a_mats[iw].get().device_data(), inv_a_mats[iw].get().data(),
                                     inv_a_mats[iw].get().size() * sizeof(TMAT), cudaMemcpyHostToDevice,
                                     cuda_handles.hstream),
                     "cudaMemcpyAsync of inv_a_mat to device failed!");
    }
  }

  /** Batched inversion and calculation of log determinants.
   *  When TMAT is full precision we can use the a_mat and inv_mat directly
   *  Side effect of this is after this call a_mats contains the LU factorization
   *  matrix.
   */
  template<typename TMAT>
  inline std::enable_if_t<std::is_same<VALUE_FP, TMAT>::value> mw_invertTranspose(
      CUDALinearAlgebraHandles& cuda_handles,
      RefVector<DualMatrix<TMAT>>& a_mats,
      RefVector<DualMatrix<TMAT>>& inv_a_mats,
      DualVector<std::complex<FullPrecReal>>& log_values,
      const std::vector<bool>& compute_mask)
  {
    assert(log_values.size() == a_mats.size());
    const int n   = a_mats[0].get().rows();
    const int lda = a_mats[0].get().cols();
    mw_computeInvertAndLog(cuda_handles, a_mats, inv_a_mats, n, lda, log_values);
  }
};

} // namespace qmcplusplus

#endif //QMCPLUSPLUS_DIRAC_MATRIX_COMPUTE_CUDA_H
