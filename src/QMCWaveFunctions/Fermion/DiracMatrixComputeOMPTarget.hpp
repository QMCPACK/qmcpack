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

#ifndef QMCPLUSPLUS_DIRAC_MATRIX_COMPUTE_OMPTARGET_H
#define QMCPLUSPLUS_DIRAC_MATRIX_COMPUTE_OMPTARGET_H

#include "CPU/Blasf.h"
#include "CPU/BlasThreadingEnv.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OMPTarget/OMPallocator.hpp"
#include "Platforms/PinnedAllocator.h"
#include "DiracMatrix.h"
#include "type_traits/scalar_traits.h"
#include "type_traits/template_types.hpp"
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
 *  this puts ownership of the scratch resources in a sensible place.
 *  Even for CPU this is better.
 */
template<typename T_FP>
class DiracMatrixComputeOMPTarget : public Resource
{
  // Why not just use QMCTraits::FullPrecRealType?
  using FullPrecReal = typename scalar_traits<T_FP>::real_type;

  template<typename T>
  using OffloadPinnedAllocator = OMPallocator<T, PinnedAlignedAllocator<T>>;

  template<typename T>
  using OffloadPinnedMatrix = Matrix<T, OffloadPinnedAllocator<T>>;

  template<typename T>
  using OffloadPinnedVector = Vector<T, OffloadPinnedAllocator<T>>;

  aligned_vector<T_FP> m_work_;
  int lwork_;
  //Unlike DiracMatrix.h these are contiguous packed representations of the Matrices
  // Contiguous Matrices for each walker, n^2 * nw  elements
  OffloadPinnedVector<T_FP> psiM_fp_;
  OffloadPinnedVector<T_FP> logdets_fp_;
  OffloadPinnedVector<T_FP> LU_diags_fp_;
  OffloadPinnedVector<int> pivots_;
  OffloadPinnedVector<int> infos_;

  Vector<char, OffloadPinnedAllocator<char>> psiM_ptrs_;
  /** reset internal work space.
   *  My understanding might be off.
   *
   *  it smells that this is so complex.
   */
  inline void reset(OffloadPinnedVector<T_FP>& psi_Ms, const int n, const int lda, const int batch_size)
  {
    int nw = batch_size;
    pivots_.resize(lda * nw);
    for (int iw = 0; iw < nw; ++iw)
    {
      lwork_ = -1;
      T_FP tmp;
      FullPrecReal lw;
      auto psi_M_ptr = psi_Ms.data() + iw * n * n;
      Xgetri(lda, psi_M_ptr, lda, pivots_.data() + iw * n, &tmp, lwork_);
      convert(tmp, lw);
      lwork_ = static_cast<int>(lw);
      m_work_.resize(lwork_);
    }
  }

  /** reset internal work space for single walker case
   *  My understanding might be off.
   *
   *  it smells that this is so complex.
   */
  template<typename TREAL>
  inline void reset(OffloadPinnedMatrix<TREAL>& psi_M, const int n, const int lda)
  {
    pivots_.resize(lda);
    LU_diags_fp_.resize(lda);
    lwork_ = -1;
    T_FP tmp;
    FullPrecReal lw;
    Xgetri(lda, psi_M.data(), lda, pivots_.data(), &tmp, lwork_);
    convert(tmp, lw);
    lwork_ = static_cast<int>(lw);
    m_work_.resize(lwork_);
  }

  /** compute the inverse of invMat (in place) and the log value of determinant
   * @tparam TREAL real type
   * @param n invMat is n x n matrix
   * @param lda the first dimension of invMat
   * @param LogDet log determinant value of invMat before inversion
   */
  template<typename TVALUE>
  inline void computeInvertAndLog(OffloadPinnedMatrix<TVALUE>& invMat,
                                  const int n,
                                  const int lda,
                                  std::complex<FullPrecReal>& log_value)
  {
    BlasThreadingEnv knob(getNextLevelNumThreads());
    if (lwork_ < lda)
      reset(invMat, n, lda);
    Xgetrf(n, n, invMat.data(), lda, pivots_.data());
    for (int i = 0; i < n; i++)
      LU_diags_fp_[i] = invMat.data()[i * lda + i];
    log_value = {0.0, 0.0};
    computeLogDet(LU_diags_fp_.data(), n, pivots_.data(), log_value);
    Xgetri(n, invMat.data(), lda, pivots_.data(), m_work_.data(), lwork_);
  }


  template<typename TVALUE>
  inline void computeInvertAndLog(OffloadPinnedVector<TVALUE>& psi_Ms,
                                  const int n,
                                  const int lda,
                                  OffloadPinnedVector<std::complex<FullPrecReal>>& log_values)
  {
    int nw = log_values.size();
    BlasThreadingEnv knob(getNextLevelNumThreads());
    if (lwork_ < lda)
      reset(psi_Ms, n, lda, nw);
    pivots_.resize(n * nw);
    LU_diags_fp_.resize(n * nw);
    for (int iw = 0; iw < nw; ++iw)
    {
      T_FP* LU_M = psi_Ms.data() + iw * n * n;
      Xgetrf(n, n, LU_M, lda, pivots_.data() + iw * n);
      for (int i = 0; i < n; i++)
        *(LU_diags_fp_.data() + iw * n + i) = LU_M[i * lda + i];
      std::complex<FullPrecReal> log_value{0.0, 0.0};
      computeLogDet(LU_diags_fp_.data() + iw * n, n, pivots_.data() + iw * n, log_value);
      log_values[iw] = log_value;
      Xgetri(n, LU_M, lda, pivots_.data() + iw * n, m_work_.data(), lwork_);
    }
  }

public:
  DiracMatrixComputeOMPTarget() : Resource("DiracMatrixComputeOMPTarget"), lwork_(0) {}

  Resource* makeClone() const override { return new DiracMatrixComputeOMPTarget(*this); }

  /** compute the inverse of the transpose of matrix A and its determinant value in log
   * when T_FP and TMAT are the same
   * @tparam TMAT matrix value type
   * @tparam TREAL real type
   */
  template<typename TMAT>
  inline std::enable_if_t<std::is_same<T_FP, TMAT>::value> invert_transpose(
      Resource& resource,
      const OffloadPinnedMatrix<TMAT>& a_mat,
      OffloadPinnedMatrix<TMAT>& inv_a_mat,
      OffloadPinnedVector<std::complex<FullPrecReal>>& log_values)
  {
    const int n   = a_mat.rows();
    const int lda = a_mat.cols();
    const int ldb = inv_a_mat.cols();
    simd::transpose(a_mat.data(), n, lda, inv_a_mat.data(), n, ldb);
    computeInvertAndLog(inv_a_mat, n, ldb, log_values[0]);
  }

  /** compute the inverse of the transpose of matrix A and its determinant value in log
   * when T_FP and TMAT are the different
   * @tparam TMAT matrix value type
   * @tparam TREAL real type
   */
  template<typename TMAT>
  inline std::enable_if_t<!std::is_same<T_FP, TMAT>::value> invert_transpose(
      Resource& resource,
      const OffloadPinnedMatrix<TMAT>& a_mat,
      OffloadPinnedMatrix<TMAT>& inv_a_mat,
      OffloadPinnedVector<std::complex<FullPrecReal>>& log_values)
  {
    const int n   = a_mat.rows();
    const int lda = a_mat.cols();
    const int ldb = inv_a_mat.cols();

    psiM_fp_.resize(n * lda);
    simd::transpose(a_mat.data(), n, lda, psiM_fp_.data(), n, lda);
    OffloadPinnedMatrix<T_FP> psiM_fp_view(psiM_fp_, psiM_fp_.data() , n, lda);
    computeInvertAndLog(psiM_fp_view, n, lda, log_values[0]);

    //Matrix<TMAT> data_ref_matrix;
    //maybe n, lda
    //data_ref_matrix.attachReference(psiM_fp_.data(), n, n);
    //Because inv_a_mat is "aligned" this is unlikely to work.
    simd::remapCopy(n, n, psiM_fp_.data(), lda, inv_a_mat.data(), ldb);
  }

  /** This covers both mixed and Full precision case.
   *  
   *  \todo measure if using the a_mats without a copy to contiguous vector is better.
   */
  template<typename TMAT, typename TREAL>
  inline void mw_invertTranspose(Resource& resource,
                                 RefVector<OffloadPinnedMatrix<TMAT>>& a_mats,
                                 RefVector<OffloadPinnedMatrix<TMAT>>& inv_a_mats,
                                 OffloadPinnedVector<std::complex<TREAL>>& log_values,
                                 const std::vector<bool>& recompute)
  {
    int nw           = a_mats.size();
    const size_t n   = a_mats[0].get().rows();
    const size_t lda = a_mats[0].get().cols();
    const size_t ldb = inv_a_mats[0].get().cols();

    size_t nsqr{n * n};
    psiM_fp_.resize(n * lda * nw);
    for (int iw = 0; iw < nw; ++iw)
      simd::transpose(a_mats[iw].get().data(), n, lda, psiM_fp_.data() + nsqr * iw, n, lda);

    computeInvertAndLog(psiM_fp_, n, lda, log_values);
    for (int iw = 0; iw < nw; ++iw)
    {
      simd::remapCopy(n, n, psiM_fp_.data() + nsqr * iw, lda, inv_a_mats[iw].get().data(), ldb);
    }
  }
};
} // namespace qmcplusplus

// template<typename TREAL>
// inline void computeInvertAndLog(T_FP* invMat, const int n, const int lda, std::complex<TREAL>& log_value);


#endif // QMCPLUSPLUS_DIRAC_MATRIX_COMPUTE_OMPTARGET_H
