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
#include "type_traits/complex_help.hpp"
#include "type_traits/template_types.hpp"
#include "Concurrency/OpenMP.h"
#include "CPU/SIMD/simd.hpp"
#include "ResourceCollection.h"

namespace qmcplusplus
{
/** class to compute matrix inversion and the log value of determinant
 *  of a batch of DiracMatrixes.
 *
 *  @tparam VALUE_FP the datatype used in the actual computation of the matrix
 *  
 *  There is one per crowd not one per MatrixUpdateEngine.
 *  this puts ownership of the scratch resources in a sensible place.
 *  
 *  Currently this is CPU only but its external API is somewhat written to
 *  enforce the passing Dual data objects as arguments.  Except for the single
 *  particle API log_value which is not Dual type but had better have an address in a OMPtarget
 *  mapped region if target is used with it. This makes this API incompatible to
 *  that used by MatrixDelayedUpdateCuda and DiracMatrixComputeCUDA.
 */
template<typename VALUE_FP>
class DiracMatrixComputeOMPTarget : public Resource
{
public:
  using FullPrecReal = RealAlias<VALUE_FP>;
  using LogValue     = std::complex<FullPrecReal>;

  // This class only works with OMPallocator so explicitly call OffloadAllocator what it
  // is and not DUAL
  template<typename T>
  using OffloadPinnedAllocator = OMPallocator<T, PinnedAlignedAllocator<T>>;
  template<typename T>
  using OffloadPinnedMatrix = Matrix<T, OffloadPinnedAllocator<T>>;
  template<typename T>
  using OffloadPinnedVector = Vector<T, OffloadPinnedAllocator<T>>;

  // maybe you'll want a resource someday, then change here.
  using HandleResource = DummyResource;

private:
  aligned_vector<VALUE_FP> m_work_;
  int lwork_;

  /// Matrices held in memory matrices n^2 * nw  elements
  OffloadPinnedVector<VALUE_FP> psiM_fp_;
  OffloadPinnedVector<VALUE_FP> LU_diags_fp_;
  OffloadPinnedVector<int> pivots_;
  OffloadPinnedVector<int> infos_;

  /** reset internal work space.
   *  My understanding might be off.
   *
   *  it smells that this is so complex.
   */
  inline void reset(OffloadPinnedVector<VALUE_FP>& psi_Ms, const int n, const int lda, const int batch_size)
  {
    const int nw = batch_size;
    pivots_.resize(lda * nw);
    for (int iw = 0; iw < nw; ++iw)
    {
      lwork_ = -1;
      VALUE_FP tmp;
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
  inline void reset(OffloadPinnedMatrix<VALUE_FP>& psi_M, const int n, const int lda)
  {
    pivots_.resize(lda);
    LU_diags_fp_.resize(lda);
    lwork_ = -1;
    VALUE_FP tmp;
    FullPrecReal lw;
    Xgetri(lda, psi_M.data(), lda, pivots_.data(), &tmp, lwork_);
    lw     = std::real(tmp);
    lwork_ = static_cast<int>(lw);
    m_work_.resize(lwork_);
  }

  /** compute the inverse of invMat (in place) and the log value of determinant
   * \tparam TMAT value type of matrix
   * \param[inout] a_mat      the matrix
   * \param[in]    n          actual dimension of square matrix (no guarantee it really has full column rank)
   * \param[in]    lda        leading dimension of Matrix container
   * \param[out]   log_value  log a_mat before inversion
   */
  template<typename TMAT>
  inline void computeInvertAndLog(OffloadPinnedMatrix<TMAT>& a_mat, const int n, const int lda, LogValue& log_value)
  {
    BlasThreadingEnv knob(getNextLevelNumThreads());
    if (lwork_ < lda)
      reset(a_mat, n, lda);
    Xgetrf(n, n, a_mat.data(), lda, pivots_.data());
    for (int i = 0; i < n; i++)
      LU_diags_fp_[i] = a_mat.data()[i * lda + i];
    log_value = {0.0, 0.0};
    computeLogDet(LU_diags_fp_.data(), n, pivots_.data(), log_value);
    Xgetri(n, a_mat.data(), lda, pivots_.data(), m_work_.data(), lwork_);
  }

  template<typename TMAT>
  inline void computeInvertAndLog(OffloadPinnedVector<TMAT>& psi_Ms,
                                  const int n,
                                  const int lda,
                                  OffloadPinnedVector<LogValue>& log_values)
  {
    const int nw = log_values.size();
    BlasThreadingEnv knob(getNextLevelNumThreads());
    if (lwork_ < lda)
      reset(psi_Ms, n, lda, nw);
    pivots_.resize(n * nw);
    LU_diags_fp_.resize(n * nw);
    for (int iw = 0; iw < nw; ++iw)
    {
      VALUE_FP* LU_M = psi_Ms.data() + iw * n * n;
      Xgetrf(n, n, LU_M, lda, pivots_.data() + iw * n);
      for (int i = 0; i < n; i++)
        *(LU_diags_fp_.data() + iw * n + i) = LU_M[i * lda + i];
      LogValue log_value{0.0, 0.0};
      computeLogDet(LU_diags_fp_.data() + iw * n, n, pivots_.data() + iw * n, log_value);
      log_values[iw] = log_value;
      Xgetri(n, LU_M, lda, pivots_.data() + iw * n, m_work_.data(), lwork_);
    }
  }

  /// matrix inversion engine
  DiracMatrix<VALUE_FP> detEng_;

public:
  DiracMatrixComputeOMPTarget() : Resource("DiracMatrixComputeOMPTarget"), lwork_(0) {}

  std::unique_ptr<Resource> makeClone() const override { return std::make_unique<DiracMatrixComputeOMPTarget>(*this); }

  /** compute the inverse of the transpose of matrix A and its determinant value in log
   * when VALUE_FP and TMAT are the same
   * @tparam TMAT matrix value type
   * @tparam TREAL real type
   * \param [in]    resource          compute resource
   * \param [in]    a_mat             matrix to be inverted
   * \param [out]   inv_a_mat         the inverted matrix
   * \param [out]   log_value         breaks compatibility of MatrixUpdateOmpTarget with
   *                                  DiracMatrixComputeCUDA but is fine for OMPTarget        
   */
  template<typename TMAT>
  inline std::enable_if_t<std::is_same<VALUE_FP, TMAT>::value> invert_transpose(HandleResource& resource,
                                                                                const OffloadPinnedMatrix<TMAT>& a_mat,
                                                                                OffloadPinnedMatrix<TMAT>& inv_a_mat,
                                                                                LogValue& log_value)
  {
    const int n   = a_mat.rows();
    const int lda = a_mat.cols();
    const int ldb = inv_a_mat.cols();
    simd::transpose(a_mat.data(), n, lda, inv_a_mat.data(), n, ldb);
    // In this case we just pass the value since
    // that makes sense for a single walker API
    computeInvertAndLog(inv_a_mat, n, ldb, log_value);
  }

  /** compute the inverse of the transpose of matrix A and its determinant value in log
   * when VALUE_FP and TMAT are the different
   * @tparam TMAT matrix value type
   * @tparam TREAL real type
   */
  template<typename TMAT>
  inline std::enable_if_t<!std::is_same<VALUE_FP, TMAT>::value> invert_transpose(HandleResource& resource,
                                                                                 const OffloadPinnedMatrix<TMAT>& a_mat,
                                                                                 OffloadPinnedMatrix<TMAT>& inv_a_mat,
                                                                                 LogValue& log_value)
  {
    const int n   = a_mat.rows();
    const int lda = a_mat.cols();
    const int ldb = inv_a_mat.cols();

    psiM_fp_.resize(n * lda);
    simd::transpose(a_mat.data(), n, lda, psiM_fp_.data(), n, lda);
    OffloadPinnedMatrix<VALUE_FP> psiM_fp_view(psiM_fp_, psiM_fp_.data(), n, lda);
    computeInvertAndLog(psiM_fp_view, n, lda, log_value);

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
  template<typename TMAT>
  inline void mw_invertTranspose(HandleResource& resource,
                                 const RefVector<const OffloadPinnedMatrix<TMAT>>& a_mats,
                                 const RefVector<OffloadPinnedMatrix<TMAT>>& inv_a_mats,
                                 OffloadPinnedVector<LogValue>& log_values)
  {
    for (int iw = 0; iw < a_mats.size(); iw++)
    {
      auto& Ainv = inv_a_mats[iw].get();
      detEng_.invert_transpose(a_mats[iw].get(), Ainv, log_values[iw]);
      Ainv.updateTo();
    }

    /* FIXME
    const int nw     = a_mats.size();
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
    */
  }
};
} // namespace qmcplusplus

#endif // QMCPLUSPLUS_DIRAC_MATRIX_COMPUTE_OMPTARGET_H
