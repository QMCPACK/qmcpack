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

#include "Configuration.h"
#include "CPU/Blasf.h"
#include "CPU/BlasThreadingEnv.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OMPTarget/AccelBLAS_OMPTarget.hpp"
#include "OMPTarget/OffloadAlignedAllocators.hpp"
#include "DiracMatrix.h"
#include "type_traits/complex_help.hpp"
#include "type_traits/template_types.hpp"
#include "Concurrency/OpenMP.h"
#include "CPU/SIMD/algorithm.hpp"
#include "DiracMatrixInverter.hpp"

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
 *  that used by MatrixDelayedUpdateCuda and DiracMatrixInverterCUDA.
 */
template<typename VALUE_FP, typename VALUE = VALUE_FP>
class DiracMatrixInverterOMPTarget : public DiracMatrixInverter<VALUE_FP, VALUE>
{
public:
  using FullPrecReal = RealAlias<VALUE_FP>;
  using LogValue     = std::complex<FullPrecReal>;

  // This class only works with OMPallocator so explicitly call OffloadAllocator what it
  // is and not DUAL
  template<typename T>
  using OffloadPinnedMatrix = Matrix<T, OffloadPinnedAllocator<T>>;
  template<typename T>
  using OffloadPinnedVector = Vector<T, OffloadPinnedAllocator<T>>;

private:
  /// matrix inversion engine
  DiracMatrix<VALUE_FP> detEng_;

public:
  DiracMatrixInverterOMPTarget() : DiracMatrixInverter<VALUE_FP, VALUE>("DiracMatrixInverterOMPTarget") {}

  std::unique_ptr<Resource> makeClone() const override { return std::make_unique<DiracMatrixInverterOMPTarget>(*this); }

  /** compute the inverse of the transpose of matrix A and its determinant value in log
   * when VALUE_FP and TMAT are the same
   * @tparam TMAT matrix value type
   * @tparam TREAL real type
   * \param [in]    resource          compute resource
   * \param [in]    a_mat             matrix to be inverted
   * \param [out]   inv_a_mat         the inverted matrix
   * \param [out]   log_value         breaks compatibility of MatrixUpdateOmpTarget with
   *                                  DiracMatrixInverterCUDA but is fine for OMPTarget        
   */
  template<typename TMAT>
  inline void invert_transpose(const OffloadPinnedMatrix<TMAT>& a_mat,
                               OffloadPinnedMatrix<TMAT>& inv_a_mat,
                               LogValue& log_value)
  {
    detEng_.invert_transpose(a_mat, inv_a_mat, log_value);
    inv_a_mat.updateTo();
  }

  /** This covers both mixed and Full precision case.
   *  
   *  \todo measure if using the a_mats without a copy to contiguous vector is better.
   */
  template<typename TMAT>
  inline void mw_invertTranspose(const RefVector<const OffloadPinnedMatrix<TMAT>>& a_mats,
                                 const RefVector<OffloadPinnedMatrix<TMAT>>& inv_a_mats,
                                 OffloadPinnedVector<LogValue>& log_values)
  {
    for (int iw = 0; iw < a_mats.size(); iw++)
    {
      auto& Ainv = inv_a_mats[iw].get();
      detEng_.invert_transpose(a_mats[iw].get(), Ainv, log_values[iw]);
      Ainv.updateTo();
    }
  }

  void mw_invert_transpose(compute::QueueBase& queue_ignored,
                           const RefVector<const OffloadPinnedMatrix<VALUE>>& a_mats,
                           const RefVector<OffloadPinnedMatrix<VALUE>>& inv_a_mats,
                           OffloadPinnedVector<LogValue>& log_values) override
  {
    mw_invertTranspose(a_mats, inv_a_mats, log_values);
  }
};

extern template class DiracMatrixInverterOMPTarget<double, float>;
extern template class DiracMatrixInverterOMPTarget<double, double>;
extern template class DiracMatrixInverterOMPTarget<std::complex<double>, std::complex<float>>;
extern template class DiracMatrixInverterOMPTarget<std::complex<double>, std::complex<double>>;
} // namespace qmcplusplus

#endif // QMCPLUSPLUS_DIRAC_MATRIX_COMPUTE_OMPTARGET_H
