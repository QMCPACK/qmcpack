//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_MATRIX_UPDATE_OMP_H
#define QMCPLUSPLUS_MATRIX_UPDATE_OMP_H

#include "OpenMP/OMPallocator.hpp"
#include "simd/allocator.hpp"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "QMCWaveFunctions/Fermion/DiracMatrix.h"

namespace qmcplusplus
{

/** implements dirac matrix update using OpenMP
 * @tparam T base precision for most computation
 * @tparam T_FP high precision for matrix inversion, T_FP >= T
 */
template<typename T, typename T_FP>
class MatrixUpdateOMP
{
  template<typename DT>
  using OffloadAllocator = OMPallocator<DT, aligned_allocator<DT>>;

  /// matrix inversion engine
  DiracMatrix<T_FP> detEng;
  /// scratch space for rank-1 update
  Vector<T, OffloadAllocator<T>> temp;
  // scratch space for keeping one row of Ainv
  Vector<T, OffloadAllocator<T>> rcopy;

public:

  /** compute the inverse of the transpose of matrix A
   * @param logdetT orbital value matrix
   * @param Ainv inverse matrix
   */
  template<typename TREAL>
  inline void invert_transpose(const Matrix<T>& logdetT, Matrix<T>& Ainv, std::complex<TREAL>& LogValue)
  {
    detEng.invert_transpose(logdetT, Ainv, LogValue);
  }

  template<typename VVT, typename RATIOT>
  inline void updateRow(Matrix<T>& Ainv, int rowchanged, const VVT& psiV, RATIOT c_ratio_in)
  {
    // update the inverse matrix
    constexpr T cone(1);
    constexpr T czero(0);
    const int norb = Ainv.rows();
    temp.resize(norb);
    rcopy.resize(norb);
    // invoke the Fahy's variant of Sherman-Morrison update.
    BLAS::gemv('T', norb, norb, cone, Ainv.data(), norb, psiV.data(), 1, czero, temp.data(), 1);
    temp[rowchanged] -= cone;
    std::copy_n(Ainv[rowchanged], norb, rcopy.data());
    BLAS::ger(norb, norb, static_cast<T>(-RATIOT(1)/c_ratio_in), rcopy.data(), 1, temp.data(), 1, Ainv.data(), norb);
  }

};
} // namespace qmcplusplus

#endif // QMCPLUSPLUS_MATRIX_UPDATE_H
