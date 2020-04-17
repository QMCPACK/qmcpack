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

#include "simd/allocator.hpp"
#include "OpenMP/OMPallocator.hpp"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "QMCWaveFunctions/Fermion/DiracMatrix.h"
#include "Platforms/OpenMP/ompBLAS.hpp"

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
  template<typename TREAL, typename OMPALLOC>
  inline void invert_transpose(const Matrix<T>& logdetT, Matrix<T, OMPALLOC>& Ainv, std::complex<TREAL>& LogValue)
  {
    Matrix<T> Ainv_host_view(Ainv.data(), Ainv.rows(), Ainv.cols());
    detEng.invert_transpose(logdetT, Ainv_host_view, LogValue);
    auto* Ainv_ptr  = Ainv.data();
    PRAGMA_OFFLOAD("omp target update to(Ainv_ptr[:Ainv.rows()*Ainv.cols()])")
  }

  template<typename VVT, typename RATIOT, typename OMPALLOC>
  inline void updateRow(Matrix<T, OMPALLOC>& Ainv, int rowchanged, const VVT& psiV, RATIOT c_ratio_in)
  {
    // update the inverse matrix
    constexpr T cone(1);
    constexpr T czero(0);
    const int norb = Ainv.rows();
    temp.resize(norb);
    rcopy.resize(norb);
    // invoke the Fahy's variant of Sherman-Morrison update.
    int dummy_handle = 0;
    int success =0;
    success = ompBLAS::gemv(dummy_handle, 'T', norb, norb, cone, Ainv.data(), norb, psiV.data(), 1, czero, temp.data(), 1);
    auto* Ainv_ptr  = Ainv.data();
    auto* temp_ptr  = temp.data();
    auto* rcopy_ptr = rcopy.data();
    PRAGMA_OFFLOAD("omp target map(to: Ainv_ptr[:norb*norb]) map(tofrom: temp_ptr[:norb]) map(from: rcopy_ptr[:norb])")
    {
      temp_ptr[rowchanged] -= cone;
      PRAGMA_OFFLOAD("omp parallel for simd")
      for(int i = 0; i < norb; i++)
        rcopy_ptr[i] = Ainv_ptr[rowchanged * norb + i];
    }
    success = ompBLAS::ger(dummy_handle, norb, norb, static_cast<T>(-RATIOT(1)/c_ratio_in), rcopy.data(), 1, temp.data(), 1, Ainv.data(), norb);
    PRAGMA_OFFLOAD("omp target update from(Ainv_ptr[:Ainv.rows()*Ainv.cols()])")
  }

};
} // namespace qmcplusplus

#endif // QMCPLUSPLUS_MATRIX_UPDATE_H
