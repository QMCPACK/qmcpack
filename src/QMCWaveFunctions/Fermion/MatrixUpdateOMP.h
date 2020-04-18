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
#include "Platforms/PinnedAllocator.h"
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
  template<typename DT>
  using OffloadPinnedAllocator = OMPallocator<DT, PinnedAlignedAllocator<DT>>;
  using OffloadValueVector_t = Vector<T, OffloadAllocator<T>>;
  using OffloadPinnedValueVector_t = Vector<T, OffloadPinnedAllocator<T>>;

  /// matrix inversion engine
  DiracMatrix<T_FP> detEng;
  /// scratch space for rank-1 update
  OffloadValueVector_t temp;
  // scratch space for keeping one row of Ainv
  OffloadValueVector_t rcopy;
  // pointer buffer
  Vector<T*, OffloadPinnedAllocator<T*>> ptr_buffer;

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
    auto* psiV_ptr = psiV.data();
    auto* Ainv_ptr  = Ainv.data();
    auto* temp_ptr  = temp.data();
    auto* rcopy_ptr = rcopy.data();
    PRAGMA_OFFLOAD("omp target data map(always, to: psiV_ptr[:norb]) map(always, from: Ainv_ptr[:Ainv.rows()*Ainv.cols()]) use_device_ptr(psiV_ptr, Ainv_ptr, temp_ptr, rcopy_ptr)")
    {
      success = ompBLAS::gemv(dummy_handle, 'T', norb, norb, cone, Ainv_ptr, norb, psiV_ptr, 1, czero, temp_ptr, 1);
      PRAGMA_OFFLOAD("omp target is_device_ptr(Ainv_ptr, temp_ptr, rcopy_ptr)")
      {
        temp_ptr[rowchanged] -= cone;
        PRAGMA_OFFLOAD("omp parallel for simd")
        for(int i = 0; i < norb; i++)
          rcopy_ptr[i] = Ainv_ptr[rowchanged * norb + i];
      }
      success = ompBLAS::ger(dummy_handle, norb, norb, static_cast<T>(-RATIOT(1)/c_ratio_in), rcopy_ptr, 1, temp_ptr, 1, Ainv_ptr, norb);
    }
  }

  template<typename DEVPTRV, typename RATIOT, typename ALLOC>
  inline void mw_updateRow(const DEVPTRV& Ainv_list, int norb, int rowchanged, const DEVPTRV& psi_v_list, const Vector<RATIOT, ALLOC>& c_ratio_v)
  {
    ptr_buffer.resize(Ainv_list.size() + psi_v_list.size());
    size_t nw = Ainv_list.size();
    for(int i = 0; i < nw; i++)
    {
      ptr_buffer[i] = Ainv_list[i];
      ptr_buffer[i + nw] = psi_v_list[i];
    }

    // update the inverse matrix
    constexpr T cone(1);
    constexpr T czero(0);
    temp.resize(norb * nw);
    rcopy.resize(norb * nw);
    auto* temp_ptr  = temp.data();
    auto* rcopy_ptr = rcopy.data();

    int dummy_handle = 0;
    int success =0;
    auto* ptr_buffer_ptr = ptr_buffer.data();
    PRAGMA_OFFLOAD("omp target data map(always, to: ptr_buffer_ptr[:ptr_buffer.size()]) use_device_ptr(temp_ptr, rcopy_ptr)")
    {
      for(int iw = 0; iw < nw; iw++)
      {
        auto* Ainv_ptr = ptr_buffer[iw];
        auto* psiV_ptr = ptr_buffer[iw + nw];
        auto c_ratio_in = c_ratio_v[iw];
        // invoke the Fahy's variant of Sherman-Morrison update.
        success = ompBLAS::gemv(dummy_handle, 'T', norb, norb, cone, Ainv_ptr, norb, psiV_ptr, 1, czero, temp_ptr + norb * iw, 1);
        PRAGMA_OFFLOAD("omp target is_device_ptr(Ainv_ptr, temp_ptr, rcopy_ptr)")
        {
          temp_ptr[rowchanged + norb * iw] -= cone;
          PRAGMA_OFFLOAD("omp parallel for simd")
          for(int i = 0; i < norb; i++)
            rcopy_ptr[i + norb * iw] = Ainv_ptr[rowchanged * norb + i];
        }
        success = ompBLAS::ger(dummy_handle, norb, norb, static_cast<T>(-RATIOT(1)/c_ratio_in), rcopy_ptr + norb * iw, 1, temp_ptr + norb * iw, 1, Ainv_ptr, norb);
      }
    }
  }

};
} // namespace qmcplusplus

#endif // QMCPLUSPLUS_MATRIX_UPDATE_H
