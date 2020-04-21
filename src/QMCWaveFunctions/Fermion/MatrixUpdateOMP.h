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
  // scratch space for keeping ratio inverse
  OffloadValueVector_t c_ratio_inv;
  // constant array value T(1)
  OffloadValueVector_t cone_vec;
  // constant array value T(0)
  OffloadValueVector_t czero_vec;
  // pointer buffer
  Vector<T*, OffloadPinnedAllocator<T*>> ptr_buffer;

  void resize_fill_constant_arrays(size_t n)
  {
    if (cone_vec.size() < n)
    {
      cone_vec.resize(n);
      czero_vec.resize(n);
      std::fill_n(cone_vec.data(), n, T(1));
      std::fill_n(czero_vec.data(), n, T(0));
      auto* cone_ptr = cone_vec.data();
      PRAGMA_OFFLOAD("omp target update to(cone_ptr[:n])")
      auto* czero_ptr = czero_vec.data();
      PRAGMA_OFFLOAD("omp target update to(czero_ptr[:n])")
    }
  }

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
    PRAGMA_OFFLOAD("omp target update to(Ainv_ptr[:Ainv.size()])")
  }

  template<typename VVT, typename RATIOT, typename OMPALLOC>
  inline void updateRow(Matrix<T, OMPALLOC>& Ainv, int rowchanged, const VVT& phiV, RATIOT c_ratio_in)
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
    auto* phiV_ptr = phiV.data();
    auto* Ainv_ptr  = Ainv.data();
    auto* temp_ptr  = temp.data();
    auto* rcopy_ptr = rcopy.data();
    PRAGMA_OFFLOAD("omp target data map(always, to: phiV_ptr[:norb]) map(always, from: Ainv_ptr[:Ainv.rows()*Ainv.cols()]) use_device_ptr(phiV_ptr, Ainv_ptr, temp_ptr, rcopy_ptr)")
    {
      success = ompBLAS::gemv(dummy_handle, 'T', norb, norb, cone, Ainv_ptr, norb, phiV_ptr, 1, czero, temp_ptr, 1);
      PRAGMA_OFFLOAD("omp target is_device_ptr(Ainv_ptr, temp_ptr, rcopy_ptr)")
      {
        temp_ptr[rowchanged] -= cone;
        PRAGMA_OFFLOAD("omp parallel for simd")
        for(int i = 0; i < norb; i++)
          rcopy_ptr[i] = Ainv_ptr[rowchanged * norb + i];
      }
      success = ompBLAS::ger(dummy_handle, norb, norb, static_cast<T>(RATIOT(-1)/c_ratio_in), rcopy_ptr, 1, temp_ptr, 1, Ainv_ptr, norb);
    }
  }

  template<typename DEVPTRV, typename VGLV, typename VGV>
  inline void mw_updateRow(const DEVPTRV& Ainv_list, int norb, int rowchanged, const VGLV& phi_vgl_v, const VGV& ratio_grads_v)
  {
    size_t nw = Ainv_list.size();
    temp.resize(norb * nw);
    rcopy.resize(norb * nw);
    c_ratio_inv.resize(nw);
    auto* temp_dev_ptr = getOffloadDevicePtr(temp.data());
    auto* rcopy_dev_ptr = getOffloadDevicePtr(rcopy.data());
    auto* phi_vgl_v_dev_ptr = getOffloadDevicePtr(phi_vgl_v.data());

    // to handle T** of Ainv, psi_v, temp, rcopy
    ptr_buffer.resize(nw * 4);
    for(int iw = 0; iw < nw; iw++)
    {
      ptr_buffer[iw]          = Ainv_list[iw];
      ptr_buffer[iw + nw]     = const_cast<T*>(phi_vgl_v_dev_ptr + norb * iw);
      ptr_buffer[iw + nw * 2] = temp_dev_ptr + norb * iw;
      ptr_buffer[iw + nw * 3] = rcopy_dev_ptr + norb * iw;
      c_ratio_inv[iw]         = T(-1) / ratio_grads_v.data(0)[iw];
    }

    // update the inverse matrix
    constexpr T cone(1);
    constexpr T czero(0);
    int dummy_handle = 0;
    int success =0;
    auto* ptr_buffer_ptr = ptr_buffer.data();
    auto* c_ratio_inv_ptr = c_ratio_inv.data();
    resize_fill_constant_arrays(nw);
    auto* cone_ptr = cone_vec.data();
    auto* czero_ptr = czero_vec.data();
    PRAGMA_OFFLOAD("omp target data map(always, to: ptr_buffer_ptr[:ptr_buffer.size()], c_ratio_inv_ptr[:c_ratio_inv.size()]) \
                    use_device_ptr(ptr_buffer_ptr, c_ratio_inv_ptr, cone_ptr, czero_ptr)")
    {
      // invoke the Fahy's variant of Sherman-Morrison update.
      success = ompBLAS::gemv_batched(dummy_handle, 'T', norb, norb, cone_ptr, ptr_buffer_ptr, norb, ptr_buffer_ptr + nw, 1, czero_ptr, ptr_buffer_ptr + nw * 2, 1, nw);
      PRAGMA_OFFLOAD("omp target teams distribute num_teams(nw) is_device_ptr(ptr_buffer_ptr)")
      for(int iw = 0; iw < nw; iw++)
      {
        auto* Ainv_ptr  = ptr_buffer_ptr[iw];
        auto* temp_ptr  = ptr_buffer_ptr[nw * 2 + iw];
        auto* rcopy_ptr = ptr_buffer_ptr[nw * 3 + iw];

        temp_ptr[rowchanged] -= cone;
        PRAGMA_OFFLOAD("omp parallel for simd")
        for(int i = 0; i < norb; i++)
          rcopy_ptr[i] = Ainv_ptr[rowchanged * norb + i];
      }

      success = ompBLAS::ger_batched(dummy_handle, norb, norb, c_ratio_inv_ptr, ptr_buffer_ptr + nw * 3, 1, ptr_buffer_ptr + nw * 2, 1, ptr_buffer_ptr, norb, nw);
    }
  }

};
} // namespace qmcplusplus

#endif // QMCPLUSPLUS_MATRIX_UPDATE_H
