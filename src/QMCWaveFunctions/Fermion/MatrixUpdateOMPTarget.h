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

#ifndef QMCPLUSPLUS_MATRIX_UPDATE_OMPTARGET_H
#define QMCPLUSPLUS_MATRIX_UPDATE_OMPTARGET_H

#include "CPU/SIMD/aligned_allocator.hpp"
#include "Platforms/PinnedAllocator.h"
#include "OMPTarget/OMPallocator.hpp"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Fermion/DiracMatrixComputeOMPTarget.hpp"
#include "OMPTarget/ompBLAS.hpp"
#include "OMPTarget/ompReduction.hpp"
#include "ResourceCollection.h"
#include "Configuration.h"

namespace qmcplusplus
{
/** Implements dirac matrix update using OpenMP offload.
 * It is used as DET_ENGINE_TYPE in DiracDeterminantBatched.
 * @tparam T base precision for most computation
 * @tparam T_FP high precision for matrix inversion, T_FP >= T
 */
template<typename T, typename T_FP>
class MatrixUpdateOMPTarget
{
  using This_t = MatrixUpdateOMPTarget<T, T_FP>;
  using Real = QMCTraits::RealType;
  template<typename DT>
  using OffloadAllocator = OMPallocator<DT, aligned_allocator<DT>>;
  template<typename DT>
  using OffloadPinnedAllocator     = OMPallocator<DT, PinnedAlignedAllocator<DT>>;
  using OffloadValueVector_t       = Vector<T, OffloadAllocator<T>>;
  using OffloadPinnedValueVector_t = Vector<T, OffloadPinnedAllocator<T>>;
  using OffloadPinnedValueMatrix_t = Matrix<T, OffloadPinnedAllocator<T>>;

  /// matrix inversion engine this crowd scope resouce and only the leader engine gets it
  UPtr<DiracMatrixComputeOMPTarget<T_FP>> det_inverter_;
  /// inverse transpose of psiM(j,i) \f$= \psi_j({\bf r}_i)\f$
  OffloadPinnedValueMatrix_t psiMinv;
  /// scratch space for rank-1 update
  OffloadValueVector_t temp;
  // scratch space for keeping one row of Ainv
  OffloadValueVector_t rcopy;
  // constant array value T(1)
  OffloadValueVector_t cone_vec;
  // constant array value T(0)
  OffloadValueVector_t czero_vec;
  // multi walker of grads for transfer needs.
  OffloadPinnedValueMatrix_t grads_value_v;
  // pointer buffer
  Vector<char, OffloadPinnedAllocator<char>> buffer_H2D;

  void resize_fill_constant_arrays(size_t nw)
  {
    if (cone_vec.size() < nw)
    {
      cone_vec.resize(nw);
      czero_vec.resize(nw);
      std::fill_n(cone_vec.data(), nw, T(1));
      std::fill_n(czero_vec.data(), nw, T(0));
      T* cone_ptr = cone_vec.data();
      PRAGMA_OFFLOAD("omp target update to(cone_ptr[:nw])")
      T* czero_ptr = czero_vec.data();
      PRAGMA_OFFLOAD("omp target update to(czero_ptr[:nw])")
    }
  }

  void resize_scratch_arrays(int norb, size_t nw)
  {
    size_t total_size = norb * nw;
    if (temp.size() < total_size)
    {
      temp.resize(total_size);
      rcopy.resize(total_size);
    }
  }

public:
  /** resize the internal storage
   * @param norb number of electrons/orbitals
   * @param delay, maximum delay 0<delay<=norb
   */
  inline void resize(int norb, int delay) { psiMinv.resize(norb, getAlignedSize<T>(norb)); }

  void createResource(ResourceCollection& collection) {}
  void acquireResource(ResourceCollection& collection) {}
  void releaseResource(ResourceCollection& collection) {}

  OffloadPinnedValueMatrix_t& get_psiMinv() { return psiMinv; }

  inline T* getRow_psiMinv_offload(int row_id) { return psiMinv.device_data() + row_id * psiMinv.cols(); }

  /** compute the inverse of the transpose of matrix logdetT, result is in psiMinv
   * @param logdetT orbital value matrix
   * @param LogValue log(det(logdetT))
   */
  template<typename TREAL>
  inline void invert_transpose(const Matrix<T>& logdetT, std::complex<TREAL>& LogValue)
  {
    auto& Ainv = psiMinv;
    Matrix<T> Ainv_host_view(Ainv.data(), Ainv.rows(), Ainv.cols());
    det_inverter_.invert_transpose(logdetT, Ainv_host_view, LogValue);
    T* Ainv_ptr = Ainv.data();
    PRAGMA_OFFLOAD("omp target update to(Ainv_ptr[:Ainv.size()])")
  }

  inline void mw_invertTranspose(const RefVectorWithLeader<MatrixUpdateOMPTarget<T, T_FP>>& engines,
                                 const RefVector<OffloadPinnedValueMatrix_t>& logdetT_list,
                                 const RefVector<OffloadPinnedValueVector_t>& LogValues)
  {
    if (!det_inverter_)
      det_inverter_ = std::make_unique<DiracMatrixComputeOMPTarget<T_FP>>();

    std::vector<Matrix<T>> Ainv_host_views;
    Ainv_host_views.reserve(engines.size());
    for (int iw = 0; iw < engines.size(); iw++)
    {
      auto& Ainv = engines[iw].get().psiMinv;
      Ainv_host_views.emplace_back(Ainv.data(), Ainv.rows(), Ainv.cols());
      T* Ainv_ptr = Ainv.data();
      PRAGMA_OFFLOAD("omp target update to(Ainv_ptr[:Ainv.size()])")
    }
    PRAGMA_OFFLOAD("omp taskwait")
    det_inverter_->mw_invert_transpose(logdetT_list, Ainv_host_views, LogValues);
  }

  template<typename GT>
  inline void mw_evalGrad(const RefVector<This_t>& engines,
                          const std::vector<const T*>& dpsiM_row_list,
                          int rowchanged,
                          std::vector<GT>& grad_now)
  {
    const int norb = psiMinv.rows();
    const int nw   = engines.size();
    buffer_H2D.resize(sizeof(T*) * 2 * nw);
    Matrix<const T*> ptr_buffer(reinterpret_cast<const T**>(buffer_H2D.data()), 2, nw);
    for (int iw = 0; iw < nw; iw++)
    {
      ptr_buffer[0][iw] = engines[iw].get().psiMinv.device_data() + rowchanged * psiMinv.cols();
      ptr_buffer[1][iw] = dpsiM_row_list[iw];
    }

    constexpr unsigned DIM = GT::Size;
    grads_value_v.resize(nw, DIM);
    auto* __restrict__ grads_value_v_ptr = grads_value_v.data();
    auto* buffer_H2D_ptr                 = buffer_H2D.data();

    PRAGMA_OFFLOAD("omp target teams distribute num_teams(nw) \
                    map(always, to: buffer_H2D_ptr[:buffer_H2D.size()]) \
                    map(always, from: grads_value_v_ptr[:grads_value_v.size()])")
    for (int iw = 0; iw < nw; iw++)
    {
      const T* __restrict__ invRow_ptr    = reinterpret_cast<const T**>(buffer_H2D_ptr)[iw];
      const T* __restrict__ dpsiM_row_ptr = reinterpret_cast<const T**>(buffer_H2D_ptr)[nw + iw];
      T grad_x(0), grad_y(0), grad_z(0);
      PRAGMA_OFFLOAD("omp parallel for reduction(+: grad_x, grad_y, grad_z)")
      for (int iorb = 0; iorb < norb; iorb++)
      {
        grad_x += invRow_ptr[iorb] * dpsiM_row_ptr[iorb * DIM];
        grad_y += invRow_ptr[iorb] * dpsiM_row_ptr[iorb * DIM + 1];
        grad_z += invRow_ptr[iorb] * dpsiM_row_ptr[iorb * DIM + 2];
      }
      grads_value_v_ptr[iw * DIM]     = grad_x;
      grads_value_v_ptr[iw * DIM + 1] = grad_y;
      grads_value_v_ptr[iw * DIM + 2] = grad_z;
    }

    for (int iw = 0; iw < nw; iw++)
      grad_now[iw] = {grads_value_v[iw][0], grads_value_v[iw][1], grads_value_v[iw][2]};
  }

  template<typename VVT, typename RATIOT>
  inline void updateRow(int rowchanged, const VVT& phiV, RATIOT c_ratio_in)
  {
    auto& Ainv = psiMinv;
    // update the inverse matrix
    constexpr T cone(1);
    constexpr T czero(0);
    const int norb = Ainv.rows();
    const int lda  = Ainv.cols();
    resize_scratch_arrays(norb, 1);
    // invoke the Fahy's variant of Sherman-Morrison update.
    int dummy_handle  = 0;
    int success       = 0;
    const T* phiV_ptr = phiV.data();
    T* Ainv_ptr       = Ainv.data();
    T* temp_ptr       = temp.data();
    T* rcopy_ptr      = rcopy.data();
    PRAGMA_OFFLOAD("omp target data map(always, to: phiV_ptr[:norb]) \
                    map(always, from: Ainv_ptr[:Ainv.size()]) \
                    use_device_ptr(phiV_ptr, Ainv_ptr, temp_ptr, rcopy_ptr)")
    {
      success = ompBLAS::gemv(dummy_handle, 'T', norb, norb, cone, Ainv_ptr, lda, phiV_ptr, 1, czero, temp_ptr, 1);
      PRAGMA_OFFLOAD("omp target is_device_ptr(Ainv_ptr, temp_ptr, rcopy_ptr)")
      {
        temp_ptr[rowchanged] -= cone;
        PRAGMA_OFFLOAD("omp parallel for simd")
        for (int i = 0; i < norb; i++)
          rcopy_ptr[i] = Ainv_ptr[rowchanged * lda + i];
      }
      success = ompBLAS::ger(dummy_handle, norb, norb, static_cast<T>(RATIOT(-1) / c_ratio_in), rcopy_ptr, 1, temp_ptr,
                             1, Ainv_ptr, lda);
    }
  }

  inline void mw_updateRow(const RefVector<This_t>& engines,
                           int rowchanged,
                           const std::vector<T*>& psiM_g_list,
                           const std::vector<T*>& psiM_l_list,
                           const std::vector<bool>& isAccepted,
                           const T* phi_vgl_v_dev_ptr,
                           const size_t phi_vgl_stride,
                           const std::vector<T>& ratios)
  {
    const int norb          = psiMinv.rows();
    const int lda           = psiMinv.cols();
    const size_t n_accepted = psiM_g_list.size();
    if (n_accepted == 0)
      return;

    resize_scratch_arrays(norb, n_accepted);

    // to handle T** of Ainv, psi_v, temp, rcopy
    buffer_H2D.resize((sizeof(T*) * 8 + sizeof(T)) * n_accepted);
    Matrix<T*> ptr_buffer(reinterpret_cast<T**>(buffer_H2D.data()), 8, n_accepted);
    T* c_ratio_inv = reinterpret_cast<T*>(buffer_H2D.data() + sizeof(T*) * 8 * n_accepted);
    for (int iw = 0, count = 0; iw < isAccepted.size(); iw++)
      if (isAccepted[iw])
      {
        ptr_buffer[0][count] = engines[iw].get().psiMinv.device_data();
        ptr_buffer[1][count] = const_cast<T*>(phi_vgl_v_dev_ptr + norb * iw);
        ptr_buffer[2][count] = temp.device_data() + norb * count;
        ptr_buffer[3][count] = rcopy.device_data() + norb * count;
        ptr_buffer[4][count] = psiM_g_list[count];
        ptr_buffer[5][count] = psiM_l_list[count];
        ptr_buffer[6][count] = const_cast<T*>(phi_vgl_v_dev_ptr + phi_vgl_stride + norb * 3 * iw);
        ptr_buffer[7][count] = const_cast<T*>(phi_vgl_v_dev_ptr + phi_vgl_stride * 4 + norb * iw);

        c_ratio_inv[count] = T(-1) / ratios[iw];
        count++;
      }

    // update the inverse matrix
    constexpr T cone(1);
    constexpr T czero(0);
    int dummy_handle     = 0;
    int success          = 0;
    auto* buffer_H2D_ptr = buffer_H2D.data();
    resize_fill_constant_arrays(n_accepted);
    T* cone_ptr  = cone_vec.data();
    T* czero_ptr = czero_vec.data();
    PRAGMA_OFFLOAD("omp target data \
                    map(always, to: buffer_H2D_ptr[:buffer_H2D.size()]) \
                    use_device_ptr(buffer_H2D_ptr, cone_ptr, czero_ptr)")
    {
      T** Ainv_mw_ptr   = reinterpret_cast<T**>(buffer_H2D_ptr);
      T** phiV_mw_ptr   = reinterpret_cast<T**>(buffer_H2D_ptr + sizeof(T*) * n_accepted);
      T** temp_mw_ptr   = reinterpret_cast<T**>(buffer_H2D_ptr + sizeof(T*) * n_accepted * 2);
      T** rcopy_mw_ptr  = reinterpret_cast<T**>(buffer_H2D_ptr + sizeof(T*) * n_accepted * 3);
      T** dpsiM_mw_out  = reinterpret_cast<T**>(buffer_H2D_ptr + sizeof(T*) * n_accepted * 4);
      T** d2psiM_mw_out = reinterpret_cast<T**>(buffer_H2D_ptr + sizeof(T*) * n_accepted * 5);
      T** dpsiM_mw_in   = reinterpret_cast<T**>(buffer_H2D_ptr + sizeof(T*) * n_accepted * 6);
      T** d2psiM_mw_in  = reinterpret_cast<T**>(buffer_H2D_ptr + sizeof(T*) * n_accepted * 7);
      T* ratio_inv_mw   = reinterpret_cast<T*>(buffer_H2D_ptr + sizeof(T*) * n_accepted * 8);

      // invoke the Fahy's variant of Sherman-Morrison update.
      success = ompBLAS::gemv_batched(dummy_handle, 'T', norb, norb, cone_ptr, Ainv_mw_ptr, lda, phiV_mw_ptr, 1,
                                      czero_ptr, temp_mw_ptr, 1, n_accepted);
      PRAGMA_OFFLOAD("omp target teams distribute num_teams(n_accepted) is_device_ptr(Ainv_mw_ptr, temp_mw_ptr, \
                     rcopy_mw_ptr, dpsiM_mw_out, d2psiM_mw_out, dpsiM_mw_in, d2psiM_mw_in)")
      for (int iw = 0; iw < n_accepted; iw++)
      {
        T* __restrict__ Ainv_ptr   = Ainv_mw_ptr[iw];
        T* __restrict__ temp_ptr   = temp_mw_ptr[iw];
        T* __restrict__ rcopy_ptr  = rcopy_mw_ptr[iw];
        T* __restrict__ dpsiM_out  = dpsiM_mw_out[iw];
        T* __restrict__ d2psiM_out = d2psiM_mw_out[iw];
        T* __restrict__ dpsiM_in   = dpsiM_mw_in[iw];
        T* __restrict__ d2psiM_in  = d2psiM_mw_in[iw];

        temp_ptr[rowchanged] -= cone;
        PRAGMA_OFFLOAD("omp parallel for simd")
        for (int i = 0; i < norb; i++)
        {
          rcopy_ptr[i] = Ainv_ptr[rowchanged * lda + i];
          // the following copying data on the device is not part of SM-1
          // it is intended to copy dpsiM and d2psiM from temporary to final without a separate kernel.
          dpsiM_out[i * 3]     = dpsiM_in[i * 3];
          dpsiM_out[i * 3 + 1] = dpsiM_in[i * 3 + 1];
          dpsiM_out[i * 3 + 2] = dpsiM_in[i * 3 + 2];
          d2psiM_out[i]        = d2psiM_in[i];
        }
      }

      success = ompBLAS::ger_batched(dummy_handle, norb, norb, ratio_inv_mw, rcopy_mw_ptr, 1, temp_mw_ptr, 1,
                                     Ainv_mw_ptr, lda, n_accepted);
    }
  }

  inline void mw_accept_rejectRow(const RefVector<This_t>& engines,
                                  const int rowchanged,
                                  const std::vector<T*>& psiM_g_list,
                                  const std::vector<T*>& psiM_l_list,
                                  const std::vector<bool>& isAccepted,
                                  const T* phi_vgl_v_dev_ptr,
                                  const size_t phi_vgl_stride,
                                  const std::vector<T>& ratios)
  {
    mw_updateRow(engines, rowchanged, psiM_g_list, psiM_l_list, isAccepted, phi_vgl_v_dev_ptr, phi_vgl_stride, ratios);
  }

  /** update the full Ainv and reset delay_count
   * @param Ainv inverse matrix
   */
  inline void mw_updateInvMat(const RefVector<This_t>& engines) {}

  std::vector<const T*> mw_getInvRow(const RefVector<This_t>& engines, const int row_id, bool on_host) const
  {
    const size_t nw = engines.size();
    std::vector<const T*> row_ptr_list;
    row_ptr_list.reserve(nw);
    if (on_host)
    {
      for (This_t& engine : engines)
      {
        auto* ptr = engine.psiMinv.data();
        PRAGMA_OFFLOAD("omp target update from(ptr[row_id * psiMinv.cols():psiMinv.cols()])")
        row_ptr_list.push_back(ptr + row_id * psiMinv.cols());
      }
    }
    else
    {
      for (This_t& engine : engines)
        row_ptr_list.push_back(engine.psiMinv.device_data() + row_id * psiMinv.cols());
    }
    return row_ptr_list;
  }

  void mw_transferAinv_D2H(const RefVector<This_t>& engines)
  {
    for (This_t& engine : engines)
    {
      auto* ptr = engine.psiMinv.data();
      PRAGMA_OFFLOAD("omp target update from(ptr[:psiMinv.size()])")
    }
  }

  DiracMatrixComputeOMPTarget<T_FP>& get_det_inverter()
  {
    if (det_inverter_)
      return *det_inverter_;
    throw std::logic_error("attempted to get null det_inverter_, this is developer logic error");
  }
};
} // namespace qmcplusplus

#endif // QMCPLUSPLUS_MATRIX_UPDATE_H
