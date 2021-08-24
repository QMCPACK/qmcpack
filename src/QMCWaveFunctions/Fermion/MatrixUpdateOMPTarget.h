//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_MATRIX_UPDATE_OMPTARGET_H
#define QMCPLUSPLUS_MATRIX_UPDATE_OMPTARGET_H

#include "OMPTarget/OffloadAlignedAllocators.hpp"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "QMCWaveFunctions/Fermion/DiracMatrix.h"
#include "OMPTarget/ompBLAS.hpp"
#include "OMPTarget/ompReduction.hpp"
#include "ResourceCollection.h"


namespace qmcplusplus
{
template<typename T>
struct MatrixUpdateOMPTargetMultiWalkerMem : public Resource
{
  using OffloadValueVector_t       = Vector<T, OffloadAllocator<T>>;
  using OffloadPinnedValueVector_t = Vector<T, OffloadPinnedAllocator<T>>;
  using OffloadPinnedValueMatrix_t = Matrix<T, OffloadPinnedAllocator<T>>;

  // constant array value T(1)
  OffloadValueVector_t cone_vec;
  // constant array value T(0)
  OffloadValueVector_t czero_vec;
  // multi walker of grads for transfer needs.
  OffloadPinnedValueMatrix_t grads_value_v;
  // pointer buffer
  Vector<char, OffloadPinnedAllocator<char>> buffer_H2D;
  /// scratch space for rank-1 update
  OffloadValueVector_t mw_temp;
  // scratch space for keeping one row of Ainv
  OffloadValueVector_t mw_rcopy;

  MatrixUpdateOMPTargetMultiWalkerMem() : Resource("MatrixUpdateOMPTargetMultiWalkerMem") {}

  MatrixUpdateOMPTargetMultiWalkerMem(const MatrixUpdateOMPTargetMultiWalkerMem&)
      : MatrixUpdateOMPTargetMultiWalkerMem()
  {}

  Resource* makeClone() const override { return new MatrixUpdateOMPTargetMultiWalkerMem(*this); }
};

/** Implements dirac matrix update using OpenMP offload.
 * It is used as DET_ENGINE in DiracDeterminantBatched.
 * @tparam T base precision for most computation
 * @tparam T_FP high precision for matrix inversion, T_FP >= T
 */
template<typename T, typename T_FP>
class MatrixUpdateOMPTarget
{
  using This_t = MatrixUpdateOMPTarget<T, T_FP>;

  using OffloadValueVector_t       = Vector<T, OffloadAllocator<T>>;
  using OffloadPinnedValueVector_t = Vector<T, OffloadPinnedAllocator<T>>;
  using OffloadPinnedValueMatrix_t = Matrix<T, OffloadPinnedAllocator<T>>;

  /// matrix inversion engine
  DiracMatrix<T_FP> detEng;
  /// inverse transpose of psiM(j,i) \f$= \psi_j({\bf r}_i)\f$
  OffloadPinnedValueMatrix_t psiMinv;
  /// scratch space for rank-1 update
  OffloadValueVector_t temp;
  // scratch space for keeping one row of Ainv
  OffloadValueVector_t rcopy;

  // multi walker memory buffers
  std::unique_ptr<MatrixUpdateOMPTargetMultiWalkerMem<T>> mw_mem_;

  void resize_fill_constant_arrays(size_t nw)
  {
    if (mw_mem_->cone_vec.size() < nw)
    {
      mw_mem_->cone_vec.resize(nw);
      mw_mem_->czero_vec.resize(nw);
      std::fill_n(mw_mem_->cone_vec.data(), nw, T(1));
      std::fill_n(mw_mem_->czero_vec.data(), nw, T(0));
      T* cone_ptr = mw_mem_->cone_vec.data();
      PRAGMA_OFFLOAD("omp target update to(cone_ptr[:nw])")
      T* czero_ptr = mw_mem_->czero_vec.data();
      PRAGMA_OFFLOAD("omp target update to(czero_ptr[:nw])")
    }
  }

  void resize_scratch_arrays(int norb, size_t nw)
  {
    size_t total_size = norb * nw;
    if (mw_mem_->mw_temp.size() < total_size)
    {
      mw_mem_->mw_temp.resize(total_size);
      mw_mem_->mw_rcopy.resize(total_size);
    }
  }

public:
  /** resize the internal storage
   * @param norb number of electrons/orbitals
   * @param delay, maximum delay 0<delay<=norb
   */
  inline void resize(int norb, int delay) { psiMinv.resize(norb, getAlignedSize<T>(norb)); }

  void createResource(ResourceCollection& collection) const
  {
    collection.addResource(std::make_unique<MatrixUpdateOMPTargetMultiWalkerMem<T>>());
  }

  void acquireResource(ResourceCollection& collection)
  {
    auto res_ptr = dynamic_cast<MatrixUpdateOMPTargetMultiWalkerMem<T>*>(collection.lendResource().release());
    if (!res_ptr)
      throw std::runtime_error(
          "MatrixUpdateOMPTarget::acquireResource dynamic_cast MatrixUpdateOMPTargetMultiWalkerMem failed");
    mw_mem_.reset(res_ptr);
  }

  void releaseResource(ResourceCollection& collection) { collection.takebackResource(std::move(mw_mem_)); }

  OffloadPinnedValueMatrix_t& get_psiMinv() { return psiMinv; }

  inline T* getRow_psiMinv_offload(int row_id) { return psiMinv.device_data() + row_id * psiMinv.cols(); }

  /** compute the inverse of the transpose of matrix logdetT, result is in psiMinv
   * @param logdetT orbital value matrix
   * @param log_value log(det(logdetT))
   */
  template<typename TREAL>
  inline void invert_transpose(const Matrix<T>& logdetT, std::complex<TREAL>& log_value)
  {
    auto& Ainv = psiMinv;
    Matrix<T> Ainv_host_view(Ainv.data(), Ainv.rows(), Ainv.cols());
    detEng.invert_transpose(logdetT, Ainv_host_view, log_value);
    T* Ainv_ptr = Ainv.data();
    PRAGMA_OFFLOAD("omp target update to(Ainv_ptr[:Ainv.size()])")
  }

  template<typename TREAL>
  static void mw_invert_transpose(const RefVectorWithLeader<This_t>& engines,
                                  const RefVector<const Matrix<T>>& logdetT_list,
                                  const RefVector<std::complex<TREAL>>& log_values)
  {
    auto& engine_leader = engines.getLeader();
    // make this class unit tests friendly without the need of setup resources.
    if (!engine_leader.mw_mem_)
    {
      app_warning() << "MatrixUpdateOMPTarget : This message should not be seen in production (performance bug) runs "
                       "but only unit tests (expected)."
                    << std::endl;
      engine_leader.mw_mem_ = std::make_unique<MatrixUpdateOMPTargetMultiWalkerMem<T>>();
    }

    for (int iw = 0; iw < engines.size(); iw++)
    {
      auto& Ainv = engines[iw].psiMinv;
      Matrix<T> Ainv_host_view(Ainv.data(), Ainv.rows(), Ainv.cols());
      engine_leader.detEng.invert_transpose(logdetT_list[iw].get(), Ainv_host_view, log_values[iw].get());
      T* Ainv_ptr = Ainv.data();
      PRAGMA_OFFLOAD("omp target update to(Ainv_ptr[:Ainv.size()])")
    }
    PRAGMA_OFFLOAD("omp taskwait")
  }

  template<typename GT>
  static void mw_evalGrad(const RefVectorWithLeader<This_t>& engines,
                          const std::vector<const T*>& dpsiM_row_list,
                          int rowchanged,
                          std::vector<GT>& grad_now)
  {
    auto& engine_leader = engines.getLeader();
    auto& buffer_H2D    = engine_leader.mw_mem_->buffer_H2D;
    auto& grads_value_v = engine_leader.mw_mem_->grads_value_v;

    const int norb = engine_leader.psiMinv.rows();
    const int nw   = engines.size();
    buffer_H2D.resize(sizeof(T*) * 2 * nw);
    Matrix<const T*> ptr_buffer(reinterpret_cast<const T**>(buffer_H2D.data()), 2, nw);
    for (int iw = 0; iw < nw; iw++)
    {
      ptr_buffer[0][iw] = engines[iw].psiMinv.device_data() + rowchanged * engine_leader.psiMinv.cols();
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
    temp.resize(norb);
    rcopy.resize(norb);
    // invoke the Fahy's variant of Sherman-Morrison update.
    int dummy_handle  = 0;
    const T* phiV_ptr = phiV.data();
    T* Ainv_ptr       = Ainv.data();
    T* temp_ptr       = temp.data();
    T* rcopy_ptr      = rcopy.data();
    PRAGMA_OFFLOAD("omp target data map(always, to: phiV_ptr[:norb]) \
                    map(always, from: Ainv_ptr[:Ainv.size()]) \
                    use_device_ptr(phiV_ptr, Ainv_ptr, temp_ptr, rcopy_ptr)")
    {
      ompBLAS::gemv(dummy_handle, 'T', norb, norb, cone, Ainv_ptr, lda, phiV_ptr, 1, czero, temp_ptr, 1);
      PRAGMA_OFFLOAD("omp target is_device_ptr(Ainv_ptr, temp_ptr, rcopy_ptr)")
      {
        temp_ptr[rowchanged] -= cone;
        PRAGMA_OFFLOAD("omp parallel for simd")
        for (int i = 0; i < norb; i++)
          rcopy_ptr[i] = Ainv_ptr[rowchanged * lda + i];
      }
      ompBLAS::ger(dummy_handle, norb, norb, static_cast<T>(RATIOT(-1) / c_ratio_in), rcopy_ptr, 1, temp_ptr, 1,
                   Ainv_ptr, lda);
    }
  }

  static void mw_updateRow(const RefVectorWithLeader<This_t>& engines,
                           int rowchanged,
                           const std::vector<T*>& psiM_g_list,
                           const std::vector<T*>& psiM_l_list,
                           const std::vector<bool>& isAccepted,
                           const T* phi_vgl_v_dev_ptr,
                           const size_t phi_vgl_stride,
                           const std::vector<T>& ratios)
  {
    const size_t n_accepted = psiM_g_list.size();
    if (n_accepted == 0)
      return;

    auto& engine_leader = engines.getLeader();
    auto& buffer_H2D    = engine_leader.mw_mem_->buffer_H2D;
    auto& grads_value_v = engine_leader.mw_mem_->grads_value_v;
    auto& cone_vec      = engine_leader.mw_mem_->cone_vec;
    auto& czero_vec     = engine_leader.mw_mem_->czero_vec;
    auto& mw_temp       = engine_leader.mw_mem_->mw_temp;
    auto& mw_rcopy      = engine_leader.mw_mem_->mw_rcopy;
    const int norb      = engine_leader.psiMinv.rows();
    const int lda       = engine_leader.psiMinv.cols();

    engine_leader.resize_scratch_arrays(norb, n_accepted);

    // to handle T** of Ainv, psi_v, temp, rcopy
    buffer_H2D.resize((sizeof(T*) * 8 + sizeof(T)) * n_accepted);
    Matrix<T*> ptr_buffer(reinterpret_cast<T**>(buffer_H2D.data()), 8, n_accepted);
    T* c_ratio_inv = reinterpret_cast<T*>(buffer_H2D.data() + sizeof(T*) * 8 * n_accepted);
    for (int iw = 0, count = 0; iw < isAccepted.size(); iw++)
      if (isAccepted[iw])
      {
        ptr_buffer[0][count] = engines[iw].psiMinv.device_data();
        ptr_buffer[1][count] = const_cast<T*>(phi_vgl_v_dev_ptr + norb * iw);
        ptr_buffer[2][count] = mw_temp.device_data() + norb * count;
        ptr_buffer[3][count] = mw_rcopy.device_data() + norb * count;
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
    auto* buffer_H2D_ptr = buffer_H2D.data();
    engine_leader.resize_fill_constant_arrays(n_accepted);
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
      ompBLAS::gemv_batched(dummy_handle, 'T', norb, norb, cone_ptr, Ainv_mw_ptr, lda, phiV_mw_ptr, 1, czero_ptr,
                            temp_mw_ptr, 1, n_accepted);
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

      ompBLAS::ger_batched(dummy_handle, norb, norb, ratio_inv_mw, rcopy_mw_ptr, 1, temp_mw_ptr, 1, Ainv_mw_ptr, lda,
                           n_accepted);
    }
  }

  static void mw_accept_rejectRow(const RefVectorWithLeader<This_t>& engines,
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
  inline static void mw_updateInvMat(const RefVectorWithLeader<This_t>& engines) {}

  std::vector<const T*> static mw_getInvRow(const RefVectorWithLeader<This_t>& engines, const int row_id, bool on_host)
  {
    const size_t nw    = engines.size();
    const size_t ncols = engines.getLeader().psiMinv.cols();
    std::vector<const T*> row_ptr_list;
    row_ptr_list.reserve(nw);
    if (on_host)
      for (This_t& engine : engines)
      {
        auto* ptr = engine.psiMinv.data();
        PRAGMA_OFFLOAD("omp target update from(ptr[row_id * ncols : ncols])")
        row_ptr_list.push_back(ptr + row_id * ncols);
      }
    else
      for (This_t& engine : engines)
        row_ptr_list.push_back(engine.psiMinv.device_data() + row_id * ncols);
    return row_ptr_list;
  }

  static void mw_transferAinv_D2H(const RefVectorWithLeader<This_t>& engines)
  {
    for (This_t& engine : engines)
    {
      auto* ptr = engine.psiMinv.data();
      PRAGMA_OFFLOAD("omp target update from(ptr[:engine.psiMinv.size()])")
    }
  }
};
} // namespace qmcplusplus

#endif // QMCPLUSPLUS_MATRIX_UPDATE_H
