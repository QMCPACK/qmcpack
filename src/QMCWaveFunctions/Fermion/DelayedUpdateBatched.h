//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_DELAYED_UPDATE_BATCHED_H
#define QMCPLUSPLUS_DELAYED_UPDATE_BATCHED_H

#include <cassert>
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "QMCWaveFunctions/Fermion/DiracMatrix.h"
#include "Platforms/OMPTarget/ompBLAS.hpp"
#include "detail/AccelMatrixUpdate.hpp"
#include "WaveFunctionTypes.hpp"
#include "QueueAliases.hpp"
#include "AccelBLAS.hpp"
#include "OMPTarget/OffloadAlignedAllocators.hpp"

namespace qmcplusplus
{

/** implements dirac matrix delayed update using OpenMP offload and CUDA.
 * It is used as DET_ENGINE in DiracDeterminantBatched.
 * This is a 1 per crowd class. It is owned by DiracDeterminantBatchedMultiWalkerResource
 * and holds the data of all the walkers in the crowd in single storages.
 * The walker index is the slowest index of every per walker buffer.
 *
 * @tparam T base precision for most computation
 * @tparam T_FP high precision for matrix inversion, T_FP >= T
 */
template<PlatformKind PL, typename VALUE>
class DelayedUpdateBatched
{
public:
  using This_t  = DelayedUpdateBatched<PL, VALUE>;
  using Value   = VALUE;
  using Real    = RealAlias<Value>;
  using Complex = std::complex<Real>;

  // containers
  template<typename DT>
  using UnpinnedDualVector = Vector<DT, OffloadAllocator<DT>>;
  template<typename DT>
  using UnpinnedDualMatrix = Matrix<DT, OffloadAllocator<DT>>;
  template<typename DT>
  using DualVector = Vector<DT, OffloadPinnedAllocator<DT>>;
  template<typename DT>
  using DualMatrix = Matrix<DT, OffloadPinnedAllocator<DT>>;
  template<typename DT>
  using DualVGLVector = VectorSoaContainer<DT, QMCTraits::DIM + 2, OffloadPinnedAllocator<DT>>;
  template<typename DT>
  using OffloadMWVGLArray = Array<DT, 3, OffloadPinnedAllocator<DT>>; // [VGL, walker, Orbs]
  template<typename DT>
  using OffloadMatrix = Matrix<DT, OffloadPinnedAllocator<DT>>;

  // CUDA stream, cublas handle object
  compute::Queue<PL> queue;
  compute::BLASHandle<PL> blas_handle;

private:
  // constant array value VALUE(1)
  UnpinnedDualVector<Value> cone_vec;
  // constant array value VALUE(-1)
  UnpinnedDualVector<Value> cminusone_vec;
  // constant array value VALUE(0)
  UnpinnedDualVector<Value> czero_vec;
  // multi walker of grads for transfer needs.
  DualMatrix<Value> grads_value_v;
  // multi walker of spingrads for transfer needs.
  DualVector<Complex> spingrads_value_v;
  // mw_updateRow pointer buffer
  Vector<char, OffloadPinnedAllocator<char>> updateRow_buffer_H2D;
  // mw_prepareInvRow pointer buffer
  Vector<char, OffloadPinnedAllocator<char>> prepare_inv_row_buffer_H2D;
  // mw_accept_rejectRow pointer buffer
  Vector<char, OffloadPinnedAllocator<char>> accept_rejectRow_buffer_H2D;
  // mw_updateInv pointer buffer
  Vector<char, OffloadPinnedAllocator<char>> updateInv_buffer_H2D;
  // mw_evalGrad pointer buffer
  Vector<char, OffloadPinnedAllocator<char>> evalGrad_buffer_H2D;
  /// scratch space for rank-1 update, [n_accepted][norb]
  UnpinnedDualVector<Value> mw_temp;
  /// scratch space for keeping one row of Ainv, [n_accepted][norb]
  UnpinnedDualVector<Value> mw_rcopy;

  /// rows of up-to-date Ainv of all the walkers, [nw][norb]
  UnpinnedDualMatrix<Value> invRow;
  /** row id correspond to the up-to-date invRow. [0 norb), invRow is ready; -1, invRow is not valid.
   *  This id is set after calling getInvRow indicating invRow has been prepared for the invRow_id row
   *  ratioGrad checks if invRow_id is consistent. If not, invRow needs to be recomputed.
   *  acceptMove and completeUpdates mark invRow invalid by setting invRow_id to -1
   */
  int invRow_id = -1;

  template<typename DT>
  using DeviceMatrix = Matrix<DT, OffloadDeviceAllocator<DT>>;
  /// orbital values of delayed electrons, [nw][max_delay * norb]
  DeviceMatrix<Value> U_gpu;
  /// rows of Ainv corresponding to delayed electrons, [nw][max_delay * norb]
  DeviceMatrix<Value> V_gpu;
  /// Matrix inverse of B, at maximum KxK, [nw][max_delay * max_delay]
  DeviceMatrix<Value> Binv_gpu;
  /// scratch space, used during inverse update, [nw][norb * max_delay]
  DeviceMatrix<Value> tempMat_gpu;
  /// new column of B, [nw][max_delay]
  DeviceMatrix<Value> p_gpu;
  /// list of delayed electrons, [nw][max_delay]
  DeviceMatrix<int> delay_list_gpu;
  /// current number of delays, increase one for each acceptance, reset to 0 after updating Ainv
  int delay_count = 0;
  /// if true, updates are not delayed.
  bool no_delayed_update_ = true;
  /// maximum delay
  size_t max_delay_ = 0;

  /// the slot-th array of n packed pointers (or values) in the device copy of a packed H2D buffer
  template<typename T>
  static T* devPtr(Vector<char, OffloadPinnedAllocator<char>>& buf, size_t slot, size_t n)
  { return reinterpret_cast<T*>(buf.device_data() + sizeof(Value*) * n * slot); }

  /// device pointer to the iw-th walker's invRow
  Value* invRow_at(size_t iw) { return invRow.device_data() + invRow.cols() * iw; }

  void resize_fill_constant_arrays(size_t nw)
  {
    if (cone_vec.size() < nw)
    {
      // cone
      cone_vec.resize(nw);
      std::fill_n(cone_vec.data(), nw, Value(1));
      cone_vec.updateTo();
      // cminusone
      cminusone_vec.resize(nw);
      std::fill_n(cminusone_vec.data(), nw, Value(-1));
      cminusone_vec.updateTo();
      // czero
      czero_vec.resize(nw);
      std::fill_n(czero_vec.data(), nw, Value(0));
      czero_vec.updateTo();
    }
  }

  /** ensure no previous delay left.
   *  This looks like it should be an assert
   */
  inline void guard_no_delay() const
  {
    if (delay_count != 0)
      throw std::runtime_error("BUG: unexpected call sequence delay_count is not 0");
  }

  /** compute the row of up-to-date Ainv
   * @param psiMinv_refs inverse matrices of all the walkers
   * @param rowchanged the row id corresponding to the proposed electron
   */
  void mw_prepareInvRow(const RefVector<DualMatrix<Value>>& psiMinv_refs, const int rowchanged)
  {
    const int nw   = invRow.rows();
    const int norb = invRow.cols();
    assert(nw == psiMinv_refs.size());

    constexpr size_t num_ptrs_packed = 7; // it must match packing and unpacking
    prepare_inv_row_buffer_H2D.resize(sizeof(Value*) * num_ptrs_packed * nw);
    resize_fill_constant_arrays(nw);

    const int lda_Binv = max_delay_;
    Matrix<Value*> ptr_buffer(reinterpret_cast<Value**>(prepare_inv_row_buffer_H2D.data()), num_ptrs_packed, nw);
    for (int iw = 0; iw < nw; iw++)
    {
      DualMatrix<Value>& psiMinv = psiMinv_refs[iw];
      ptr_buffer[0][iw]          = psiMinv.device_data() + rowchanged * psiMinv.cols();
      ptr_buffer[1][iw]          = invRow_at(iw);
      ptr_buffer[2][iw]          = U_gpu[iw];
      ptr_buffer[3][iw]          = p_gpu[iw];
      ptr_buffer[4][iw]          = Binv_gpu[iw];
      ptr_buffer[5][iw]          = Binv_gpu[iw] + delay_count * lda_Binv;
      ptr_buffer[6][iw]          = V_gpu[iw];
    }

    queue.enqueueH2D(prepare_inv_row_buffer_H2D);

    Value** oldRow_mw_ptr  = devPtr<Value*>(prepare_inv_row_buffer_H2D, 0, nw);
    Value** invRow_mw_ptr  = devPtr<Value*>(prepare_inv_row_buffer_H2D, 1, nw);
    Value** U_mw_ptr       = devPtr<Value*>(prepare_inv_row_buffer_H2D, 2, nw);
    Value** p_mw_ptr       = devPtr<Value*>(prepare_inv_row_buffer_H2D, 3, nw);
    Value** Binv_mw_ptr    = devPtr<Value*>(prepare_inv_row_buffer_H2D, 4, nw);
    Value** BinvRow_mw_ptr = devPtr<Value*>(prepare_inv_row_buffer_H2D, 5, nw);
    Value** V_mw_ptr       = devPtr<Value*>(prepare_inv_row_buffer_H2D, 6, nw);

    // save Ainv[rowchanged] to invRow
    //std::copy_n(Ainv[rowchanged], norb, invRow.data());
    compute::BLAS::copy_batched(blas_handle, norb, oldRow_mw_ptr, 1, invRow_mw_ptr, 1, nw);
    // multiply V (NxK) Binv(KxK) U(KxN) invRow right to the left
    //BLAS::gemv('T', norb, delay_count, cone, U_gpu.data(), norb, invRow.data(), 1, czero, p_gpu.data(), 1);
    //BLAS::gemv('N', delay_count, delay_count, -cone, Binv.data(), lda_Binv, p.data(), 1, czero, Binv[delay_count], 1);
    //BLAS::gemv('N', norb, delay_count, cone, V.data(), norb, Binv[delay_count], 1, cone, invRow.data(), 1);
    compute::BLAS::gemv_batched(blas_handle, 'T', norb, delay_count, cone_vec.device_data(), U_mw_ptr, norb,
                                invRow_mw_ptr, 1, czero_vec.device_data(), p_mw_ptr, 1, nw);
    compute::BLAS::gemv_batched(blas_handle, 'N', delay_count, delay_count, cminusone_vec.device_data(), Binv_mw_ptr,
                                lda_Binv, p_mw_ptr, 1, czero_vec.device_data(), BinvRow_mw_ptr, 1, nw);
    compute::BLAS::gemv_batched(blas_handle, 'N', norb, delay_count, cone_vec.device_data(), V_mw_ptr, norb,
                                BinvRow_mw_ptr, 1, cone_vec.device_data(), invRow_mw_ptr, 1, nw);
    // mark row prepared
    invRow_id = rowchanged;
  }

  /** Do complete row updates
   *  many of these const arguments provide pointers or references
   *  somewhere in here is an update that doesn't get where it belongs resulting in a 0
   *  gradient later.
   *  Sad example of OpenMP target code that is far from clear and a poor substitute for a
   *  clear CPU reference implementation.
   *
   *  \param[in] psiMinv_refs
   *  \param[in] rowchanged
   *  \param[in] psiM_g_list        device ptrs
   *  \param[in] psiM_l_list        device ptrs
   *  \param[in] isAccepted         bool but wait some lists are also filtered
   *  \param[in] phi_vgl_v          multiple walker orbital VGL
   *  \param[inout] ratios
   */
  void mw_updateRow(const RefVector<DualMatrix<Value>>& psiMinv_refs,
                    const int rowchanged,
                    const std::vector<Value*>& psiM_g_list,
                    const std::vector<Value*>& psiM_l_list,
                    const std::vector<bool>& isAccepted,
                    const OffloadMWVGLArray<Value>& phi_vgl_v,
                    const std::vector<Value>& ratios)
  {
    const int nw   = invRow.rows();
    const int norb = invRow.cols();
    assert(nw == psiMinv_refs.size());
    assert(nw == isAccepted.size());
    assert(nw == ratios.size());
    assert(nw == phi_vgl_v.size(1));
    guard_no_delay();

    const size_t n_accepted = psiM_g_list.size();
#ifndef NDEBUG
    size_t n_true = std::count_if(isAccepted.begin(), isAccepted.end(), [](bool accepted) { return accepted; });
    assert(n_accepted == n_true);
#endif
    if (n_accepted == 0)
      return;

    const int lda               = psiMinv_refs[0].get().cols();
    const size_t phi_vgl_stride = nw * norb;
    mw_temp.resize(norb * n_accepted);
    mw_rcopy.resize(norb * n_accepted);

    constexpr size_t num_ptrs_packed = 6; // it must match packing and unpacking
    updateRow_buffer_H2D.resize((sizeof(Value*) * num_ptrs_packed + sizeof(Value)) * n_accepted);

    // to handle T** of Ainv, psi_v, temp, rcopy
    Matrix<Value*> ptr_buffer(reinterpret_cast<Value**>(updateRow_buffer_H2D.data()), num_ptrs_packed, n_accepted);
    Value* c_ratio_inv =
        reinterpret_cast<Value*>(updateRow_buffer_H2D.data() + sizeof(Value*) * num_ptrs_packed * n_accepted);
    for (int iw = 0, count = 0; iw < isAccepted.size(); iw++)
      if (isAccepted[iw])
      {
        ptr_buffer[0][count] = psiMinv_refs[iw].get().device_data();
        ptr_buffer[1][count] = const_cast<Value*>(phi_vgl_v.device_data_at(0, iw, 0));
        ptr_buffer[2][count] = mw_temp.device_data() + norb * count;
        ptr_buffer[3][count] = mw_rcopy.device_data() + norb * count;
        ptr_buffer[4][count] = psiM_g_list[count];
        ptr_buffer[5][count] = psiM_l_list[count];

        c_ratio_inv[count] = Value(-1) / ratios[iw];
        count++;
      }

    // update the inverse matrix
    resize_fill_constant_arrays(n_accepted);

    queue.enqueueH2D(updateRow_buffer_H2D);

    {
      Value** Ainv_mw_ptr   = devPtr<Value*>(updateRow_buffer_H2D, 0, n_accepted);
      Value** phiVGL_mw_ptr = devPtr<Value*>(updateRow_buffer_H2D, 1, n_accepted);
      Value** temp_mw_ptr   = devPtr<Value*>(updateRow_buffer_H2D, 2, n_accepted);
      Value** rcopy_mw_ptr  = devPtr<Value*>(updateRow_buffer_H2D, 3, n_accepted);
      Value** dpsiM_mw_out  = devPtr<Value*>(updateRow_buffer_H2D, 4, n_accepted);
      Value** d2psiM_mw_out = devPtr<Value*>(updateRow_buffer_H2D, 5, n_accepted);
      Value* ratio_inv_mw   = devPtr<Value>(updateRow_buffer_H2D, 6, n_accepted);


      // invoke the Fahy's variant of Sherman-Morrison update.
      compute::BLAS::gemv_batched(blas_handle, 'T', norb, norb, cone_vec.device_data(), Ainv_mw_ptr, lda, phiVGL_mw_ptr,
                                  1, czero_vec.device_data(), temp_mw_ptr, 1, n_accepted);

      compute::copyAinvRow_saveGL_batched(queue, rowchanged, norb, Ainv_mw_ptr, lda, temp_mw_ptr, rcopy_mw_ptr,
                                          phiVGL_mw_ptr, phi_vgl_stride, dpsiM_mw_out, d2psiM_mw_out, n_accepted);


      compute::BLAS::ger_batched(blas_handle, norb, norb, ratio_inv_mw, rcopy_mw_ptr, 1, temp_mw_ptr, 1, Ainv_mw_ptr,
                                 lda, n_accepted);
    }
  }

public:
  DelayedUpdateBatched() : blas_handle(queue) {}

  DelayedUpdateBatched(const DelayedUpdateBatched&) = delete;

  /** resize the internal storage
   * @param nw number of walkers in the crowd
   * @param norb number of electrons/orbitals
   * @param delay, maximum delay 0<delay<=norb
   */
  void resize(size_t nw, size_t norb, size_t delay)
  {
    max_delay_         = delay;
    no_delayed_update_ = (delay == 1);
    delay_count        = 0;
    invRow_id          = -1;
    V_gpu.resize(nw, delay * norb);
    U_gpu.resize(nw, delay * norb);
    p_gpu.resize(nw, delay);
    tempMat_gpu.resize(nw, norb * delay);
    Binv_gpu.resize(nw, delay * delay);
    delay_list_gpu.resize(nw, delay);
    invRow.resize(nw, norb);
  }

  // prepare invRow and compute the old gradients.
  template<typename GT>
  void mw_evalGrad(const RefVector<DualMatrix<Value>>& psiMinv_refs,
                   const std::vector<const Value*>& dpsiM_row_list,
                   const int rowchanged,
                   std::vector<GT>& grad_now)
  {
    const int nw   = invRow.rows();
    const int norb = invRow.cols();
    assert(nw == psiMinv_refs.size());
    assert(nw == dpsiM_row_list.size());
    assert(nw == grad_now.size());
    if (!no_delayed_update_)
      mw_prepareInvRow(psiMinv_refs, rowchanged);

    constexpr size_t num_ptrs_packed = 2; // it must match packing and unpacking
    evalGrad_buffer_H2D.resize(sizeof(Value*) * num_ptrs_packed * nw);
    Matrix<const Value*> ptr_buffer(reinterpret_cast<const Value**>(evalGrad_buffer_H2D.data()), num_ptrs_packed, nw);
    for (int iw = 0; iw < nw; iw++)
    {
      if (no_delayed_update_)
      {
        DualMatrix<Value>& psiMinv = psiMinv_refs[iw];
        ptr_buffer[0][iw]          = psiMinv.device_data() + rowchanged * psiMinv.cols();
      }
      else
        ptr_buffer[0][iw] = invRow_at(iw);
      ptr_buffer[1][iw] = dpsiM_row_list[iw];
    }

    queue.enqueueH2D(evalGrad_buffer_H2D);

    if (grads_value_v.rows() != nw || grads_value_v.cols() != GT::Size)
      grads_value_v.resize(nw, GT::Size);

    const Value** invRow_ptr    = devPtr<const Value*>(evalGrad_buffer_H2D, 0, nw);
    const Value** dpsiM_row_ptr = devPtr<const Value*>(evalGrad_buffer_H2D, 1, nw);

    compute::calcGradients_batched(queue, norb, invRow_ptr, dpsiM_row_ptr, grads_value_v.device_data(), nw);
    queue.enqueueD2H(grads_value_v);
    queue.sync();

    for (int iw = 0; iw < nw; iw++)
      grad_now[iw] = {grads_value_v[iw][0], grads_value_v[iw][1], grads_value_v[iw][2]};
  }

  template<typename GT>
  void mw_evalGradWithSpin(const RefVector<DualMatrix<Value>>& psiMinv_refs,
                           const std::vector<const Value*>& dpsiM_row_list,
                           OffloadMatrix<Complex>& mw_dspin,
                           const int rowchanged,
                           std::vector<GT>& grad_now,
                           std::vector<Complex>& spingrad_now)
  {
    const int nw     = invRow.rows();
    const int norb   = invRow.cols();
    assert(nw == psiMinv_refs.size());
    assert(nw == dpsiM_row_list.size());
    assert(nw == mw_dspin.rows());
    assert(nw == grad_now.size());
    assert(nw == spingrad_now.size());
    auto& buffer_H2D = evalGrad_buffer_H2D;

    //Need to pack these into a transfer buffer since psiMinv and dpsiM_row_list are not multiwalker data
    //i.e. each engine has its own psiMinv which is an OffloadMatrix instead of the leader having the data for all the walkers in the crowd.
    //Wouldn't have to do this if dpsiM and psiMinv were part of the mw_rsc_handle_ with data across all walkers in the crowd and could just use use_device_ptr for the offload.
    //That is how mw_dspin is handled below
    constexpr size_t num_ptrs_packed = 2; // it must match packing and unpacking
    buffer_H2D.resize(sizeof(Value*) * num_ptrs_packed * nw);
    Matrix<const Value*> ptr_buffer(reinterpret_cast<const Value**>(buffer_H2D.data()), num_ptrs_packed, nw);
    for (int iw = 0; iw < nw; iw++)
    {
      DualMatrix<Value>& psiMinv = psiMinv_refs[iw];
      ptr_buffer[0][iw]          = psiMinv.device_data() + rowchanged * psiMinv.cols();
      ptr_buffer[1][iw]          = dpsiM_row_list[iw];
    }
    constexpr unsigned DIM = GT::Size;
    grads_value_v.resize(nw, DIM);
    spingrads_value_v.resize(nw);
    auto* __restrict__ grads_value_v_ptr     = grads_value_v.data();
    auto* __restrict__ spingrads_value_v_ptr = spingrads_value_v.data();
    auto* buffer_H2D_ptr                     = buffer_H2D.data();
    auto* mw_dspin_ptr                       = mw_dspin.data();

    //Note that mw_dspin should already be in sync between device and host...updateTo was called in
    //SPOSet::mw_evaluateVGLWithSpin to sync
    //Also note that since mw_dspin is Dual, I can just use mw_dpsin.data() above and then use directly inside
    //then offload region. OMP will figure out the correct translation to the device address, i.e. no
    //need to include in the PRAGMA_OFFLOAD below
    PRAGMA_OFFLOAD("omp target teams distribute num_teams(nw) \
                    map(always, to: buffer_H2D_ptr[:buffer_H2D.size()]) \
                    map(always, from: grads_value_v_ptr[:grads_value_v.size()]) \
                    map(always, from: spingrads_value_v_ptr[:spingrads_value_v.size()])")
    for (int iw = 0; iw < nw; iw++)
    {
      const Value* __restrict__ invRow_ptr    = reinterpret_cast<const Value**>(buffer_H2D_ptr)[iw];
      const Value* __restrict__ dpsiM_row_ptr = reinterpret_cast<const Value**>(buffer_H2D_ptr)[nw + iw];
      Value grad_x(0), grad_y(0), grad_z(0);
      Complex spingrad(0);
#if defined(QMC_COMPLEX)
      // COMPILER WORKAROUND
      // This was causing a llvm-link error in icpx due to the lack of declare reduction on complex datatypes.
      // Keep real builds free of any reduction on a complex datatype. Just serialize the reduction.
      // Because mw_evalGradWithSpin is only being called in complex builds in simulations, the impact of this workaround is basically zero.
      // It is still beneficial to keep it functional in real builds.
      PRAGMA_OFFLOAD("omp parallel for reduction(+: grad_x, grad_y, grad_z, spingrad)")
#endif
      for (int iorb = 0; iorb < norb; iorb++)
      {
        grad_x += invRow_ptr[iorb] * dpsiM_row_ptr[iorb * DIM];
        grad_y += invRow_ptr[iorb] * dpsiM_row_ptr[iorb * DIM + 1];
        grad_z += invRow_ptr[iorb] * dpsiM_row_ptr[iorb * DIM + 2];
        spingrad += invRow_ptr[iorb] * mw_dspin_ptr[iw * norb + iorb];
      }
      grads_value_v_ptr[iw * DIM]     = grad_x;
      grads_value_v_ptr[iw * DIM + 1] = grad_y;
      grads_value_v_ptr[iw * DIM + 2] = grad_z;
      spingrads_value_v_ptr[iw]       = spingrad;
    }

    for (int iw = 0; iw < nw; iw++)
    {
      grad_now[iw]     = {grads_value_v[iw][0], grads_value_v[iw][1], grads_value_v[iw][2]};
      spingrad_now[iw] = spingrads_value_v[iw];
    }
  }

  /** Accept or Reject row updates
   *  many of these const arguments provide pointers or references
   *  to objects that do get modified.
   *  \param[in] psiMinv_refs
   *  \param[in] rowchanged
   *  \param[in] psiM_g_list
   *  \param[in] psiM_l_list
   *  \param[in] isAccepted
   *  \param[in] phi_vgl_v          multiple walker orbital VGL
   *  \param[inout] ratios
   */
  void mw_accept_rejectRow(const RefVector<DualMatrix<Value>>& psiMinv_refs,
                           const int rowchanged,
                           const std::vector<Value*>& psiM_g_list,
                           const std::vector<Value*>& psiM_l_list,
                           const std::vector<bool>& isAccepted,
                           const OffloadMWVGLArray<Value>& phi_vgl_v,
                           const std::vector<Value>& ratios)
  {
    const int nw   = invRow.rows();
    const int norb = invRow.cols();
    assert(nw == psiMinv_refs.size());
    assert(nw == isAccepted.size());
    assert(nw == ratios.size());
    assert(nw == phi_vgl_v.size(1));
    // invRow consumed, mark invRow_id unset
    invRow_id = -1;

    if (no_delayed_update_)
    {
      mw_updateRow(psiMinv_refs, rowchanged, psiM_g_list, psiM_l_list, isAccepted, phi_vgl_v, ratios);
      return;
    }

    const int lda_Binv          = max_delay_;
    const int n_accepted        = psiM_g_list.size();
    const size_t phi_vgl_stride = nw * norb;

    constexpr size_t num_ptrs_packed = 12; // it must match packing and unpacking
    accept_rejectRow_buffer_H2D.resize((sizeof(Value*) * num_ptrs_packed + sizeof(Value)) * nw);
    resize_fill_constant_arrays(nw);

    Matrix<Value*> ptr_buffer(reinterpret_cast<Value**>(accept_rejectRow_buffer_H2D.data()), num_ptrs_packed, nw);
    Value* c_ratio_inv =
        reinterpret_cast<Value*>(accept_rejectRow_buffer_H2D.data() + sizeof(Value*) * num_ptrs_packed * nw);
    for (int iw = 0, count_accepted = 0, count_rejected = 0; iw < nw; iw++)
    {
      DualMatrix<Value>& psiMinv = psiMinv_refs[iw];
      const int lda              = psiMinv.cols();
      Value* const V_iw          = V_gpu[iw];
      Value* const U_iw          = U_gpu[iw];
      Value* const p_iw          = p_gpu[iw];
      Value* const Binv_iw       = Binv_gpu[iw];
      int* const delay_list_iw   = delay_list_gpu[iw];
      if (isAccepted[iw])
      {
        ptr_buffer[0][count_accepted]  = psiMinv.device_data() + lda * rowchanged;
        ptr_buffer[1][count_accepted]  = V_iw;
        ptr_buffer[2][count_accepted]  = U_iw + norb * delay_count;
        ptr_buffer[3][count_accepted]  = p_iw;
        ptr_buffer[4][count_accepted]  = Binv_iw;
        ptr_buffer[5][count_accepted]  = Binv_iw + delay_count * lda_Binv;
        ptr_buffer[6][count_accepted]  = Binv_iw + delay_count;
        ptr_buffer[7][count_accepted]  = reinterpret_cast<Value*>(delay_list_iw);
        ptr_buffer[8][count_accepted]  = V_iw + norb * delay_count;
        ptr_buffer[9][count_accepted]  = const_cast<Value*>(phi_vgl_v.device_data_at(0, iw, 0));
        ptr_buffer[10][count_accepted] = psiM_g_list[count_accepted];
        ptr_buffer[11][count_accepted] = psiM_l_list[count_accepted];
        c_ratio_inv[count_accepted]    = Value(1) / ratios[iw];
        count_accepted++;
      }
      else
      {
        ptr_buffer[0][n_accepted + count_rejected] = psiMinv.device_data() + lda * rowchanged;
        ptr_buffer[1][n_accepted + count_rejected] = V_iw;
        ptr_buffer[2][n_accepted + count_rejected] = U_iw + norb * delay_count;
        ptr_buffer[3][n_accepted + count_rejected] = p_iw;
        ptr_buffer[4][n_accepted + count_rejected] = Binv_iw;
        ptr_buffer[5][n_accepted + count_rejected] = Binv_iw + delay_count * lda_Binv;
        ptr_buffer[6][n_accepted + count_rejected] = Binv_iw + delay_count;
        ptr_buffer[7][n_accepted + count_rejected] = reinterpret_cast<Value*>(delay_list_iw);
        ptr_buffer[8][n_accepted + count_rejected] = V_iw + norb * delay_count;
        count_rejected++;
      }
    }

    queue.enqueueH2D(accept_rejectRow_buffer_H2D);

    Value** invRow_mw_ptr   = devPtr<Value*>(accept_rejectRow_buffer_H2D, 0, nw);
    Value** V_mw_ptr        = devPtr<Value*>(accept_rejectRow_buffer_H2D, 1, nw);
    Value** U_row_mw_ptr    = devPtr<Value*>(accept_rejectRow_buffer_H2D, 2, nw);
    Value** p_mw_ptr        = devPtr<Value*>(accept_rejectRow_buffer_H2D, 3, nw);
    Value** Binv_mw_ptr     = devPtr<Value*>(accept_rejectRow_buffer_H2D, 4, nw);
    Value** BinvRow_mw_ptr  = devPtr<Value*>(accept_rejectRow_buffer_H2D, 5, nw);
    Value** BinvCol_mw_ptr  = devPtr<Value*>(accept_rejectRow_buffer_H2D, 6, nw);
    int** delay_list_mw_ptr = devPtr<int*>(accept_rejectRow_buffer_H2D, 7, nw);
    Value** V_row_mw_ptr    = devPtr<Value*>(accept_rejectRow_buffer_H2D, 8, nw);
    Value** phiVGL_mw_ptr   = devPtr<Value*>(accept_rejectRow_buffer_H2D, 9, nw);
    Value** dpsiM_mw_out    = devPtr<Value*>(accept_rejectRow_buffer_H2D, 10, nw);
    Value** d2psiM_mw_out   = devPtr<Value*>(accept_rejectRow_buffer_H2D, 11, nw);
    Value* ratio_inv_mw_ptr = devPtr<Value>(accept_rejectRow_buffer_H2D, 12, nw);

    //std::copy_n(Ainv[rowchanged], norb, V[delay_count]);
    compute::BLAS::copy_batched(blas_handle, norb, invRow_mw_ptr, 1, V_row_mw_ptr, 1, nw);
    // handle accepted walkers
    // the new Binv is [[X y] [z sigma]]
    //BLAS::gemv('T', norb, delay_count + 1, cminusone, V.data(), norb, psiV.data(), 1, czero, p.data(), 1);
    compute::BLAS::gemv_batched(blas_handle, 'T', norb, delay_count, cminusone_vec.device_data(), V_mw_ptr, norb,
                                phiVGL_mw_ptr, 1, czero_vec.device_data(), p_mw_ptr, 1, n_accepted);
    // y
    //BLAS::gemv('T', delay_count, delay_count, sigma, Binv.data(), lda_Binv, p.data(), 1, czero, Binv.data() + delay_count,
    //           lda_Binv);
    compute::BLAS::gemv_batched(blas_handle, 'T', delay_count, delay_count, ratio_inv_mw_ptr, Binv_mw_ptr, lda_Binv,
                                p_mw_ptr, 1, czero_vec.device_data(), BinvCol_mw_ptr, lda_Binv, n_accepted);
    // X
    //BLAS::ger(delay_count, delay_count, cone, Binv[delay_count], 1, Binv.data() + delay_count, lda_Binv,
    //          Binv.data(), lda_Binv);
    compute::BLAS::ger_batched(blas_handle, delay_count, delay_count, cone_vec.device_data(), BinvRow_mw_ptr, 1,
                               BinvCol_mw_ptr, lda_Binv, Binv_mw_ptr, lda_Binv, n_accepted);
    // sigma and Z
    compute::add_delay_list_save_sigma_VGL_batched(queue, delay_list_mw_ptr, rowchanged, delay_count, Binv_mw_ptr,
                                                   lda_Binv, ratio_inv_mw_ptr, phiVGL_mw_ptr, phi_vgl_stride,
                                                   U_row_mw_ptr, dpsiM_mw_out, d2psiM_mw_out, norb, n_accepted, nw);
    delay_count++;
    // update Ainv when maximal delay is reached
    if (delay_count == lda_Binv)
      mw_updateInvMat(psiMinv_refs);
  }

  /** update the full Ainv and reset delay_count
   * @param psiMinv_refs inverse matrices of all the walkers
   */
  void mw_updateInvMat(const RefVector<DualMatrix<Value>>& psiMinv_refs)
  {
    const int nw   = invRow.rows();
    const int norb = invRow.cols();
    assert(nw == psiMinv_refs.size());
    if (delay_count == 0)
      return;
    // update the inverse matrix
    const int lda = psiMinv_refs[0].get().cols();

    constexpr size_t num_ptrs_packed = 6; // it must match packing and unpacking
    updateInv_buffer_H2D.resize(sizeof(Value*) * num_ptrs_packed * nw);
    resize_fill_constant_arrays(nw);

    Matrix<Value*> ptr_buffer(reinterpret_cast<Value**>(updateInv_buffer_H2D.data()), num_ptrs_packed, nw);
    for (int iw = 0; iw < nw; iw++)
    {
      ptr_buffer[0][iw] = U_gpu[iw];
      ptr_buffer[1][iw] = psiMinv_refs[iw].get().device_data();
      ptr_buffer[2][iw] = tempMat_gpu[iw];
      ptr_buffer[3][iw] = reinterpret_cast<Value*>(delay_list_gpu[iw]);
      ptr_buffer[4][iw] = V_gpu[iw];
      ptr_buffer[5][iw] = Binv_gpu[iw];
    }

    queue.enqueueH2D(updateInv_buffer_H2D);

    Value** U_mw_ptr        = devPtr<Value*>(updateInv_buffer_H2D, 0, nw);
    Value** Ainv_mw_ptr     = devPtr<Value*>(updateInv_buffer_H2D, 1, nw);
    Value** tempMat_mw_ptr  = devPtr<Value*>(updateInv_buffer_H2D, 2, nw);
    int** delay_list_mw_ptr = devPtr<int*>(updateInv_buffer_H2D, 3, nw);
    Value** V_mw_ptr        = devPtr<Value*>(updateInv_buffer_H2D, 4, nw);
    Value** Binv_mw_ptr     = devPtr<Value*>(updateInv_buffer_H2D, 5, nw);

    /*
    if (delay_count == 1)
    {
      // this is a special case invoking the Fahy's variant of Sherman-Morrison update.
      // Only use the first norb elements of tempMat as a temporal array
      BLAS::gemv('T', norb, norb, cone, Ainv.data(), norb, U[0], 1, czero, temp.data(), 1);
      temp[delay_list[0]] -= cone;
      BLAS::ger(norb, norb, -Binv[0][0], V[0], 1, temp.data(), 1, Ainv.data(), norb);
    }
    else
*/
    {
      const int lda_Binv = max_delay_;
      compute::BLAS::gemm_batched(blas_handle, 'T', 'N', delay_count, norb, norb, Value(1), U_mw_ptr, norb, Ainv_mw_ptr,
                                  lda, Value(0), tempMat_mw_ptr, lda_Binv, nw);
      compute::applyW_batched(queue, delay_list_mw_ptr, delay_count, tempMat_mw_ptr, lda_Binv, nw);
      compute::BLAS::gemm_batched(blas_handle, 'N', 'N', norb, delay_count, delay_count, Value(1), V_mw_ptr, norb,
                                  Binv_mw_ptr, lda_Binv, Value(0), U_mw_ptr, norb, nw);
      compute::BLAS::gemm_batched(blas_handle, 'N', 'N', norb, norb, delay_count, Value(-1), U_mw_ptr, norb,
                                  tempMat_mw_ptr, lda_Binv, Value(1), Ainv_mw_ptr, lda, nw);
    }
    delay_count = 0;
  }

  /** return invRow host or device pointers based on on_host request
   * prepare invRow if not already.
   */
  std::vector<const Value*> mw_getInvRow(const RefVector<DualMatrix<Value>>& psiMinv_refs,
                                         const int row_id,
                                         bool on_host)
  {
    const size_t nw = invRow.rows();
    assert(nw == psiMinv_refs.size());
    if (no_delayed_update_)
      queue.sync();
    else if (invRow_id != row_id)
    {
      // this can be skipped if mw_evalGrad gets called already.
      mw_prepareInvRow(psiMinv_refs, row_id);
      queue.sync();
    }

    std::vector<const Value*> row_ptr_list;
    row_ptr_list.reserve(psiMinv_refs.size());
    if (on_host)
    {
      // copy values to host and return host pointer
      if (no_delayed_update_)
        for (DualMatrix<Value>& psiMinv : psiMinv_refs)
        {
          const size_t ncols = psiMinv.cols();
          psiMinv.updateFrom(ncols, row_id * ncols);
          row_ptr_list.push_back(psiMinv.data() + row_id * ncols);
        }
      else
      {
        invRow.updateFrom();
        for (size_t iw = 0; iw < nw; iw++)
          row_ptr_list.push_back(invRow[iw]);
      }
    }
    else
    {
      // return device pointer
      if (no_delayed_update_)
        for (DualMatrix<Value>& psiMinv : psiMinv_refs)
          row_ptr_list.push_back(psiMinv.device_data() + row_id * psiMinv.cols());
      else
        for (size_t iw = 0; iw < nw; iw++)
          row_ptr_list.push_back(invRow_at(iw));
    }
    return row_ptr_list;
  }

  /// transfer Ainv to the host
  void mw_transferAinv_D2H(const RefVector<DualMatrix<Value>>& psiMinv_refs)
  {
    assert(invRow.rows() == psiMinv_refs.size());
    guard_no_delay();

    for (DualMatrix<Value>& psiMinv : psiMinv_refs)
      queue.enqueueD2H(psiMinv);
    queue.sync();
  }
};
} // namespace qmcplusplus

#endif // QMCPLUSPLUS_DELAYED_UPDATE_BATCHED_H
