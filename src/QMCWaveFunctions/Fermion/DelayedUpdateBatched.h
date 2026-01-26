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

#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "QMCWaveFunctions/Fermion/DiracMatrix.h"
#include "Platforms/OMPTarget/ompBLAS.hpp"
#include "detail/AccelMatrixUpdate.hpp"
#include "WaveFunctionTypes.hpp"
#include "QueueAliases.hpp"
#include "AccelBLAS.hpp"

namespace qmcplusplus
{

/** implements dirac matrix delayed update using OpenMP offload and CUDA.
 * It is used as DET_ENGINE in DiracDeterminantBatched.
 * This is a 1 per walker class
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
  using DualVector = Vector<DT, OffloadPinnedAllocator<DT>>;
  template<typename DT>
  using DualMatrix = Matrix<DT, OffloadPinnedAllocator<DT>>;
  template<typename DT>
  using DualVGLVector = VectorSoaContainer<DT, QMCTraits::DIM + 2, OffloadPinnedAllocator<DT>>;
  template<typename DT>
  using OffloadMWVGLArray = Array<DT, 3, OffloadPinnedAllocator<DT>>; // [VGL, walker, Orbs]
  template<typename DT>
  using OffloadMatrix = Matrix<DT, OffloadPinnedAllocator<DT>>;

  struct MultiWalkerResource
  {
    // CUDA stream, cublas handle object
    compute::Queue<PL> queue;
    compute::BLASHandle<PL> blas_handle;

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
    /// scratch space for rank-1 update
    UnpinnedDualVector<Value> mw_temp;
    // scratch space for keeping one row of Ainv
    UnpinnedDualVector<Value> mw_rcopy;

    MultiWalkerResource() : blas_handle(queue) {}

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
  };

private:
  /// scratch space for rank-1 update
  UnpinnedDualVector<Value> temp;
  /// row of up-to-date Ainv
  UnpinnedDualVector<Value> invRow;
  /** row id correspond to the up-to-date invRow. [0 norb), invRow is ready; -1, invRow is not valid.
   *  This id is set after calling getInvRow indicating invRow has been prepared for the invRow_id row
   *  ratioGrad checks if invRow_id is consistent. If not, invRow needs to be recomputed.
   *  acceptMove and completeUpdates mark invRow invalid by setting invRow_id to -1
   */
  int invRow_id;
  // scratch space for keeping one row of Ainv
  UnpinnedDualVector<Value> rcopy;

  template<typename DT>
  using DeviceMatrix = Matrix<DT, OffloadDeviceAllocator<DT>>;
  template<typename DT>
  using DeviceVector = Vector<DT, OffloadDeviceAllocator<DT>>;
  /// orbital values of delayed electrons
  DeviceMatrix<Value> U_gpu;
  /// rows of Ainv corresponding to delayed electrons
  DeviceMatrix<Value> V_gpu;
  /// Matrix inverse of B, at maximum KxK
  DeviceMatrix<Value> Binv_gpu;
  /// scratch space, used during inverse update
  DeviceMatrix<Value> tempMat_gpu;
  /// new column of B
  DeviceVector<Value> p_gpu;
  /// list of delayed electrons
  DeviceVector<int> delay_list_gpu;
  /// current number of delays, increase one for each acceptance, reset to 0 after updating Ainv
  int delay_count;
  /// if true, updates are not delayed.
  const bool no_delayed_update_;

  /** resize the internal storage
   * @param norb number of electrons/orbitals
   * @param delay, maximum delay 0<delay<=norb
   */
  inline void resize(int norb, int delay)
  {
    V_gpu.resize(delay, norb);
    U_gpu.resize(delay, norb);
    p_gpu.resize(delay);
    tempMat_gpu.resize(norb, delay);
    Binv_gpu.resize(delay, delay);
    delay_list_gpu.resize(delay);
    invRow.resize(norb);
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
   * @param Ainv inverse matrix
   * @param rowchanged the row id corresponding to the proposed electron
   */
  static void mw_prepareInvRow(const RefVectorWithLeader<This_t>& engines,
                               MultiWalkerResource& mw_rsc,
                               const RefVector<DualMatrix<Value>>& psiMinv_refs,
                               const int rowchanged)
  {
    auto& engine_leader              = engines.getLeader();
    auto& blas_handle                = mw_rsc.blas_handle;
    auto& queue                      = mw_rsc.queue;
    auto& cminusone_vec              = mw_rsc.cminusone_vec;
    auto& cone_vec                   = mw_rsc.cone_vec;
    auto& czero_vec                  = mw_rsc.czero_vec;
    auto& prepare_inv_row_buffer_H2D = mw_rsc.prepare_inv_row_buffer_H2D;
    const int norb                   = engine_leader.invRow.size();
    const int nw                     = engines.size();
    int& delay_count                 = engine_leader.delay_count;

    constexpr size_t num_ptrs_packed = 7; // it must match packing and unpacking
    prepare_inv_row_buffer_H2D.resize(sizeof(Value*) * num_ptrs_packed * nw);
    mw_rsc.resize_fill_constant_arrays(nw);

    const int lda_Binv = engine_leader.Binv_gpu.cols();
    Matrix<Value*> ptr_buffer(reinterpret_cast<Value**>(prepare_inv_row_buffer_H2D.data()), num_ptrs_packed, nw);
    for (int iw = 0; iw < nw; iw++)
    {
      This_t& engine             = engines[iw];
      DualMatrix<Value>& psiMinv = psiMinv_refs[iw];
      ptr_buffer[0][iw]          = psiMinv.device_data() + rowchanged * psiMinv.cols();
      ptr_buffer[1][iw]          = engine.invRow.device_data();
      ptr_buffer[2][iw]          = engine.U_gpu.data();
      ptr_buffer[3][iw]          = engine.p_gpu.data();
      ptr_buffer[4][iw]          = engine.Binv_gpu.data();
      ptr_buffer[5][iw]          = engine.Binv_gpu.data() + delay_count * lda_Binv;
      ptr_buffer[6][iw]          = engine.V_gpu.data();
    }

    queue.enqueueH2D(prepare_inv_row_buffer_H2D);

    Value** oldRow_mw_ptr = reinterpret_cast<Value**>(prepare_inv_row_buffer_H2D.device_data());
    Value** invRow_mw_ptr = reinterpret_cast<Value**>(prepare_inv_row_buffer_H2D.device_data() + sizeof(Value*) * nw);
    Value** U_mw_ptr    = reinterpret_cast<Value**>(prepare_inv_row_buffer_H2D.device_data() + sizeof(Value*) * nw * 2);
    Value** p_mw_ptr    = reinterpret_cast<Value**>(prepare_inv_row_buffer_H2D.device_data() + sizeof(Value*) * nw * 3);
    Value** Binv_mw_ptr = reinterpret_cast<Value**>(prepare_inv_row_buffer_H2D.device_data() + sizeof(Value*) * nw * 4);
    Value** BinvRow_mw_ptr =
        reinterpret_cast<Value**>(prepare_inv_row_buffer_H2D.device_data() + sizeof(Value*) * nw * 5);
    Value** V_mw_ptr = reinterpret_cast<Value**>(prepare_inv_row_buffer_H2D.device_data() + sizeof(Value*) * nw * 6);

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
    engine_leader.invRow_id = rowchanged;
  }

  /** Do complete row updates
   *  many of these const arguments provide pointers or references
   *  somewhere in here is an update that doesn't get where it belongs resulting in a 0
   *  gradient later.
   *  Sad example of OpenMP target code that is far from clear and a poor substitute for a
   *  clear CPU reference implementation.
   *
   *  \param[in] engines
   *  \param[in] rowchanged
   *  \param[in] psiM_g_list        device ptrs
   *  \param[in] psiM_l_list        device ptrs
   *  \param[in] isAccepted         bool but wait some lists are also filtered
   *  \param[in] phi_vgl_v          multiple walker orbital VGL
   *  \param[inout] ratios
   */
  static void mw_updateRow(const RefVectorWithLeader<This_t>& engines,
                           MultiWalkerResource& mw_rsc,
                           const RefVector<DualMatrix<Value>>& psiMinv_refs,
                           const int rowchanged,
                           const std::vector<Value*>& psiM_g_list,
                           const std::vector<Value*>& psiM_l_list,
                           const std::vector<bool>& isAccepted,
                           const OffloadMWVGLArray<Value>& phi_vgl_v,
                           const std::vector<Value>& ratios)
  {
    auto& engine_leader = engines.getLeader();
    engine_leader.guard_no_delay();

    const size_t n_accepted = psiM_g_list.size();
#ifndef NDEBUG
    size_t n_true = std::count_if(isAccepted.begin(), isAccepted.end(), [](bool accepted) { return accepted; });
    assert(n_accepted == n_true);
#endif
    if (n_accepted == 0)
      return;

    auto& queue                 = mw_rsc.queue;
    auto& blas_handle           = mw_rsc.blas_handle;
    auto& updateRow_buffer_H2D  = mw_rsc.updateRow_buffer_H2D;
    auto& mw_temp               = mw_rsc.mw_temp;
    auto& mw_rcopy              = mw_rsc.mw_rcopy;
    auto& cone_vec              = mw_rsc.cone_vec;
    auto& czero_vec             = mw_rsc.czero_vec;
    const int norb              = engine_leader.invRow.size();
    const int lda               = psiMinv_refs[0].get().cols();
    const int nw                = engines.size();
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
    mw_rsc.resize_fill_constant_arrays(n_accepted);

    queue.enqueueH2D(updateRow_buffer_H2D);

    {
      Value** Ainv_mw_ptr = reinterpret_cast<Value**>(updateRow_buffer_H2D.device_data());
      Value** phiVGL_mw_ptr =
          reinterpret_cast<Value**>(updateRow_buffer_H2D.device_data() + sizeof(Value*) * n_accepted);
      Value** temp_mw_ptr =
          reinterpret_cast<Value**>(updateRow_buffer_H2D.device_data() + sizeof(Value*) * n_accepted * 2);
      Value** rcopy_mw_ptr =
          reinterpret_cast<Value**>(updateRow_buffer_H2D.device_data() + sizeof(Value*) * n_accepted * 3);
      Value** dpsiM_mw_out =
          reinterpret_cast<Value**>(updateRow_buffer_H2D.device_data() + sizeof(Value*) * n_accepted * 4);
      Value** d2psiM_mw_out =
          reinterpret_cast<Value**>(updateRow_buffer_H2D.device_data() + sizeof(Value*) * n_accepted * 5);
      Value* ratio_inv_mw =
          reinterpret_cast<Value*>(updateRow_buffer_H2D.device_data() + sizeof(Value*) * n_accepted * 6);


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
  /// default constructor
  DelayedUpdateBatched(size_t norb, size_t max_delay)
      : invRow_id(-1), delay_count(0), no_delayed_update_(max_delay == 1)
  {
    resize(norb, max_delay);
  }

  DelayedUpdateBatched(const DelayedUpdateBatched&) = delete;

  // prepare invRow and compute the old gradients.
  template<typename GT>
  static void mw_evalGrad(const RefVectorWithLeader<This_t>& engines,
                          MultiWalkerResource& mw_rsc,
                          const RefVector<DualMatrix<Value>>& psiMinv_refs,
                          const std::vector<const Value*>& dpsiM_row_list,
                          const int rowchanged,
                          std::vector<GT>& grad_now)
  {
    auto& engine_leader = engines.getLeader();
    if (!engine_leader.no_delayed_update_)
      mw_prepareInvRow(engines, mw_rsc, psiMinv_refs, rowchanged);

    auto& queue               = mw_rsc.queue;
    auto& evalGrad_buffer_H2D = mw_rsc.evalGrad_buffer_H2D;
    auto& grads_value_v       = mw_rsc.grads_value_v;

    const int nw                     = engines.size();
    constexpr size_t num_ptrs_packed = 2; // it must match packing and unpacking
    evalGrad_buffer_H2D.resize(sizeof(Value*) * num_ptrs_packed * nw);
    Matrix<const Value*> ptr_buffer(reinterpret_cast<const Value**>(evalGrad_buffer_H2D.data()), num_ptrs_packed, nw);
    for (int iw = 0; iw < nw; iw++)
    {
      if (engine_leader.no_delayed_update_)
      {
        DualMatrix<Value>& psiMinv = psiMinv_refs[iw];
        ptr_buffer[0][iw]          = psiMinv.device_data() + rowchanged * psiMinv.cols();
      }
      else
        ptr_buffer[0][iw] = engines[iw].invRow.device_data();
      ptr_buffer[1][iw] = dpsiM_row_list[iw];
    }

    queue.enqueueH2D(evalGrad_buffer_H2D);

    if (grads_value_v.rows() != nw || grads_value_v.cols() != GT::Size)
      grads_value_v.resize(nw, GT::Size);

    const Value** invRow_ptr    = reinterpret_cast<const Value**>(evalGrad_buffer_H2D.device_data());
    const Value** dpsiM_row_ptr = reinterpret_cast<const Value**>(evalGrad_buffer_H2D.device_data()) + nw;

    compute::calcGradients_batched(queue, engine_leader.invRow.size(), invRow_ptr, dpsiM_row_ptr,
                                   grads_value_v.device_data(), nw);
    queue.enqueueD2H(grads_value_v);
    queue.sync();

    for (int iw = 0; iw < nw; iw++)
      grad_now[iw] = {grads_value_v[iw][0], grads_value_v[iw][1], grads_value_v[iw][2]};
  }

  template<typename GT>
  static void mw_evalGradWithSpin(const RefVectorWithLeader<This_t>& engines,
                                  MultiWalkerResource& mw_rsc,
                                  const RefVector<DualMatrix<Value>>& psiMinv_refs,
                                  const std::vector<const Value*>& dpsiM_row_list,
                                  OffloadMatrix<Complex>& mw_dspin,
                                  const int rowchanged,
                                  std::vector<GT>& grad_now,
                                  std::vector<Complex>& spingrad_now)
  {
    auto& engine_leader     = engines.getLeader();
    auto& buffer_H2D        = mw_rsc.evalGrad_buffer_H2D;
    auto& grads_value_v     = mw_rsc.grads_value_v;
    auto& spingrads_value_v = mw_rsc.spingrads_value_v;

    //Need to pack these into a transfer buffer since psiMinv and dpsiM_row_list are not multiwalker data
    //i.e. each engine has its own psiMinv which is an OffloadMatrix instead of the leader having the data for all the walkers in the crowd.
    //Wouldn't have to do this if dpsiM and psiMinv were part of the mw_rsc_handle_ with data across all walkers in the crowd and could just use use_device_ptr for the offload.
    //That is how mw_dspin is handled below
    const int norb                   = psiMinv_refs[0].get().rows();
    const int nw                     = engines.size();
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

  /** Update the "local" psiMinv_ on the device.
   *  Side Effect Transfers:
   *  * psiMinv_ is transferred back to host since single calls from QMCHamitonian and others
   *  * expect it to be.
   *
   *  Forced to use OpenMP target since resources are banned for single walker functions APIs
   *  and the acquireRelease pattern for a single DDB was removed by #3324
   */
  template<typename VVT, typename FPVT>
  void updateRow(DualMatrix<Value>& Ainv, int rowchanged, const VVT& phiV, FPVT c_ratio_in)
  {
    guard_no_delay();
    // update the inverse matrix
    constexpr Value cone(1), czero(0);
    const int norb = Ainv.rows();
    const int lda  = Ainv.cols();
    temp.resize(norb);
    rcopy.resize(norb);
    // invoke the Fahy's variant of Sherman-Morrison update.
    int dummy_handle      = 0;
    const Value* phiV_ptr = phiV.data();
    Value* Ainv_ptr       = Ainv.data();
    Value* temp_ptr       = temp.data();
    Value* rcopy_ptr      = rcopy.data();
    // psiMinv_ must be update-to-date on both host and device
    PRAGMA_OFFLOAD("omp target data map(always, tofrom: Ainv_ptr[:Ainv.size()]) \
                    use_device_ptr(phiV_ptr, Ainv_ptr, temp_ptr, rcopy_ptr)")
    {
      int success = ompBLAS::gemv(dummy_handle, 'T', norb, norb, cone, Ainv_ptr, lda, phiV_ptr, 1, czero, temp_ptr, 1);
      if (success != 0)
        throw std::runtime_error("ompBLAS::gemv failed.");

      PRAGMA_OFFLOAD("omp target parallel for simd is_device_ptr(Ainv_ptr, temp_ptr, rcopy_ptr)")
      for (int i = 0; i < norb; i++)
      {
        rcopy_ptr[i] = Ainv_ptr[rowchanged * lda + i];
        if (i == 0)
          temp_ptr[rowchanged] -= cone;
      }

      success = ompBLAS::ger(dummy_handle, norb, norb, static_cast<Value>(FPVT(-1) / c_ratio_in), rcopy_ptr, 1,
                             temp_ptr, 1, Ainv_ptr, lda);
      if (success != 0)
        throw std::runtime_error("ompBLAS::ger failed.");
    }
  }

  /** Accept or Reject row updates
   *  many of these const arguments provide pointers or references
   *  to objects that do get modified.
   *  \param[in] engines
   *  \param[in] rowchanged
   *  \param[in] psiM_g_list
   *  \param[in] psiM_l_list
   *  \param[in] isAccepted
   *  \param[in] phi_vgl_v          multiple walker orbital VGL
   *  \param[inout] ratios
   */
  static void mw_accept_rejectRow(const RefVectorWithLeader<This_t>& engines,
                                  MultiWalkerResource& mw_rsc,
                                  const RefVector<DualMatrix<Value>>& psiMinv_refs,
                                  const int rowchanged,
                                  const std::vector<Value*>& psiM_g_list,
                                  const std::vector<Value*>& psiM_l_list,
                                  const std::vector<bool>& isAccepted,
                                  const OffloadMWVGLArray<Value>& phi_vgl_v,
                                  const std::vector<Value>& ratios)
  {
    auto& engine_leader = engines.getLeader();
    // invRow consumed, mark invRow_id unset
    engine_leader.invRow_id = -1;

    if (engine_leader.no_delayed_update_)
    {
      mw_updateRow(engines, mw_rsc, psiMinv_refs, rowchanged, psiM_g_list, psiM_l_list, isAccepted, phi_vgl_v, ratios);
      return;
    }

    auto& queue                       = mw_rsc.queue;
    auto& blas_handle                 = mw_rsc.blas_handle;
    auto& cminusone_vec               = mw_rsc.cminusone_vec;
    auto& cone_vec                    = mw_rsc.cone_vec;
    auto& czero_vec                   = mw_rsc.czero_vec;
    auto& accept_rejectRow_buffer_H2D = mw_rsc.accept_rejectRow_buffer_H2D;
    int& delay_count                  = engine_leader.delay_count;
    const int lda_Binv                = engine_leader.Binv_gpu.cols();
    const int norb                    = engine_leader.invRow.size();
    const int nw                      = engines.size();
    const int n_accepted              = psiM_g_list.size();
    const size_t phi_vgl_stride       = nw * norb;

    constexpr size_t num_ptrs_packed = 12; // it must match packing and unpacking
    accept_rejectRow_buffer_H2D.resize((sizeof(Value*) * num_ptrs_packed + sizeof(Value)) * nw);
    mw_rsc.resize_fill_constant_arrays(nw);

    Matrix<Value*> ptr_buffer(reinterpret_cast<Value**>(accept_rejectRow_buffer_H2D.data()), num_ptrs_packed, nw);
    Value* c_ratio_inv =
        reinterpret_cast<Value*>(accept_rejectRow_buffer_H2D.data() + sizeof(Value*) * num_ptrs_packed * nw);
    for (int iw = 0, count_accepted = 0, count_rejected = 0; iw < nw; iw++)
    {
      DualMatrix<Value>& psiMinv = psiMinv_refs[iw];
      const int lda              = psiMinv.cols();
      This_t& engine             = engines[iw];
      if (isAccepted[iw])
      {
        ptr_buffer[0][count_accepted]  = psiMinv.device_data() + lda * rowchanged;
        ptr_buffer[1][count_accepted]  = engine.V_gpu.data();
        ptr_buffer[2][count_accepted]  = engine.U_gpu.data() + norb * delay_count;
        ptr_buffer[3][count_accepted]  = engine.p_gpu.data();
        ptr_buffer[4][count_accepted]  = engine.Binv_gpu.data();
        ptr_buffer[5][count_accepted]  = engine.Binv_gpu.data() + delay_count * lda_Binv;
        ptr_buffer[6][count_accepted]  = engine.Binv_gpu.data() + delay_count;
        ptr_buffer[7][count_accepted]  = reinterpret_cast<Value*>(engine.delay_list_gpu.data());
        ptr_buffer[8][count_accepted]  = engine.V_gpu.data() + norb * delay_count;
        ptr_buffer[9][count_accepted]  = const_cast<Value*>(phi_vgl_v.device_data_at(0, iw, 0));
        ptr_buffer[10][count_accepted] = psiM_g_list[count_accepted];
        ptr_buffer[11][count_accepted] = psiM_l_list[count_accepted];
        c_ratio_inv[count_accepted]    = Value(1) / ratios[iw];
        count_accepted++;
      }
      else
      {
        ptr_buffer[0][n_accepted + count_rejected] = psiMinv.device_data() + lda * rowchanged;
        ptr_buffer[1][n_accepted + count_rejected] = engine.V_gpu.data();
        ptr_buffer[2][n_accepted + count_rejected] = engine.U_gpu.data() + norb * delay_count;
        ptr_buffer[3][n_accepted + count_rejected] = engine.p_gpu.data();
        ptr_buffer[4][n_accepted + count_rejected] = engine.Binv_gpu.data();
        ptr_buffer[5][n_accepted + count_rejected] = engine.Binv_gpu.data() + delay_count * lda_Binv;
        ptr_buffer[6][n_accepted + count_rejected] = engine.Binv_gpu.data() + delay_count;
        ptr_buffer[7][n_accepted + count_rejected] = reinterpret_cast<Value*>(engine.delay_list_gpu.data());
        ptr_buffer[8][n_accepted + count_rejected] = engine.V_gpu.data() + norb * delay_count;
        count_rejected++;
      }
    }

    queue.enqueueH2D(accept_rejectRow_buffer_H2D);

    Value** invRow_mw_ptr = reinterpret_cast<Value**>(accept_rejectRow_buffer_H2D.device_data());
    Value** V_mw_ptr      = reinterpret_cast<Value**>(accept_rejectRow_buffer_H2D.device_data() + sizeof(Value*) * nw);
    Value** U_row_mw_ptr =
        reinterpret_cast<Value**>(accept_rejectRow_buffer_H2D.device_data() + sizeof(Value*) * nw * 2);
    Value** p_mw_ptr = reinterpret_cast<Value**>(accept_rejectRow_buffer_H2D.device_data() + sizeof(Value*) * nw * 3);
    Value** Binv_mw_ptr =
        reinterpret_cast<Value**>(accept_rejectRow_buffer_H2D.device_data() + sizeof(Value*) * nw * 4);
    Value** BinvRow_mw_ptr =
        reinterpret_cast<Value**>(accept_rejectRow_buffer_H2D.device_data() + sizeof(Value*) * nw * 5);
    Value** BinvCol_mw_ptr =
        reinterpret_cast<Value**>(accept_rejectRow_buffer_H2D.device_data() + sizeof(Value*) * nw * 6);
    int** delay_list_mw_ptr =
        reinterpret_cast<int**>(accept_rejectRow_buffer_H2D.device_data() + sizeof(Value*) * nw * 7);
    Value** V_row_mw_ptr =
        reinterpret_cast<Value**>(accept_rejectRow_buffer_H2D.device_data() + sizeof(Value*) * nw * 8);
    Value** phiVGL_mw_ptr =
        reinterpret_cast<Value**>(accept_rejectRow_buffer_H2D.device_data() + sizeof(Value*) * nw * 9);
    Value** dpsiM_mw_out =
        reinterpret_cast<Value**>(accept_rejectRow_buffer_H2D.device_data() + sizeof(Value*) * nw * 10);
    Value** d2psiM_mw_out =
        reinterpret_cast<Value**>(accept_rejectRow_buffer_H2D.device_data() + sizeof(Value*) * nw * 11);
    Value* ratio_inv_mw_ptr =
        reinterpret_cast<Value*>(accept_rejectRow_buffer_H2D.device_data() + sizeof(Value*) * nw * 12);

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
      mw_updateInvMat(engines, mw_rsc, psiMinv_refs);
  }

  /** update the full Ainv and reset delay_count
   * @param Ainv inverse matrix
   */
  static void mw_updateInvMat(const RefVectorWithLeader<This_t>& engines,
                              MultiWalkerResource& mw_rsc,
                              const RefVector<DualMatrix<Value>>& psiMinv_refs)
  {
    auto& engine_leader = engines.getLeader();
    int& delay_count    = engine_leader.delay_count;
    if (delay_count == 0)
      return;
    // update the inverse matrix
    auto& queue                = mw_rsc.queue;
    auto& blas_handle          = mw_rsc.blas_handle;
    auto& updateInv_buffer_H2D = mw_rsc.updateInv_buffer_H2D;
    const int norb             = engine_leader.invRow.size();
    const int lda              = psiMinv_refs[0].get().cols();
    const int nw               = engines.size();

    constexpr size_t num_ptrs_packed = 6; // it must match packing and unpacking
    updateInv_buffer_H2D.resize(sizeof(Value*) * num_ptrs_packed * nw);
    mw_rsc.resize_fill_constant_arrays(nw);

    Matrix<Value*> ptr_buffer(reinterpret_cast<Value**>(updateInv_buffer_H2D.data()), num_ptrs_packed, nw);
    for (int iw = 0; iw < nw; iw++)
    {
      This_t& engine    = engines[iw];
      ptr_buffer[0][iw] = engine.U_gpu.data();
      ptr_buffer[1][iw] = psiMinv_refs[iw].get().device_data();
      ptr_buffer[2][iw] = engine.tempMat_gpu.data();
      ptr_buffer[3][iw] = reinterpret_cast<Value*>(engine.delay_list_gpu.data());
      ptr_buffer[4][iw] = engine.V_gpu.data();
      ptr_buffer[5][iw] = engine.Binv_gpu.data();
    }

    queue.enqueueH2D(updateInv_buffer_H2D);

    Value** U_mw_ptr        = reinterpret_cast<Value**>(updateInv_buffer_H2D.device_data());
    Value** Ainv_mw_ptr     = reinterpret_cast<Value**>(updateInv_buffer_H2D.device_data() + sizeof(Value*) * nw);
    Value** tempMat_mw_ptr  = reinterpret_cast<Value**>(updateInv_buffer_H2D.device_data() + sizeof(Value*) * nw * 2);
    int** delay_list_mw_ptr = reinterpret_cast<int**>(updateInv_buffer_H2D.device_data() + sizeof(Value*) * nw * 3);
    Value** V_mw_ptr        = reinterpret_cast<Value**>(updateInv_buffer_H2D.device_data() + sizeof(Value*) * nw * 4);
    Value** Binv_mw_ptr     = reinterpret_cast<Value**>(updateInv_buffer_H2D.device_data() + sizeof(Value*) * nw * 5);

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
      const int lda_Binv = engine_leader.Binv_gpu.cols();
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

  /*
  inline void print_Ainv(const RefVector<This_t>& engines)
  {
    for (This_t& engine : engines)
    {
      std::cout << "debug Ainv host  " << engine.get_psiMinv()[0][0] << " " << engine.get_psiMinv()[0][32] << " "
                << engine.get_psiMinv()[1][0] << " " << engine.get_psiMinv()[1][32] << std::endl;
      auto* temp_ptr = engine.get_psiMinv().data();
      PRAGMA_OFFLOAD("omp target update from(temp_ptr[:psiMinv_.size()])")
      std::cout << "debug Ainv devi";
      for (int row = 0; row < psiMinv_.rows(); row++)
      {
        for (int col = 0; col < psiMinv_.cols(); col++)
          std::cout << " " << row << " " << col << " " << engine.get_psiMinv()[row][col];
        std::cout << std::endl;
      }
    }
  }
*/

  /** return invRow host or device pointers based on on_host request
   * prepare invRow if not already.
   */
  static std::vector<const Value*> mw_getInvRow(const RefVectorWithLeader<This_t>& engines,
                                                MultiWalkerResource& mw_rsc,
                                                const RefVector<DualMatrix<Value>>& psiMinv_refs,
                                                const int row_id,
                                                bool on_host)
  {
    auto& engine_leader = engines.getLeader();
    auto& queue         = mw_rsc.queue;
    if (engine_leader.no_delayed_update_)
      queue.sync();
    else if (engine_leader.invRow_id != row_id)
    {
      // this can be skipped if mw_evalGrad gets called already.
      mw_prepareInvRow(engines, mw_rsc, psiMinv_refs, row_id);
      queue.sync();
    }

    std::vector<const Value*> row_ptr_list;
    row_ptr_list.reserve(psiMinv_refs.size());
    if (on_host)
    {
      // copy values to host and return host pointer
      if (engine_leader.no_delayed_update_)
        for (DualMatrix<Value>& psiMinv : psiMinv_refs)
        {
          const size_t ncols = psiMinv.cols();
          psiMinv.updateFrom(ncols, row_id * ncols);
          row_ptr_list.push_back(psiMinv.data() + row_id * ncols);
        }
      else
        for (This_t& engine : engines)
        {
          engine.invRow.updateFrom();
          row_ptr_list.push_back(engine.invRow.data());
        }
    }
    else
    {
      // return device pointer
      if (engine_leader.no_delayed_update_)
        for (DualMatrix<Value>& psiMinv : psiMinv_refs)
          row_ptr_list.push_back(psiMinv.device_data() + row_id * psiMinv.cols());
      else
        for (This_t& engine : engines)
          row_ptr_list.push_back(engine.invRow.device_data());
    }
    return row_ptr_list;
  }

  /// transfer Ainv to the host
  static void mw_transferAinv_D2H(const RefVectorWithLeader<This_t>& engines,
                                  MultiWalkerResource& mw_rsc,
                                  const RefVector<DualMatrix<Value>>& psiMinv_refs)
  {
    auto& engine_leader = engines.getLeader();
    auto& queue         = mw_rsc.queue;
    engine_leader.guard_no_delay();

    for (DualMatrix<Value>& psiMinv : psiMinv_refs)
      queue.enqueueD2H(psiMinv);
    queue.sync();
  }
};
} // namespace qmcplusplus

#endif // QMCPLUSPLUS_DELAYED_UPDATE_BATCHED_H
