//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_MATRIX_UPDATE_OMPTARGET_H
#define QMCPLUSPLUS_MATRIX_UPDATE_OMPTARGET_H

#include "OMPTarget/OffloadAlignedAllocators.hpp"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "OMPTarget/ompBLAS.hpp"
#include "OMPTarget/ompReductionComplex.hpp"
#include "ResourceCollection.h"
#include "DiracMatrixComputeOMPTarget.hpp"
#include "WaveFunctionTypes.hpp"

namespace qmcplusplus
{
/** Implements dirac matrix update using OpenMP offload.
 * It is used as DET_ENGINE in DiracDeterminantBatched.
 * @tparam VALUE base precision for most computation
 * @tparam VALUE_FP high precision for matrix inversion, T_FP >= T
 */
template<typename VALUE, typename VALUE_FP>
class MatrixUpdateOMPTarget
{
public:
  using WFT           = WaveFunctionTypes<VALUE, VALUE_FP>;
  using Value         = typename WFT::Value;
  using Complex       = typename WFT::Complex;
  using FullPrecValue = typename WFT::FullPrecValue;
  using LogValue      = typename WFT::LogValue;
  using This_t        = MatrixUpdateOMPTarget<VALUE, VALUE_FP>;
  using DetInverter   = DiracMatrixComputeOMPTarget<VALUE_FP>;

  // There is no expectation that this OMP only code can possibly work with anything but OMPallocator based
  // Allocators.  So this file should never include DualAllocatorAliases.hpp this violation of YAGNI isn't.
  // I do this for symmetry with MatrixDelayedUpdateCUDA only.
  template<typename DT>
  using PinnedDualAllocator = OffloadPinnedAllocator<DT>;

  template<typename DT>
  using UnpinnedOffloadVector = Vector<DT, OffloadAllocator<DT>>;
  template<typename DT>
  using OffloadVector = Vector<DT, OffloadPinnedAllocator<DT>>;
  template<typename DT>
  using OffloadMatrix = Matrix<DT, OffloadPinnedAllocator<DT>>;
  template<typename DT>
  using OffloadVGLVector = VectorSoaContainer<DT, QMCTraits::DIM + 2, OffloadPinnedAllocator<DT>>;
  template<typename DT>
  using OffloadMWVGLArray = Array<DT, 3, OffloadPinnedAllocator<DT>>; // [VGL, walker, Orbs]
  template<typename DT>
  using DualMatrix = Matrix<DT, PinnedDualAllocator<DT>>;

  struct MatrixUpdateOMPTargetMultiWalkerMem : public Resource
  {
    // constant array value T(1)
    UnpinnedOffloadVector<Value> cone_vec;
    // constant array value T(0)
    UnpinnedOffloadVector<Value> czero_vec;
    // multi walker of grads for transfer needs.
    OffloadMatrix<Value> grads_value_v;
    // multi walker of spingrads for transfer needs.
    OffloadVector<Complex> spingrads_value_v;
    // pointer buffer
    Vector<char, PinnedAllocator<char>> buffer_H2D;
    /// scratch space for rank-1 update
    UnpinnedOffloadVector<Value> mw_temp;
    // scratch space for keeping one row of Ainv
    UnpinnedOffloadVector<Value> mw_rcopy;


    MatrixUpdateOMPTargetMultiWalkerMem() : Resource("MatrixUpdateOMPTargetMultiWalkerMem") {}

    MatrixUpdateOMPTargetMultiWalkerMem(const MatrixUpdateOMPTargetMultiWalkerMem&)
        : MatrixUpdateOMPTargetMultiWalkerMem()
    {}

    std::unique_ptr<Resource> makeClone() const override
    {
      return std::make_unique<MatrixUpdateOMPTargetMultiWalkerMem>(*this);
    }
  };

  /// inverse transpose of psiM(j,i) \f$= \psi_j({\bf r}_i)\f$
  OffloadMatrix<Value> psiMinv_;
  /// scratch space for rank-1 update
  UnpinnedOffloadVector<Value> temp;
  // scratch space for keeping one row of Ainv
  UnpinnedOffloadVector<Value> rcopy;

  // multi walker memory buffers
  ResourceHandle<MatrixUpdateOMPTargetMultiWalkerMem> mw_mem_handle_;

  typename DetInverter::HandleResource dummy;

  void resize_fill_constant_arrays(size_t nw)
  {
    auto& mw_mem = mw_mem_handle_.getResource();
    if (mw_mem.cone_vec.size() < nw)
    {
      mw_mem.cone_vec.resize(nw);
      mw_mem.czero_vec.resize(nw);
      std::fill_n(mw_mem.cone_vec.data(), nw, Value(1));
      std::fill_n(mw_mem.czero_vec.data(), nw, Value(0));
      Value* cone_ptr = mw_mem.cone_vec.data();
      PRAGMA_OFFLOAD("omp target update to(cone_ptr[:nw])")
      Value* czero_ptr = mw_mem.czero_vec.data();
      PRAGMA_OFFLOAD("omp target update to(czero_ptr[:nw])")
    }
  }

  void resize_scratch_arrays(int norb, size_t nw)
  {
    size_t total_size = norb * nw;
    auto& mw_mem      = mw_mem_handle_.getResource();
    if (mw_mem.mw_temp.size() < total_size)
    {
      mw_mem.mw_temp.resize(total_size);
      mw_mem.mw_rcopy.resize(total_size);
    }
  }

public:
  /** resize the internal storage
   * @param norb number of electrons/orbitals
   * @param delay, maximum delay 0<delay<=norb
   */
  inline void resize(int norb, int delay) { psiMinv_.resize(norb, getAlignedSize<Value>(norb)); }

  void createResource(ResourceCollection& collection) const
  {
    collection.addResource(std::make_unique<MatrixUpdateOMPTargetMultiWalkerMem>());
  }

  void acquireResource(ResourceCollection& collection)
  {
    mw_mem_handle_ = collection.lendResource<MatrixUpdateOMPTargetMultiWalkerMem>();
  }

  void releaseResource(ResourceCollection& collection) { collection.takebackResource(mw_mem_handle_); }

  const OffloadMatrix<Value>& get_psiMinv() const { return psiMinv_; }
  OffloadMatrix<Value>& get_ref_psiMinv() { return psiMinv_; }

  inline Value* getRow_psiMinv_offload(int row_id) { return psiMinv_.device_data() + row_id * psiMinv_.cols(); }

  template<typename GT>
  static void mw_evalGrad(const RefVectorWithLeader<This_t>& engines,
                          const std::vector<const Value*>& dpsiM_row_list,
                          int rowchanged,
                          std::vector<GT>& grad_now)
  {
    auto& engine_leader = engines.getLeader();
    auto& mw_mem        = engine_leader.mw_mem_handle_.getResource();
    auto& buffer_H2D    = mw_mem.buffer_H2D;
    auto& grads_value_v = mw_mem.grads_value_v;

    const int norb                   = engine_leader.get_psiMinv().rows();
    const int nw                     = engines.size();
    constexpr size_t num_ptrs_packed = 2; // it must match packing and unpacking
    buffer_H2D.resize(sizeof(Value*) * num_ptrs_packed * nw);
    Matrix<const Value*> ptr_buffer(reinterpret_cast<const Value**>(buffer_H2D.data()), num_ptrs_packed, nw);
    for (int iw = 0; iw < nw; iw++)
    {
      ptr_buffer[0][iw] = engines[iw].get_psiMinv().device_data() + rowchanged * engine_leader.get_psiMinv().cols();
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
      const Value* __restrict__ invRow_ptr    = reinterpret_cast<const Value**>(buffer_H2D_ptr)[iw];
      const Value* __restrict__ dpsiM_row_ptr = reinterpret_cast<const Value**>(buffer_H2D_ptr)[nw + iw];
      Value grad_x(0), grad_y(0), grad_z(0);
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

  template<typename GT>
  static void mw_evalGradWithSpin(const RefVectorWithLeader<This_t>& engines,
                                  const std::vector<const Value*>& dpsiM_row_list,
                                  OffloadMatrix<Complex>& mw_dspin,
                                  int rowchanged,
                                  std::vector<GT>& grad_now,
                                  std::vector<Complex>& spingrad_now)
  {
    auto& engine_leader     = engines.getLeader();
    auto& mw_mem            = engine_leader.mw_mem_handle_.getResource();
    auto& buffer_H2D        = mw_mem.buffer_H2D;
    auto& grads_value_v     = mw_mem.grads_value_v;
    auto& spingrads_value_v = mw_mem.spingrads_value_v;

    //Need to pack these into a transfer buffer since psiMinv and dpsiM_row_list are not multiwalker data
    //i.e. each engine has its own psiMinv which is an OffloadMatrix instead of the leader having the data for all the walkers in the crowd.
    //Wouldn't have to do this if dpsiM and psiMinv were part of the mw_mem_handle_ with data across all walkers in the crowd and could just use use_device_ptr for the offload.
    //That is how mw_dspin is handled below
    const int norb                   = engine_leader.get_psiMinv().rows();
    const int nw                     = engines.size();
    constexpr size_t num_ptrs_packed = 2; // it must match packing and unpacking
    buffer_H2D.resize(sizeof(Value*) * num_ptrs_packed * nw);
    Matrix<const Value*> ptr_buffer(reinterpret_cast<const Value**>(buffer_H2D.data()), num_ptrs_packed, nw);
    for (int iw = 0; iw < nw; iw++)
    {
      ptr_buffer[0][iw] = engines[iw].get_psiMinv().device_data() + rowchanged * engine_leader.get_psiMinv().cols();
      ptr_buffer[1][iw] = dpsiM_row_list[iw];
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

  template<typename VVT>
  inline void updateRow(int rowchanged, const VVT& phiV, FullPrecValue c_ratio_in)
  {
    auto& Ainv = psiMinv_;
    // update the inverse matrix
    constexpr Value cone(1);
    constexpr Value czero(0);
    const int norb = Ainv.rows();
    const int lda  = Ainv.cols();
    temp.resize(norb);
    rcopy.resize(norb);
    // invoke the Fahy's variant of Sherman-Morrison update.
    int dummy_handle      = 0;
    int success           = 0;
    const Value* phiV_ptr = phiV.data();
    Value* Ainv_ptr       = Ainv.data();
    Value* temp_ptr       = temp.data();
    Value* rcopy_ptr      = rcopy.data();
    PRAGMA_OFFLOAD("omp target data map(always, to: phiV_ptr[:norb]) \
                    map(always, from: Ainv_ptr[:Ainv.size()]) \
                    use_device_ptr(phiV_ptr, Ainv_ptr, temp_ptr, rcopy_ptr)")
    {
      success = ompBLAS::gemv(dummy_handle, 'T', norb, norb, cone, Ainv_ptr, lda, phiV_ptr, 1, czero, temp_ptr, 1);
      if (success != 0)
        throw std::runtime_error("ompBLAS::gemv failed.");

      PRAGMA_OFFLOAD("omp target parallel for simd is_device_ptr(Ainv_ptr, temp_ptr, rcopy_ptr)")
      for (int i = 0; i < norb; i++)
      {
        rcopy_ptr[i] = Ainv_ptr[rowchanged * lda + i];
        if (i == 0)
          temp_ptr[rowchanged] -= cone;
      }

      success = ompBLAS::ger(dummy_handle, norb, norb, static_cast<Value>(-1.0 / c_ratio_in), rcopy_ptr, 1, temp_ptr, 1,
                             Ainv_ptr, lda);
      if (success != 0)
        throw std::runtime_error("ompBLAS::ger failed.");
    }
  }

  /** The potential for mayhem here without a unit test is great.
   */
  static void mw_updateRow(const RefVectorWithLeader<This_t>& engines,
                           int rowchanged,
                           const std::vector<Value*>& psiM_g_list,
                           const std::vector<Value*>& psiM_l_list,
                           const std::vector<bool>& isAccepted,
                           const OffloadMWVGLArray<Value>& phi_vgl_v,
                           const std::vector<Value>& ratios)
  {
    const size_t n_accepted = psiM_g_list.size();
    if (n_accepted == 0)
      return;

    auto& engine_leader         = engines.getLeader();
    auto& mw_mem                = engine_leader.mw_mem_handle_.getResource();
    auto& buffer_H2D            = mw_mem.buffer_H2D;
    auto& grads_value_v         = mw_mem.grads_value_v;
    auto& cone_vec              = mw_mem.cone_vec;
    auto& czero_vec             = mw_mem.czero_vec;
    auto& mw_temp               = mw_mem.mw_temp;
    auto& mw_rcopy              = mw_mem.mw_rcopy;
    const int norb              = engine_leader.get_psiMinv().rows();
    const int lda               = engine_leader.get_psiMinv().cols();
    const size_t nw             = isAccepted.size();
    const size_t phi_vgl_stride = nw * norb;

    engine_leader.resize_scratch_arrays(norb, n_accepted);

    // to handle Value** of Ainv, psi_v, temp, rcopy
    constexpr size_t num_ptrs_packed = 6; // it must match packing and unpacking
    buffer_H2D.resize((sizeof(Value*) * num_ptrs_packed + sizeof(Value)) * n_accepted);
    Matrix<Value*> ptr_buffer(reinterpret_cast<Value**>(buffer_H2D.data()), num_ptrs_packed, n_accepted);
    Value* c_ratio_inv = reinterpret_cast<Value*>(buffer_H2D.data() + sizeof(Value*) * num_ptrs_packed * n_accepted);
    for (int iw = 0, count = 0; iw < nw; iw++)
      if (isAccepted[iw])
      {
        ptr_buffer[0][count] = engines[iw].get_ref_psiMinv().device_data();
        ptr_buffer[1][count] = const_cast<Value*>(phi_vgl_v.device_data_at(0, iw, 0));
        ptr_buffer[2][count] = mw_temp.device_data() + norb * count;
        ptr_buffer[3][count] = mw_rcopy.device_data() + norb * count;
        ptr_buffer[4][count] = psiM_g_list[count];
        ptr_buffer[5][count] = psiM_l_list[count];

        c_ratio_inv[count] = Value(-1) / ratios[iw];
        count++;
      }

    // update the inverse matrix
    constexpr Value cone(1);
    constexpr Value czero(0);
    int dummy_handle     = 0;
    int success          = 0;
    auto* buffer_H2D_ptr = buffer_H2D.data();
    engine_leader.resize_fill_constant_arrays(n_accepted);
    Value* cone_ptr  = cone_vec.data();
    Value* czero_ptr = czero_vec.data();
    PRAGMA_OFFLOAD("omp target data \
                    map(always, to: buffer_H2D_ptr[:buffer_H2D.size()]) \
                    use_device_ptr(buffer_H2D_ptr, cone_ptr, czero_ptr)")
    {
      Value** Ainv_mw_ptr   = reinterpret_cast<Value**>(buffer_H2D_ptr);
      Value** phiVGL_mw_ptr = reinterpret_cast<Value**>(buffer_H2D_ptr + sizeof(Value*) * n_accepted);
      Value** temp_mw_ptr   = reinterpret_cast<Value**>(buffer_H2D_ptr + sizeof(Value*) * n_accepted * 2);
      Value** rcopy_mw_ptr  = reinterpret_cast<Value**>(buffer_H2D_ptr + sizeof(Value*) * n_accepted * 3);
      Value** dpsiM_mw_out  = reinterpret_cast<Value**>(buffer_H2D_ptr + sizeof(Value*) * n_accepted * 4);
      Value** d2psiM_mw_out = reinterpret_cast<Value**>(buffer_H2D_ptr + sizeof(Value*) * n_accepted * 5);
      Value* ratio_inv_mw   = reinterpret_cast<Value*>(buffer_H2D_ptr + sizeof(Value*) * n_accepted * 6);

      // invoke the Fahy's variant of Sherman-Morrison update.
      success = ompBLAS::gemv_batched(dummy_handle, 'T', norb, norb, cone_ptr, Ainv_mw_ptr, lda, phiVGL_mw_ptr, 1,
                                      czero_ptr, temp_mw_ptr, 1, n_accepted);
      if (success != 0)
        throw std::runtime_error("ompBLAS::gemv_batched failed.");

      PRAGMA_OFFLOAD("omp target teams distribute num_teams(n_accepted) is_device_ptr(Ainv_mw_ptr, temp_mw_ptr, \
                     rcopy_mw_ptr, dpsiM_mw_out, d2psiM_mw_out, phiVGL_mw_ptr)")
      for (int iw = 0; iw < n_accepted; iw++)
      {
        Value* __restrict__ Ainv_ptr   = Ainv_mw_ptr[iw];
        Value* __restrict__ temp_ptr   = temp_mw_ptr[iw];
        Value* __restrict__ rcopy_ptr  = rcopy_mw_ptr[iw];
        Value* __restrict__ dpsiM_out  = dpsiM_mw_out[iw];
        Value* __restrict__ d2psiM_out = d2psiM_mw_out[iw];
        Value* __restrict__ dpsiM_in   = phiVGL_mw_ptr[iw] + phi_vgl_stride;
        Value* __restrict__ d2psiM_in  = phiVGL_mw_ptr[iw] + phi_vgl_stride * 4;

        temp_ptr[rowchanged] -= cone;
        PRAGMA_OFFLOAD("omp parallel for simd")
        for (int i = 0; i < norb; i++)
        {
          rcopy_ptr[i] = Ainv_ptr[rowchanged * lda + i];
          // the following copying data on the device is not part of SM-1
          // it is intended to copy dpsiM and d2psiM from temporary to final without a separate kernel.
          dpsiM_out[i * 3]     = dpsiM_in[i];
          dpsiM_out[i * 3 + 1] = dpsiM_in[i + phi_vgl_stride];
          dpsiM_out[i * 3 + 2] = dpsiM_in[i + phi_vgl_stride * 2];
          d2psiM_out[i]        = d2psiM_in[i];
        }
      }

      success = ompBLAS::ger_batched(dummy_handle, norb, norb, ratio_inv_mw, rcopy_mw_ptr, 1, temp_mw_ptr, 1,
                                     Ainv_mw_ptr, lda, n_accepted);
      if (success != 0)
        throw std::runtime_error("ompBLAS::ger failed.");
    }
  }

  static void mw_accept_rejectRow(const RefVectorWithLeader<This_t>& engines,
                                  const int rowchanged,
                                  const std::vector<Value*>& psiM_g_list,
                                  const std::vector<Value*>& psiM_l_list,
                                  const std::vector<bool>& isAccepted,
                                  const OffloadMWVGLArray<Value>& phi_vgl_v,
                                  const std::vector<Value>& ratios)
  {
    mw_updateRow(engines, rowchanged, psiM_g_list, psiM_l_list, isAccepted, phi_vgl_v, ratios);
  }

  /** update the full Ainv and reset delay_count
   * @param Ainv inverse matrix
   */
  inline static void mw_updateInvMat(const RefVectorWithLeader<This_t>& engines) {}

  std::vector<const Value*> static mw_getInvRow(const RefVectorWithLeader<This_t>& engines,
                                                const int row_id,
                                                bool on_host)
  {
    const size_t nw    = engines.size();
    const size_t ncols = engines.getLeader().get_psiMinv().cols();
    std::vector<const Value*> row_ptr_list;
    row_ptr_list.reserve(nw);
    if (on_host)
      for (This_t& engine : engines)
      {
        auto* ptr = engine.get_psiMinv().data();
        PRAGMA_OFFLOAD("omp target update from(ptr[row_id * ncols : ncols])")
        row_ptr_list.push_back(ptr + row_id * ncols);
      }
    else
      for (This_t& engine : engines)
        row_ptr_list.push_back(engine.get_psiMinv().device_data() + row_id * ncols);
    return row_ptr_list;
  }

  /// transfer Ainv to the host
  static void mw_transferAinv_D2H(const RefVectorWithLeader<This_t>& engines)
  {
    for (This_t& engine : engines)
    {
      auto* ptr = engine.get_psiMinv().data();
      PRAGMA_OFFLOAD("omp target update from(ptr[:engine.get_psiMinv().size()])")
    }
  }

  /** transfer psiM_vgl to the host. psiM_vgl has 5 rows, V(1) G(3) L(1)
   * @param engine_leader for accessing shared resource
   * @param psiM_vgl_list list of psiM_vgl
   * @param row_begin first row to copy
   * @param row_size the number of rows to be copied
   */
  static void mw_transferVGL_D2H(This_t& engine_leader,
                                 const RefVector<OffloadVGLVector<Value>>& psiM_vgl_list,
                                 size_t row_begin,
                                 size_t row_size)
  {
    for (OffloadVGLVector<Value>& psiM_vgl : psiM_vgl_list)
    {
      auto* psiM_vgl_ptr  = psiM_vgl.data();
      const size_t stride = psiM_vgl.capacity();
      PRAGMA_OFFLOAD("omp target update from(psiM_vgl_ptr[row_begin * stride: row_size * stride]) nowait")
    }

    // The only tasks being waited on here are the psiM_vgl updates
    PRAGMA_OFFLOAD("omp taskwait")
  }

  /** transfer psiM_vgl to the device. psiM_vgl has 5 rows, V(1) G(3) L(1)
   * @param engine_leader for accessing shared resource
   * @param psiM_vgl_list list of psiM_vgl
   * @param row_begin first row to copy
   * @param row_size the number of rows to be copied
   */
  static void mw_transferVGL_H2D(This_t& engine_leader,
                                 const RefVector<OffloadVGLVector<Value>>& psiM_vgl_list,
                                 size_t row_begin,
                                 size_t row_size)
  {
    for (OffloadVGLVector<Value>& psiM_vgl : psiM_vgl_list)
    {
      auto* psiM_vgl_ptr  = psiM_vgl.data();
      const size_t stride = psiM_vgl.capacity();
      PRAGMA_OFFLOAD("omp target update to(psiM_vgl_ptr[row_begin * stride: row_size * stride]) nowait")
    }

    // The only tasks being waited on here are the psiM_vgl updates
    PRAGMA_OFFLOAD("omp taskwait")
  }

  auto& getLAhandles() { return dummy; }
};
} // namespace qmcplusplus

#endif // QMCPLUSPLUS_MATRIX_UPDATE_H
