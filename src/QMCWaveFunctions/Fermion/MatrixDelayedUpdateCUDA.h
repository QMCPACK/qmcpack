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

#ifndef QMCPLUSPLUS_MATRIX_DELAYED_UPDATE_CUDA_H
#define QMCPLUSPLUS_MATRIX_DELAYED_UPDATE_CUDA_H

#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Fermion/DiracMatrixComputeCUDA.hpp"
#include "Fermion/DiracMatrix.h"
#include "Platforms/OMPTarget/ompBLAS.hpp"
#include "CUDA/CUDAruntime.hpp"
#include "CUDA/cuBLAS.hpp"
#include "CUDA/cuBLAS_missing_functions.hpp"
#include "QMCWaveFunctions/detail/CUDA/matrix_update_helper.hpp"
#include "DualAllocatorAliases.hpp"
#include "CUDA/CUDALinearAlgebraHandles.h"
#include "ResourceCollection.h"
#include "WaveFunctionTypes.hpp"

namespace qmcplusplus
{

namespace testing
{
class DiracDeterminantBatchedTest;
}

/** Implements dirac matrix delayed update using OpenMP offload and CUDA.
 * It is used as DET_ENGINE in DiracDeterminantBatched.
 * This is a 1 per walker class unlike DiracDeterminantBatched and
 * DiracMatrixComputeCUDA which are 1 per many walkers
 *
 * @tparam T base precision for most computation
 * @tparam T_FP high precision for matrix inversion, T_FP >= T
 */
template<typename T, typename T_FP>
class MatrixDelayedUpdateCUDA
{
public:
  using WFT           = WaveFunctionTypes<T, T_FP>;
  using Value         = typename WFT::Value;
  using FullPrecValue = typename WFT::FullPrecValue;
  using LogValue      = typename WFT::LogValue;
  using This_t        = MatrixDelayedUpdateCUDA<T, T_FP>;

  // allows enforcement of OMPallocator matched between DiracDeterminantBatched and MatrixUpdateOMPTarget
  template<typename DT>
  using PinnedDualAllocator = PinnedDualAllocator<DT>;
  template<typename DT>
  using UnpinnedDualAllocator = UnpinnedDualAllocator<DT>;

  // Want to emphasize these because at least for cuda they can't be transferred async, which is bad.
  template<typename DT>
  using UnpinnedDualVector = Vector<DT, UnpinnedDualAllocator<DT>>;
  template<typename DT>
  using DualVector = Vector<DT, PinnedDualAllocator<DT>>;
  template<typename DT>
  using DualMatrix = Matrix<DT, PinnedDualAllocator<DT>>;

  struct MatrixDelayedUpdateCUDAMultiWalkerMem : public Resource
  {
    // constant array value T(1)
    UnpinnedDualVector<T> cone_vec;
    // constant array value T(-1)
    UnpinnedDualVector<T> cminusone_vec;
    // constant array value T(0)
    UnpinnedDualVector<T> czero_vec;
    // multi walker of grads for transfer needs.
    DualMatrix<T> grads_value_v;
    // mw_updateRow pointer buffer
    Vector<char, PinnedDualAllocator<char>> updateRow_buffer_H2D;
    // mw_prepareInvRow pointer buffer
    Vector<char, PinnedDualAllocator<char>> prepare_inv_row_buffer_H2D;
    // mw_accept_rejectRow pointer buffer
    Vector<char, PinnedDualAllocator<char>> accept_rejectRow_buffer_H2D;
    // mw_updateInv pointer buffer
    Vector<char, PinnedDualAllocator<char>> updateInv_buffer_H2D;
    // mw_evalGrad pointer buffer
    Vector<char, PinnedDualAllocator<char>> evalGrad_buffer_H2D;
    /// scratch space for rank-1 update
    UnpinnedDualVector<T> mw_temp;
    // scratch space for keeping one row of Ainv
    UnpinnedDualVector<T> mw_rcopy;

    MatrixDelayedUpdateCUDAMultiWalkerMem() : Resource("MatrixDelayedUpdateCUDAMultiWalkerMem") {}

    MatrixDelayedUpdateCUDAMultiWalkerMem(const MatrixDelayedUpdateCUDAMultiWalkerMem&)
        : MatrixDelayedUpdateCUDAMultiWalkerMem()
    {}

    Resource* makeClone() const override { return new MatrixDelayedUpdateCUDAMultiWalkerMem(*this); }
  };

private:
  /// legacy single walker matrix inversion engine
  DiracMatrix<T_FP> detEng;

  /// inverse transpose of psiM(j,i) \f$= \psi_j({\bf r}_i)\f$
  DualMatrix<Value> psiMinv_;
  /// scratch space for rank-1 update
  UnpinnedDualVector<Value> temp;
  /// row of up-to-date Ainv
  UnpinnedDualVector<Value> invRow;
  /** row id correspond to the up-to-date invRow. [0 norb), invRow is ready; -1, invRow is not valid.
   *  This id is set after calling getInvRow indicating invRow has been prepared for the invRow_id row
   *  ratioGrad checks if invRow_id is consistent. If not, invRow needs to be recomputed.
   *  acceptMove and completeUpdates mark invRow invalid by setting invRow_id to -1
   */
  int inv_row_id_;
  // scratch space for keeping one row of Ainv
  UnpinnedDualVector<Value> rcopy;

  template<typename DT>
  using DeviceMatrix = Matrix<DT, CUDAAllocator<DT>>;
  template<typename DT>
  using DeviceVector = Vector<DT, CUDAAllocator<DT>>;
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
  Vector<int, CUDAAllocator<int>> delay_list_gpu;
  /// current number of delays, increase one for each acceptance, reset to 0 after updating Ainv
  int delay_count;

  // psi(r')/psi(r) during a PbyP move
  FullPrecValue cur_ratio_;

  /** @ingroup Resources
   *  @{ */
  // CUDA stream, cublas handle object
  std::unique_ptr<CUDALinearAlgebraHandles> cuda_handles_;
  /// matrix inversion engine this a crowd scope resource and only the leader engine gets it
  UPtr<DiracMatrixComputeCUDA<T_FP>> det_inverter_;
  /**}@ */

  std::unique_ptr<MatrixDelayedUpdateCUDAMultiWalkerMem> mw_mem_;

  inline void waitStream()
  {
    cudaErrorCheck(cudaStreamSynchronize(cuda_handles_->hstream), "cudaStreamSynchronize failed!");
  }

  /** ensure no previous delay left.
   *  This looks like it should be an assert
   */
  inline void guard_no_delay() const
  {
    if (delay_count != 0)
      throw std::runtime_error("BUG: unexpected call sequence delay_count is not 0");
  }

  // check if the number of maximal delay is 1 (SM-1)
  // \todo rename this something containing delay.
  inline bool isSM1() const { return Binv_gpu.rows() == 1; }

  /** a bad smell */
  void resize_fill_constant_arrays(size_t nw)
  {
    if (mw_mem_->cone_vec.size() < nw)
    {
      // cone
      mw_mem_->cone_vec.resize(nw);
      std::fill_n(mw_mem_->cone_vec.data(), nw, T(1));
      T* cone_ptr = mw_mem_->cone_vec.data();
      PRAGMA_OFFLOAD("omp target update to(cone_ptr[:nw])")
      // cminusone
      mw_mem_->cminusone_vec.resize(nw);
      std::fill_n(mw_mem_->cminusone_vec.data(), nw, T(-1));
      T* cminusone_ptr = mw_mem_->cminusone_vec.data();
      PRAGMA_OFFLOAD("omp target update to(cminusone_ptr[:nw])")
      // czero
      mw_mem_->czero_vec.resize(nw);
      std::fill_n(mw_mem_->czero_vec.data(), nw, T(0));
      T* czero_ptr = mw_mem_->czero_vec.data();
      PRAGMA_OFFLOAD("omp target update to(czero_ptr[:nw])")
    }
  }

  /** compute the row of up-to-date Ainv
   *  What the hell is this?
   * @param Ainv inverse matrix
   * @param rowchanged the row id corresponding to the proposed electron
   */
  static void mw_prepareInvRow(const RefVectorWithLeader<This_t>& engines, const int rowchanged)
  {
    auto& engine_leader              = engines.getLeader();
    auto& hstream                    = engine_leader.cuda_handles_->hstream;
    auto& h_cublas                   = engine_leader.cuda_handles_->h_cublas;
    auto& cminusone_vec              = engine_leader.mw_mem_->cminusone_vec;
    auto& cone_vec                   = engine_leader.mw_mem_->cone_vec;
    auto& czero_vec                  = engine_leader.mw_mem_->czero_vec;
    auto& prepare_inv_row_buffer_H2D = engine_leader.mw_mem_->prepare_inv_row_buffer_H2D;
    const int norb                   = engine_leader.psiMinv_.rows();
    const int nw                     = engines.size();
    int& delay_count                 = engine_leader.delay_count;
    prepare_inv_row_buffer_H2D.resize(sizeof(T*) * 7 * nw);
    engine_leader.resize_fill_constant_arrays(nw);

    const int lda_Binv = engine_leader.Binv_gpu.cols();
    Matrix<T*> ptr_buffer(reinterpret_cast<T**>(prepare_inv_row_buffer_H2D.data()), 7, nw);
    for (int iw = 0; iw < nw; iw++)
    {
      This_t& engine    = engines[iw];
      ptr_buffer[0][iw] = engine.psiMinv_.device_data() + rowchanged * engine.psiMinv_.cols();
      ptr_buffer[1][iw] = engine.invRow.device_data();
      ptr_buffer[2][iw] = engine.U_gpu.data();
      ptr_buffer[3][iw] = engine.p_gpu.data();
      ptr_buffer[4][iw] = engine.Binv_gpu.data();
      ptr_buffer[5][iw] = engine.Binv_gpu.data() + delay_count * lda_Binv;
      ptr_buffer[6][iw] = engine.V_gpu.data();
    }

    cudaErrorCheck(cudaMemcpyAsync(prepare_inv_row_buffer_H2D.device_data(), prepare_inv_row_buffer_H2D.data(),
                                   prepare_inv_row_buffer_H2D.size(), cudaMemcpyHostToDevice, hstream),
                   "cudaMemcpyAsync prepare_inv_row_buffer_H2D failed!");

    T** oldRow_mw_ptr  = reinterpret_cast<T**>(prepare_inv_row_buffer_H2D.device_data());
    T** invRow_mw_ptr  = reinterpret_cast<T**>(prepare_inv_row_buffer_H2D.device_data() + sizeof(T*) * nw);
    T** U_mw_ptr       = reinterpret_cast<T**>(prepare_inv_row_buffer_H2D.device_data() + sizeof(T*) * nw * 2);
    T** p_mw_ptr       = reinterpret_cast<T**>(prepare_inv_row_buffer_H2D.device_data() + sizeof(T*) * nw * 3);
    T** Binv_mw_ptr    = reinterpret_cast<T**>(prepare_inv_row_buffer_H2D.device_data() + sizeof(T*) * nw * 4);
    T** BinvRow_mw_ptr = reinterpret_cast<T**>(prepare_inv_row_buffer_H2D.device_data() + sizeof(T*) * nw * 5);
    T** V_mw_ptr       = reinterpret_cast<T**>(prepare_inv_row_buffer_H2D.device_data() + sizeof(T*) * nw * 6);

    // save Ainv[rowchanged] to invRow
    //std::copy_n(Ainv[rowchanged], norb, invRow.data());
    cudaErrorCheck(cuBLAS_MFs::copy_batched(hstream, norb, oldRow_mw_ptr, 1, invRow_mw_ptr, 1, nw),
                   "cuBLAS_MFs::copy_batched failed!");
    // multiply V (NxK) Binv(KxK) U(KxN) invRow right to the left
    //BLAS::gemv('T', norb, delay_count, cone, U_gpu.data(), norb, invRow.data(), 1, czero, p_gpu.data(), 1);
    //BLAS::gemv('N', delay_count, delay_count, -cone, Binv.data(), lda_Binv, p.data(), 1, czero, Binv[delay_count], 1);
    //BLAS::gemv('N', norb, delay_count, cone, V.data(), norb, Binv[delay_count], 1, cone, invRow.data(), 1);
    cudaErrorCheck(cuBLAS_MFs::gemv_batched(hstream, 'T', norb, delay_count, cone_vec.device_data(), U_mw_ptr, norb,
                                            invRow_mw_ptr, 1, czero_vec.device_data(), p_mw_ptr, 1, nw),
                   "cuBLAS_MFs::gemv_batched failed!");
    cudaErrorCheck(cuBLAS_MFs::gemv_batched(hstream, 'N', delay_count, delay_count, cminusone_vec.device_data(),
                                            Binv_mw_ptr, lda_Binv, p_mw_ptr, 1, czero_vec.device_data(), BinvRow_mw_ptr,
                                            1, nw),
                   "cuBLAS_MFs::gemv_batched failed!");
    cudaErrorCheck(cuBLAS_MFs::gemv_batched(hstream, 'N', norb, delay_count, cone_vec.device_data(), V_mw_ptr, norb,
                                            BinvRow_mw_ptr, 1, cone_vec.device_data(), invRow_mw_ptr, 1, nw),
                   "cuBLAS_MFs::gemv_batched failed!");
    // mark row prepared
    engine_leader.inv_row_id() = rowchanged;
  }

  /** Do complete row updates
   *  many of these const arguments provide pointers or references
   *  somwhere in here is an update that doesn't get where it belongs resulting in a 0
   *  gradient later.
   *  Sad example of OpenMP target code that is far from clear and a poor substitute for a
   *  clear CPU reference implementation.
   *
   *  \param[in] engines
   *  \param[in] rowchanged
   *  \param[in] psiM_g_list        device ptrs
   *  \param[in] psiM_l_list        device ptrs
   *  \param[in] isAccepted         bool but wait some lists are also filtered
   *  \param[in] phi_vgl_v_dev_ptr  device ptr
   *  \param[in] phi_vgl_stride     size of each "vector" in phi_vgl_v
   *  \param[inout] ratios
   */
  static void mw_updateRow(const RefVectorWithLeader<This_t>& engines,
                           const int rowchanged,
                           const std::vector<T*>& psiM_g_list,
                           const std::vector<T*>& psiM_l_list,
                           const std::vector<bool>& isAccepted,
                           const T* phi_vgl_v_dev_ptr,
                           const size_t phi_vgl_stride,
                           const std::vector<T>& ratios)
  {
    auto& engine_leader = engines.getLeader();
    engine_leader.guard_no_delay();

    // I think this should be the number of true's in isAccepted
    // But actually everything except phi_vgl_v needs to be filtered to this size.
    // I think it would be much simpler to not filter and just have everything be size of isAccepted.
    const size_t n_accepted = psiM_g_list.size();
#ifndef NDEBUG
    size_t n_true = std::count_if(isAccepted.begin(), isAccepted.end(), [](bool accepted) { return accepted; });
    assert(n_accepted == n_true);
#endif
    if (n_accepted == 0)
      return;
    engine_leader.resize_fill_constant_arrays(n_accepted);

    auto& hstream              = engine_leader.cuda_handles_->hstream;
    auto& updateRow_buffer_H2D = engine_leader.mw_mem_->updateRow_buffer_H2D;
    auto& mw_temp              = engine_leader.mw_mem_->mw_temp;
    auto& mw_rcopy             = engine_leader.mw_mem_->mw_rcopy;
    auto& cone_vec             = engine_leader.mw_mem_->cone_vec;
    auto& czero_vec            = engine_leader.mw_mem_->czero_vec;
    const int norb             = engine_leader.psiMinv_.rows();
    const int lda              = engine_leader.psiMinv_.cols();
    mw_temp.resize(norb * n_accepted);
    mw_rcopy.resize(norb * n_accepted);
    updateRow_buffer_H2D.resize((sizeof(T*) * 8 * n_accepted + sizeof(T)) * n_accepted);

    // to handle T** of Ainv, psi_v, temp, rcopy
    Matrix<T*> ptr_buffer(reinterpret_cast<T**>(updateRow_buffer_H2D.data()), 8, n_accepted);
    T* c_ratio_inv = reinterpret_cast<T*>(updateRow_buffer_H2D.data() + sizeof(T*) * 8 * n_accepted);
    for (int iw = 0, count = 0; iw < isAccepted.size(); iw++)
      if (isAccepted[iw])
      {
        ptr_buffer[0][count] = engines[iw].psiMinv_.device_data();
        ptr_buffer[1][count] = const_cast<T*>(phi_vgl_v_dev_ptr + norb * iw);
        ptr_buffer[2][count] = mw_temp.device_data() + norb * count;
        ptr_buffer[3][count] = mw_rcopy.device_data() + norb * count;
        ptr_buffer[4][count] = psiM_g_list[count];
        ptr_buffer[5][count] = psiM_l_list[count];
        ptr_buffer[6][count] = const_cast<T*>(phi_vgl_v_dev_ptr + phi_vgl_stride + norb * 3 * iw);
        ptr_buffer[7][count] = const_cast<T*>(phi_vgl_v_dev_ptr + phi_vgl_stride * 4 + norb * iw);
        c_ratio_inv[count]   = T(-1) / ratios[iw];
        count++;
      }

    // update the inverse matrix
    cudaErrorCheck(cudaMemcpyAsync(updateRow_buffer_H2D.device_data(), updateRow_buffer_H2D.data(),
                                   updateRow_buffer_H2D.size(), cudaMemcpyHostToDevice, hstream),
                   "cudaMemcpyAsync updateRow_buffer_H2D failed!");

    {
      T** Ainv_mw_ptr   = reinterpret_cast<T**>(updateRow_buffer_H2D.device_data());
      T** phiV_mw_ptr   = reinterpret_cast<T**>(updateRow_buffer_H2D.device_data() + sizeof(T*) * n_accepted);
      T** temp_mw_ptr   = reinterpret_cast<T**>(updateRow_buffer_H2D.device_data() + sizeof(T*) * n_accepted * 2);
      T** rcopy_mw_ptr  = reinterpret_cast<T**>(updateRow_buffer_H2D.device_data() + sizeof(T*) * n_accepted * 3);
      T** dpsiM_mw_out  = reinterpret_cast<T**>(updateRow_buffer_H2D.device_data() + sizeof(T*) * n_accepted * 4);
      T** d2psiM_mw_out = reinterpret_cast<T**>(updateRow_buffer_H2D.device_data() + sizeof(T*) * n_accepted * 5);
      T** dpsiM_mw_in   = reinterpret_cast<T**>(updateRow_buffer_H2D.device_data() + sizeof(T*) * n_accepted * 6);
      T** d2psiM_mw_in  = reinterpret_cast<T**>(updateRow_buffer_H2D.device_data() + sizeof(T*) * n_accepted * 7);
      T* ratio_inv_mw   = reinterpret_cast<T*>(updateRow_buffer_H2D.device_data() + sizeof(T*) * n_accepted * 8);

      // invoke the Fahy's variant of Sherman-Morrison update.
      cudaErrorCheck(cuBLAS_MFs::gemv_batched(hstream, 'T', norb, norb, cone_vec.device_data(), Ainv_mw_ptr, lda,
                                              phiV_mw_ptr, 1, czero_vec.device_data(), temp_mw_ptr, 1, n_accepted),
                     "cuBLAS_MFs::gemv_batched failed!");

      cudaErrorCheck(CUDA::copyAinvRow_saveGL_cuda(hstream, rowchanged, norb, Ainv_mw_ptr, lda, temp_mw_ptr,
                                                   rcopy_mw_ptr, dpsiM_mw_in, d2psiM_mw_in, dpsiM_mw_out, d2psiM_mw_out,
                                                   n_accepted),
                     "CUDA::copyAinvRow_saveGL_cuda failed!");


      cudaErrorCheck(cuBLAS_MFs::ger_batched(hstream, norb, norb, ratio_inv_mw, rcopy_mw_ptr, 1, temp_mw_ptr, 1,
                                             Ainv_mw_ptr, lda, n_accepted),
                     "cuBLAS_MFs::ger_batched failed!");
    }
  }

public:
  /// default constructor
  MatrixDelayedUpdateCUDA() : inv_row_id_(-1), delay_count(0) {}

  MatrixDelayedUpdateCUDA(const MatrixDelayedUpdateCUDA&) = delete;

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
    psiMinv_.resize(norb, getAlignedSize<T>(norb));
  }

  void createResource(ResourceCollection& collection) const
  {
    //the semantics of the ResourceCollection are such that we don't want to add a Resource that we need
    //later in the chain of resource creation.
    auto clah_ptr = std::make_unique<CUDALinearAlgebraHandles>();
    auto dmcc_ptr = std::make_unique<DiracMatrixComputeCUDA<T_FP>>(clah_ptr->hstream);
    collection.addResource(std::move(clah_ptr));
    collection.addResource(std::move(dmcc_ptr));
    collection.addResource(std::make_unique<MatrixDelayedUpdateCUDAMultiWalkerMem>());
  }

  void acquireResource(ResourceCollection& collection)
  {
    auto res_ptr = dynamic_cast<CUDALinearAlgebraHandles*>(collection.lendResource().release());
    if (!res_ptr)
      throw std::runtime_error("MatrixDelayedUpdateCUDA::acquireResource dynamic_cast CUDALinearAlgebraHandles failed");
    cuda_handles_.reset(res_ptr);
    auto det_eng_ptr = dynamic_cast<DiracMatrixComputeCUDA<T_FP>*>(collection.lendResource().release());
    if (!det_eng_ptr)
      throw std::runtime_error(
          "MatrixDelayedUpdateCUDA::acquireResource dynamic_cast to DiracMatrixComputeCUDA<T_FP>* failed");
    det_inverter_.reset(det_eng_ptr);
    auto res2_ptr = dynamic_cast<MatrixDelayedUpdateCUDAMultiWalkerMem*>(collection.lendResource().release());
    if (!res2_ptr)
      throw std::runtime_error(
          "MatrixDelayedUpdateCUDA::acquireResource dynamic_cast MatrixDelayedUpdateCUDAMultiWalkerMem failed");
    mw_mem_.reset(res2_ptr);
  }

  void releaseResource(ResourceCollection& collection)
  {
    collection.takebackResource(std::move(cuda_handles_));
    collection.takebackResource(std::move(det_inverter_));
    collection.takebackResource(std::move(mw_mem_));
  }

  inline const DualMatrix<FullPrecValue>& get_psiMinv() const { return psiMinv_; }

  inline DualMatrix<FullPrecValue>& psiMinv() { return psiMinv_; }
  const int get_inv_row_id() const { return inv_row_id_; }
  int& inv_row_id() { return inv_row_id_; }

  inline T* getRow_psiMinv_offload(int row_id) { return psiMinv_.device_data() + row_id * psiMinv_.cols(); }

  QMCTraits::FullPrecRealType& cur_ratio() { return cur_ratio_; }
  const QMCTraits::FullPrecRealType get_cur_ratio() const { return cur_ratio_; }

  /** make this class unit tests friendly without the need of setup resources.
   *  belongs in a friend class in test
   */
  inline void checkResourcesForTest()
  {
    if (!cuda_handles_)
    {
      throw std::logic_error(
          "Null cuda_handles_, Even for testing proper resource creation and acquisition must be made.");
    }

    if (!det_inverter_)
    {
      throw std::logic_error(
          "Null det_inverter_, Even for testing proper resource creation and acquisition must be made.");
    }
  }


  inline void invert_transpose(const Matrix<Value>& log_det, Matrix<Value>& a_inv, LogValue& log_value)
  {
    guard_no_delay();
    detEng.invert_transpose(log_det, a_inv, log_value);
  }


  /** compute the inverse of the transpose of matrix logdetT, result is in psiMinv
   *
   *  This does not get called constantly so to get real benchmark data that redirection to mw
   *  is a big deal before optimizing.
   *  psiMinv is copied to the Host as a side effect
   * @param logdetT orbital value matrix
   * @param log_value log(det(logdetT))
   */
  inline void invert_transpose(DualMatrix<Value>& log_det, DualMatrix<Value>& a_inv, DualVector<LogValue>& log_values)
  {
    guard_no_delay();
    det_inverter_->invert_transpose(*cuda_handles_, log_det, a_inv, log_values);
    T* a_inv_ptr = a_inv.data();
    // This is likely better optional
    PRAGMA_OFFLOAD("omp target update to(a_inv_ptr[:a_inv.size()])")
    // For the above transfer
    PRAGMA_OFFLOAD("omp taskwait")
  }

  /** Compute the inversions of the transpose of matrices logdetT_list and calculate
   *  the log determinants of the logdetTs.
   *  if logdetT is of element type T_FP it will be returned fille with the LU matrix
   *  compute_mask is in the API to reserve the right to reduce transfers to/from device.
   */
  static void mw_invertTranspose(RefVectorWithLeader<MatrixDelayedUpdateCUDA<T, T_FP>>& engines,
                                 RefVector<DualMatrix<Value>>& logdetT_list,
                                 DualVector<LogValue>& log_values,
                                 const std::vector<bool>& compute_mask)
  {
    auto& engine_leader = engines.getLeader();

    engine_leader.guard_no_delay();

    RefVector<DualMatrix<Value>> a_inv_refs;
    a_inv_refs.reserve(engines.size());
    for (int iw = 0; iw < engines.size(); iw++)
      a_inv_refs.emplace_back(engines[iw].psiMinv_);
    engine_leader.get_det_inverter().mw_invertTranspose(*(engine_leader.cuda_handles_), logdetT_list, a_inv_refs,
                                                        log_values, compute_mask);
  }

  static void mw_invertTranspose(RefVectorWithLeader<MatrixDelayedUpdateCUDA<T, T_FP>>& engines,
                                 RefVector<DualMatrix<Value>>& logdetT_list,
                                 RefVector<DualMatrix<Value>>& a_inv_list,
                                 DualVector<LogValue>& log_values,
                                 const std::vector<bool>& compute_mask)
  {
    auto& engine_leader = engines.getLeader();
    engine_leader.guard_no_delay();
    engine_leader.get_det_inverter().mw_invertTranspose(*(engine_leader.cuda_handles_), logdetT_list, a_inv_list,
                                                        log_values, compute_mask);
  }

  // prepare invRow and compute the old gradients.
  template<typename GT>
  static void mw_evalGrad(const RefVectorWithLeader<This_t>& engines,
                          const std::vector<const T*>& dpsiM_row_list,
                          const int rowchanged,
                          std::vector<GT>& grad_now)
  {
    auto& engine_leader = engines.getLeader();
    if (!engine_leader.isSM1())
      mw_prepareInvRow(engines, rowchanged);

    auto& hstream             = engine_leader.cuda_handles_->hstream;
    auto& evalGrad_buffer_H2D = engine_leader.mw_mem_->evalGrad_buffer_H2D;
    auto& grads_value_v       = engine_leader.mw_mem_->grads_value_v;

    const int nw = engines.size();
    evalGrad_buffer_H2D.resize(sizeof(T*) * 2 * nw);
    Matrix<const T*> ptr_buffer(reinterpret_cast<const T**>(evalGrad_buffer_H2D.data()), 2, nw);
    for (int iw = 0; iw < nw; iw++)
    {
      if (engine_leader.isSM1())
        ptr_buffer[0][iw] = engines[iw].psiMinv_.device_data() + rowchanged * engine_leader.psiMinv_.cols();
      else
        ptr_buffer[0][iw] = engines[iw].invRow.device_data();
      ptr_buffer[1][iw] = dpsiM_row_list[iw];
    }

    cudaErrorCheck(cudaMemcpyAsync(evalGrad_buffer_H2D.device_data(), evalGrad_buffer_H2D.data(),
                                   evalGrad_buffer_H2D.size(), cudaMemcpyHostToDevice, hstream),
                   "cudaMemcpyAsync evalGrad_buffer_H2D failed!");

    if (grads_value_v.rows() != nw || grads_value_v.cols() != GT::Size)
      grads_value_v.resize(nw, GT::Size);

    const T** invRow_ptr    = reinterpret_cast<const T**>(evalGrad_buffer_H2D.device_data());
    const T** dpsiM_row_ptr = reinterpret_cast<const T**>(evalGrad_buffer_H2D.device_data()) + nw;

    const int norb = engine_leader.psiMinv_.rows();
    cudaErrorCheck(CUDA::calcGradients_cuda(hstream, norb, invRow_ptr, dpsiM_row_ptr, grads_value_v.device_data(), nw),
                   "CUDA::calcGradients_cuda failed!");
    cudaErrorCheck(cudaMemcpyAsync(grads_value_v.data(), grads_value_v.device_data(), grads_value_v.size() * sizeof(T),
                                   cudaMemcpyDeviceToHost, hstream),
                   "cudaMemcpyAsync grads_value_v failed!");
    engine_leader.waitStream();

    for (int iw = 0; iw < nw; iw++)
      grad_now[iw] = {grads_value_v[iw][0], grads_value_v[iw][1], grads_value_v[iw][2]};
  }

  /** Update the "local" psiMinv_ on the device.
   *  Side Effect Transfers:
   *  * phiV is left on host side in the single methods so it must be transferred to device
   *  * psiMinv_ is transferred back to host since other single calls from QMCHamitonian and others
   *  * expect it to be.
   *
   *  Forced to use OpenMP target since resources are banned for single walker functions APIs
   *  and the acquireRelease pattern for a single DDB was removed by #3324
   */
  template<typename VVT, typename RATIOT>
  void updateRow(int rowchanged, VVT& phiV, RATIOT c_ratio_in)
  {
    auto& Ainv = psiMinv_;
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
    // This must be Ainv must be tofrom due to NonlocalEcpComponent and possibly
    // other modules assumptions about the state of psiMinv. :(
    PRAGMA_OFFLOAD("omp target data map(always, to: phiV_ptr[:norb]) \
                    map(always, tofrom: Ainv_ptr[:Ainv.size()]) \
                    use_device_ptr(phiV_ptr, Ainv_ptr, temp_ptr, rcopy_ptr)")
    {
      int success = ompBLAS::gemv(dummy_handle, 'T', norb, norb, cone, Ainv_ptr, lda, phiV_ptr, 1, czero, temp_ptr, 1);
      if (success != 0)
        throw std::runtime_error("ompBLAS::gemv failed.");
      PRAGMA_OFFLOAD("omp target is_device_ptr(Ainv_ptr, temp_ptr, rcopy_ptr)")
      {
        temp_ptr[rowchanged] -= cone;
        PRAGMA_OFFLOAD("omp parallel for simd")
        for (int i = 0; i < norb; i++)
          rcopy_ptr[i] = Ainv_ptr[rowchanged * lda + i];
      }
      success = ompBLAS::ger(dummy_handle, norb, norb, static_cast<T>(RATIOT(-1) / c_ratio_in), rcopy_ptr, 1, temp_ptr,
                             1, Ainv_ptr, lda);
      if (success != 0)
        throw std::runtime_error("ompBLAS::ger failed.");
    }

    // this is a working CUDA only implementation. but only if DDB is allowed to access resources
    // implicitly available to openmp Target code through the openMP runtime.

    // guard_no_delay();

    // auto& Ainv = psiMinv_;
    // // update the inverse matrix
    // constexpr T cone(1), czero(0);
    // const int norb = Ainv.rows();
    // // This is Binv.cols() in DelayedUpdate
    // const int lda  = Ainv.cols();
    // temp.resize(norb);
    // rcopy.resize(norb);

    // // invoke the Fahy's variant of Sherman-Morrison update.
    // cudaErrorCheck(cudaMemcpyAsync(phiV.device_data(), phiV.data(), sizeof(T) * phiV.size(), cudaMemcpyHostToDevice, cuda_handles_->hstream), "cuda copy failed in MatrixDelayedUpdateCUDA::updateRow");

    // cublasErrorCheck(cuBLAS::gemv(cuda_handles_->h_cublas, CUBLAS_OP_T, norb, norb, &cone, psiMinv_.device_data(), lda, phiV.device_data(), 1, &czero, temp.device_data(), 1), "cuBLAS::gemv failed in MatrixDelayedUpdateCUDA::updateRow");
    // cublasErrorCheck(cuBLAS::copy(cuda_handles_->h_cublas, norb, Ainv.device_data() + rowchanged * lda, 1, rcopy.device_data(), 1), "cuBLAS::copy failed in MatrixDelayedUpdateCUDA::updateRow");
    // T alpha = static_cast<T>(RATIOT(-1) / c_ratio_in);
    // cublasErrorCheck(cuBLAS::ger(cuda_handles_->h_cublas, norb, norb, &alpha, rcopy.device_data(), 1, temp.device_data(), 1, psiMinv_.device_data(), lda), "cuBLAS::ger failed in MatrixDelayedUpdateCUDA::updateRow");
    // cudaErrorCheck(cudaMemcpyAsync(psiMinv_.data(), psiMinv_.device_data(), sizeof(T) * norb * norb, cudaMemcpyDeviceToHost, cuda_handles_->hstream), "cuda copy failed in MatrixDelayedUpdateCUDA::updateRow");
    // cudaErrorCheck(cudaStreamSynchronize(cuda_handles_->hstream), "cudaStreamSynchronize failed!");
  }

  /** Accept or Reject row updates
   *  many of these const arguments provide pointers or references
   *  to objects that do get modified.
   *  \param[in] engines
   *  \param[in] rowchanged
   *  \param[in] psiM_g_list
   *  \param[in] psiM_l_list
   *  \param[in] isAccepted
   *  \param[in] phi_vgl_v_dev_ptr
   *  \param[in] phi_vgl_stride     size of each "vector" in phi_vgl_v
   *  \param[inout] ratios
   */
  static void mw_accept_rejectRow(const RefVectorWithLeader<This_t>& engines,
                                  const int rowchanged,
                                  const std::vector<T*>& psiM_g_list,
                                  const std::vector<T*>& psiM_l_list,
                                  const std::vector<bool>& isAccepted,
                                  const T* phi_vgl_v_dev_ptr,
                                  const size_t phi_vgl_stride,
                                  const std::vector<T>& ratios)
  {
    auto& engine_leader = engines.getLeader();
    // invRow consumed, mark invRow_id unset
    engine_leader.inv_row_id() = -1;

    if (engine_leader.isSM1())
    {
      mw_updateRow(engines, rowchanged, psiM_g_list, psiM_l_list, isAccepted, phi_vgl_v_dev_ptr, phi_vgl_stride,
                   ratios);
      return;
    }

    auto& hstream                     = engine_leader.cuda_handles_->hstream;
    auto& cminusone_vec               = engine_leader.mw_mem_->cminusone_vec;
    auto& cone_vec                    = engine_leader.mw_mem_->cone_vec;
    auto& czero_vec                   = engine_leader.mw_mem_->czero_vec;
    auto& accept_rejectRow_buffer_H2D = engine_leader.mw_mem_->accept_rejectRow_buffer_H2D;
    int& delay_count                  = engine_leader.delay_count;
    const int lda_Binv                = engine_leader.Binv_gpu.cols();
    const int norb                    = engine_leader.psiMinv_.rows();
    const int lda                     = engine_leader.psiMinv_.cols();
    const int nw                      = engines.size();
    const int n_accepted              = psiM_g_list.size();
    accept_rejectRow_buffer_H2D.resize((sizeof(T*) * 14 + sizeof(T)) * nw);
    engine_leader.resize_fill_constant_arrays(nw);

    Matrix<T*> ptr_buffer(reinterpret_cast<T**>(accept_rejectRow_buffer_H2D.data()), 14, nw);
    T* c_ratio_inv = reinterpret_cast<T*>(accept_rejectRow_buffer_H2D.data() + sizeof(T*) * 14 * nw);
    for (int iw = 0, count_accepted = 0, count_rejected = 0; iw < nw; iw++)
    {
      This_t& engine = engines[iw];
      if (isAccepted[iw])
      {
        ptr_buffer[0][count_accepted]  = engine.psiMinv_.device_data() + lda * rowchanged;
        ptr_buffer[1][count_accepted]  = engine.V_gpu.data();
        ptr_buffer[2][count_accepted]  = engine.U_gpu.data() + norb * delay_count;
        ptr_buffer[3][count_accepted]  = engine.p_gpu.data();
        ptr_buffer[4][count_accepted]  = engine.Binv_gpu.data();
        ptr_buffer[5][count_accepted]  = engine.Binv_gpu.data() + delay_count * lda_Binv;
        ptr_buffer[6][count_accepted]  = engine.Binv_gpu.data() + delay_count;
        ptr_buffer[7][count_accepted]  = reinterpret_cast<T*>(engine.delay_list_gpu.data());
        ptr_buffer[8][count_accepted]  = engine.V_gpu.data() + norb * delay_count;
        ptr_buffer[9][count_accepted]  = const_cast<T*>(phi_vgl_v_dev_ptr + norb * iw);
        ptr_buffer[10][count_accepted] = const_cast<T*>(phi_vgl_v_dev_ptr + phi_vgl_stride + norb * 3 * iw);
        ptr_buffer[11][count_accepted] = const_cast<T*>(phi_vgl_v_dev_ptr + phi_vgl_stride * 4 + norb * iw);
        ptr_buffer[12][count_accepted] = psiM_g_list[count_accepted];
        ptr_buffer[13][count_accepted] = psiM_l_list[count_accepted];
        c_ratio_inv[count_accepted]    = T(1) / ratios[iw];
        count_accepted++;
      }
      else
      {
        ptr_buffer[0][n_accepted + count_rejected] = engine.psiMinv_.device_data() + lda * rowchanged;
        ptr_buffer[1][n_accepted + count_rejected] = engine.V_gpu.data();
        ptr_buffer[2][n_accepted + count_rejected] = engine.U_gpu.data() + norb * delay_count;
        ptr_buffer[3][n_accepted + count_rejected] = engine.p_gpu.data();
        ptr_buffer[4][n_accepted + count_rejected] = engine.Binv_gpu.data();
        ptr_buffer[5][n_accepted + count_rejected] = engine.Binv_gpu.data() + delay_count * lda_Binv;
        ptr_buffer[6][n_accepted + count_rejected] = engine.Binv_gpu.data() + delay_count;
        ptr_buffer[7][n_accepted + count_rejected] = reinterpret_cast<T*>(engine.delay_list_gpu.data());
        ptr_buffer[8][n_accepted + count_rejected] = engine.V_gpu.data() + norb * delay_count;
        count_rejected++;
      }
    }

    cudaErrorCheck(cudaMemcpyAsync(accept_rejectRow_buffer_H2D.device_data(), accept_rejectRow_buffer_H2D.data(),
                                   accept_rejectRow_buffer_H2D.size(), cudaMemcpyHostToDevice, hstream),
                   "cudaMemcpyAsync prepare_inv_row_buffer_H2D failed!");

    T** invRow_mw_ptr       = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D.device_data());
    T** V_mw_ptr            = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D.device_data() + sizeof(T*) * nw);
    T** U_row_mw_ptr        = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D.device_data() + sizeof(T*) * nw * 2);
    T** p_mw_ptr            = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D.device_data() + sizeof(T*) * nw * 3);
    T** Binv_mw_ptr         = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D.device_data() + sizeof(T*) * nw * 4);
    T** BinvRow_mw_ptr      = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D.device_data() + sizeof(T*) * nw * 5);
    T** BinvCol_mw_ptr      = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D.device_data() + sizeof(T*) * nw * 6);
    int** delay_list_mw_ptr = reinterpret_cast<int**>(accept_rejectRow_buffer_H2D.device_data() + sizeof(T*) * nw * 7);
    T** V_row_mw_ptr        = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D.device_data() + sizeof(T*) * nw * 8);
    T** phiV_mw_ptr         = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D.device_data() + sizeof(T*) * nw * 9);
    T** dpsiM_mw_in         = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D.device_data() + sizeof(T*) * nw * 10);
    T** d2psiM_mw_in        = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D.device_data() + sizeof(T*) * nw * 11);
    T** dpsiM_mw_out        = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D.device_data() + sizeof(T*) * nw * 12);
    T** d2psiM_mw_out       = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D.device_data() + sizeof(T*) * nw * 13);
    T* ratio_inv_mw_ptr     = reinterpret_cast<T*>(accept_rejectRow_buffer_H2D.device_data() + sizeof(T*) * nw * 14);

    //std::copy_n(Ainv[rowchanged], norb, V[delay_count]);
    cudaErrorCheck(cuBLAS_MFs::copy_batched(hstream, norb, invRow_mw_ptr, 1, V_row_mw_ptr, 1, nw),
                   "cuBLAS_MFs::copy_batched failed!");
    // handle accepted walkers
    // the new Binv is [[X Y] [Z sigma]]
    //BLAS::gemv('T', norb, delay_count + 1, cminusone, V.data(), norb, psiV.data(), 1, czero, p.data(), 1);
    cudaErrorCheck(cuBLAS_MFs::gemv_batched(hstream, 'T', norb, delay_count, cminusone_vec.device_data(), V_mw_ptr,
                                            norb, phiV_mw_ptr, 1, czero_vec.device_data(), p_mw_ptr, 1, n_accepted),
                   "cuBLAS_MFs::gemv_batched failed!");
    // Y
    //BLAS::gemv('T', delay_count, delay_count, sigma, Binv.data(), lda_Binv, p.data(), 1, czero, Binv.data() + delay_count,
    //           lda_Binv);
    cudaErrorCheck(cuBLAS_MFs::gemv_batched(hstream, 'T', delay_count, delay_count, ratio_inv_mw_ptr, Binv_mw_ptr,
                                            lda_Binv, p_mw_ptr, 1, czero_vec.device_data(), BinvCol_mw_ptr, lda_Binv,
                                            n_accepted),
                   "cuBLAS_MFs::gemv_batched failed!");
    // X
    //BLAS::ger(delay_count, delay_count, cone, Binv[delay_count], 1, Binv.data() + delay_count, lda_Binv,
    //          Binv.data(), lda_Binv);
    cudaErrorCheck(cuBLAS_MFs::ger_batched(hstream, delay_count, delay_count, cone_vec.device_data(), BinvRow_mw_ptr, 1,
                                           BinvCol_mw_ptr, lda_Binv, Binv_mw_ptr, lda_Binv, n_accepted),
                   "cuBLAS_MFs::ger_batched failed!");
    // sigma and Z
    cudaErrorCheck(CUDA::add_delay_list_save_sigma_VGL_batched(hstream, delay_list_mw_ptr, rowchanged, delay_count,
                                                               Binv_mw_ptr, lda_Binv, ratio_inv_mw_ptr, phiV_mw_ptr,
                                                               dpsiM_mw_in, d2psiM_mw_in, U_row_mw_ptr, dpsiM_mw_out,
                                                               d2psiM_mw_out, norb, n_accepted, nw),
                   "CUDA::add_delay_list_save_y_VGL_batched failed!");
    delay_count++;
    // update Ainv when maximal delay is reached
    if (delay_count == lda_Binv)
      mw_updateInvMat(engines);
  }

  /** update the full Ainv and reset delay_count
   * @param Ainv inverse matrix
   */
  static void mw_updateInvMat(const RefVectorWithLeader<This_t>& engines)
  {
    auto& engine_leader = engines.getLeader();
    int& delay_count    = engine_leader.delay_count;
    if (delay_count == 0)
      return;
    // update the inverse matrix
    auto& hstream              = engine_leader.cuda_handles_->hstream;
    auto& h_cublas             = engine_leader.cuda_handles_->h_cublas;
    auto& updateInv_buffer_H2D = engine_leader.mw_mem_->updateInv_buffer_H2D;
    const int norb             = engine_leader.psiMinv_.rows();
    const int lda              = engine_leader.psiMinv_.cols();
    const int nw               = engines.size();
    updateInv_buffer_H2D.resize(sizeof(T*) * 6 * nw);
    engine_leader.resize_fill_constant_arrays(nw);

    Matrix<T*> ptr_buffer(reinterpret_cast<T**>(updateInv_buffer_H2D.data()), 6, nw);
    for (int iw = 0; iw < nw; iw++)
    {
      This_t& engine    = engines[iw];
      ptr_buffer[0][iw] = engine.U_gpu.data();
      ptr_buffer[1][iw] = engine.psiMinv_.device_data();
      ptr_buffer[2][iw] = engine.tempMat_gpu.data();
      ptr_buffer[3][iw] = reinterpret_cast<T*>(engine.delay_list_gpu.data());
      ptr_buffer[4][iw] = engine.V_gpu.data();
      ptr_buffer[5][iw] = engine.Binv_gpu.data();
    }

    cudaErrorCheck(cudaMemcpyAsync(updateInv_buffer_H2D.device_data(), updateInv_buffer_H2D.data(),
                                   updateInv_buffer_H2D.size(), cudaMemcpyHostToDevice, hstream),
                   "cudaMemcpyAsync updateInv_buffer_H2D failed!");

    T** U_mw_ptr            = reinterpret_cast<T**>(updateInv_buffer_H2D.device_data());
    T** Ainv_mw_ptr         = reinterpret_cast<T**>(updateInv_buffer_H2D.device_data() + sizeof(T*) * nw);
    T** tempMat_mw_ptr      = reinterpret_cast<T**>(updateInv_buffer_H2D.device_data() + sizeof(T*) * nw * 2);
    int** delay_list_mw_ptr = reinterpret_cast<int**>(updateInv_buffer_H2D.device_data() + sizeof(T*) * nw * 3);
    T** V_mw_ptr            = reinterpret_cast<T**>(updateInv_buffer_H2D.device_data() + sizeof(T*) * nw * 4);
    T** Binv_mw_ptr         = reinterpret_cast<T**>(updateInv_buffer_H2D.device_data() + sizeof(T*) * nw * 5);

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
      constexpr T cone(1), czero(0), cminusone(-1);
      cublasErrorCheck(cuBLAS::gemm_batched(h_cublas, CUBLAS_OP_T, CUBLAS_OP_N, delay_count, norb, norb, &cone,
                                            U_mw_ptr, norb, Ainv_mw_ptr, lda, &czero, tempMat_mw_ptr, lda_Binv, nw),
                       "cuBLAS::gemm_batched failed!");
      cudaErrorCheck(CUDA::applyW_batched(hstream, delay_list_mw_ptr, delay_count, tempMat_mw_ptr, lda_Binv, nw),
                     "CUDA::applyW_batched failed!");
      cublasErrorCheck(cuBLAS::gemm_batched(h_cublas, CUBLAS_OP_N, CUBLAS_OP_N, norb, delay_count, delay_count, &cone,
                                            V_mw_ptr, norb, Binv_mw_ptr, lda_Binv, &czero, U_mw_ptr, norb, nw),
                       "cuBLAS::gemm_batched failed!");
      cublasErrorCheck(cuBLAS::gemm_batched(h_cublas, CUBLAS_OP_N, CUBLAS_OP_N, norb, norb, delay_count, &cminusone,
                                            U_mw_ptr, norb, tempMat_mw_ptr, lda_Binv, &cone, Ainv_mw_ptr, lda, nw),
                       "cuBLAS::gemm_batched failed!");
    }

    for (This_t& engine : engines)
    {
      cudaErrorCheck(cudaMemcpyAsync(engine.psiMinv_.data(), engine.psiMinv_.device_data(), engine.psiMinv_.size(),
                                     cudaMemcpyDeviceToHost, hstream),
                     "Failed to updated copy psiMinv back to host in mw_updateInvMat.");
    }
    delay_count = 0;
  }

  inline void print_Ainv(const RefVector<This_t>& engines)
  {
    for (This_t& engine : engines)
    {
      std::cout << "debug Ainv host  " << engine.psiMinv_[0][0] << " " << engine.psiMinv_[0][1] << " "
                << engine.psiMinv_[1][0] << " " << engine.psiMinv_[1][1] << std::endl;
      auto* temp_ptr = engine.psiMinv_.data();
      PRAGMA_OFFLOAD("omp target update from(temp_ptr[:engine.psiMinv_.size()])")
      std::cout << "debug Ainv devi  " << engine.psiMinv_[0][0] << " " << engine.psiMinv_[0][1] << " "
                << engine.psiMinv_[1][0] << " " << engine.psiMinv_[1][1] << std::endl;
    }
  }

  /** return invRow host or device pointers based on on_host request
   * prepare invRow if not already.
   */
  static std::vector<const T*> mw_getInvRow(const RefVectorWithLeader<This_t>& engines, const int row_id, bool on_host)
  {
    auto& engine_leader = engines.getLeader();
    if (engine_leader.isSM1())
      engine_leader.waitStream();
    else if (engine_leader.get_inv_row_id() != row_id)
    {
      // this can be skipped if mw_evalGrad gets called already.
      mw_prepareInvRow(engines, row_id);
      engine_leader.waitStream();
    }

    const size_t ncols = engines.getLeader().psiMinv_.cols();
    const size_t nw    = engines.size();
    std::vector<const T*> row_ptr_list;
    row_ptr_list.reserve(nw);
    if (on_host)
    {
      // copy values to host and return host pointer
      for (This_t& engine : engines)
        if (engine_leader.isSM1())
        {
          auto* ptr = engine.psiMinv_.data();
          PRAGMA_OFFLOAD("omp target update from(ptr[row_id * ncols : ncols])")
          row_ptr_list.push_back(ptr + row_id * ncols);
        }
        else
        {
          auto* ptr = engine.invRow.data();
          PRAGMA_OFFLOAD("omp target update from(ptr[:engine.invRow.size()])")
          row_ptr_list.push_back(ptr);
        }
    }
    else
    {
      // return device pointer
      for (This_t& engine : engines)
        if (engine_leader.isSM1())
          row_ptr_list.push_back(engine.psiMinv_.device_data() + row_id * ncols);
        else
          row_ptr_list.push_back(engine.invRow.device_data());
    }
    return row_ptr_list;
  }

  static void mw_transferAinv_D2H(const RefVectorWithLeader<This_t>& engines)
  {
    auto& engine_leader = engines.getLeader();
    auto& hstream       = engine_leader.cuda_handles_->hstream;
    engine_leader.guard_no_delay();

    for (This_t& engine : engines)
      cudaErrorCheck(cudaMemcpyAsync(engine.psiMinv_.data(), engine.psiMinv_.device_data(),
                                     engine.psiMinv_.size() * sizeof(T), cudaMemcpyDeviceToHost, hstream),
                     "cudaMemcpyAsync Ainv failed!");
    engine_leader.waitStream();
  }

  DiracMatrixComputeCUDA<T_FP>& get_det_inverter()
  {
    if (det_inverter_)
      return *det_inverter_;
    throw std::logic_error("attempted to get null det_inverter_, this is developer logic error");
  }

  friend class qmcplusplus::testing::DiracDeterminantBatchedTest;
};
} // namespace qmcplusplus

#endif // QMCPLUSPLUS_MATRIX_DELAYED_UPDATE_CUDA_H
