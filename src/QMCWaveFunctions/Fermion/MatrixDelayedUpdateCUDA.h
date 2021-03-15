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

#include "CPU/SIMD/aligned_allocator.hpp"
#include "Platforms/PinnedAllocator.h"
#include "OMPTarget/OMPallocator.hpp"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Fermion/DiracMatrixComputeCUDA.hpp"
#include "Platforms/OMPTarget/ompBLAS.hpp"
#include <cuda_runtime_api.h>
#include "CUDA/cuBLAS.hpp"
#include "CUDA/cuBLAS_missing_functions.hpp"
#include "QMCWaveFunctions/detail/CUDA/matrix_update_helper.hpp"
#include "CUDA/CUDAallocator.hpp"
#include "Platforms/CUDA/CUDALinearAlgebraHandles.h"

namespace qmcplusplus
{

namespace testing
{
  class DiracDeterminantBatchedTest;
}
/** implements dirac matrix delayed update using OpenMP offload and CUDA.
 * It is used as DET_ENGINE_TYPE in DiracDeterminantBatched.
 * @tparam T base precision for most computation
 * @tparam T_FP high precision for matrix inversion, T_FP >= T
 */
template<typename T, typename T_FP>
class MatrixDelayedUpdateCUDA
{
public:
  using This_t = MatrixDelayedUpdateCUDA<T, T_FP>;

  template<typename DT>
  using OffloadAllocator = OMPallocator<DT, aligned_allocator<DT>>;
  template<typename DT>
  using OffloadPinnedAllocator     = OMPallocator<DT, PinnedAlignedAllocator<DT>>;
  using OffloadValueVector_t       = Vector<T, OffloadAllocator<T>>;
  using OffloadPinnedLogValueVector_t = Vector<std::complex<T>, OffloadPinnedAllocator<std::complex<T>>>;
  using OffloadPinnedValueVector_t = Vector<T, OffloadPinnedAllocator<T>>;
  using OffloadPinnedValueMatrix_t = Matrix<T, OffloadPinnedAllocator<T>>;

  using DiracMatrixCompute = DiracMatrixComputeCUDA<T_FP>;
  using Handles = CUDALinearAlgebraHandles;
private:
  /// inverse transpose of psiM(j,i) \f$= \psi_j({\bf r}_i)\f$
  OffloadPinnedValueMatrix_t psiMinv;
  /// scratch space for rank-1 update
  OffloadValueVector_t temp;
  // row of up-to-date Ainv
  OffloadValueVector_t invRow;
  /** row id correspond to the up-to-date invRow. [0 norb), invRow is ready; -1, invRow is not valid.
   *  This id is set after calling getInvRow indicating invRow has been prepared for the invRow_id row
   *  ratioGrad checks if invRow_id is consistent. If not, invRow needs to be recomputed.
   *  acceptMove and completeUpdates mark invRow invalid by setting invRow_id to -1
   */
  int invRow_id;
  // scratch space for keeping one row of Ainv
  OffloadValueVector_t rcopy;
  // constant array value T(1)
  OffloadValueVector_t cone_vec;
  // constant array value T(-1)
  OffloadValueVector_t cminusone_vec;
  // constant array value T(0)
  OffloadValueVector_t czero_vec;
  // multi walker of grads for transfer needs.
  OffloadPinnedValueMatrix_t grads_value_v;
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
  // mw_invert_transpose pointer buffer
  Vector<char, OffloadPinnedAllocator<char>> invert_transpose_buffer_H2D;
  // mw_invert transpose result pointer buffer
  Vector<char, OffloadPinnedAllocator<char>> invert_transpose_result_buffer_H2D;

  using DeviceValueMatrix_t = Matrix<T, CUDAAllocator<T>>;
  using DeviceValueVector_t = Vector<T, CUDAAllocator<T>>;
  /// orbital values of delayed electrons
  DeviceValueMatrix_t U_gpu;
  /// rows of Ainv corresponding to delayed electrons
  DeviceValueMatrix_t V_gpu;
  /// Matrix inverse of B, at maximum KxK
  DeviceValueMatrix_t Binv_gpu;
  /// scratch space, used during inverse update
  DeviceValueMatrix_t tempMat_gpu;
  /// new column of B
  DeviceValueVector_t p_gpu;
  /// list of delayed electrons
  Vector<int, CUDAAllocator<int>> delay_list_gpu;
  /// current number of delays, increase one for each acceptance, reset to 0 after updating Ainv
  int delay_count;

  /** @ingroup Resources
   *  @{ */
  // CUDA stream, cublas handle object
  std::unique_ptr<Handles> cuda_handles_;
  /// matrix inversion engine this a crowd scope resource and only the leader engine gets it
  UPtr<DiracMatrixCompute> det_inverter_;
  /**}@ */
  
  inline void waitStream()
  {
    cudaErrorCheck(cudaStreamSynchronize(cuda_handles_->hstream), "cudaStreamSynchronize failed!");
  }
  // ensure no previous delay left
  inline void guard_no_delay() const
  {
    if (delay_count != 0)
      throw std::runtime_error("BUG: unexpected call sequence delay_count is not 0");
  }

  // check if the number of maximal delay is 1 (SM-1)
  inline bool isSM1() const { return Binv_gpu.rows() == 1; }

  void resize_fill_constant_arrays(size_t nw)
  {
    if (cone_vec.size() < nw)
    {
      // cone
      cone_vec.resize(nw);
      std::fill_n(cone_vec.data(), nw, T(1));
      T* cone_ptr = cone_vec.data();
      PRAGMA_OFFLOAD("omp target update to(cone_ptr[:nw])")
      // cminusone
      cminusone_vec.resize(nw);
      std::fill_n(cminusone_vec.data(), nw, T(-1));
      T* cminusone_ptr = cminusone_vec.data();
      PRAGMA_OFFLOAD("omp target update to(cminusone_ptr[:nw])")
      // czero
      czero_vec.resize(nw);
      std::fill_n(czero_vec.data(), nw, T(0));
      T* czero_ptr = czero_vec.data();
      PRAGMA_OFFLOAD("omp target update to(czero_ptr[:nw])")
    }
  }

  void resize_updateRow_scratch_arrays(int norb, size_t nw)
  {
    size_t total_size = norb * nw;
    if (temp.size() < total_size)
    {
      temp.resize(total_size);
      rcopy.resize(total_size);
    }
  }

  /** compute the row of up-to-date Ainv
   * @param Ainv inverse matrix
   * @param rowchanged the row id corresponding to the proposed electron
   */
  void mw_prepareInvRow(const RefVector<This_t>& engines, const int rowchanged)
  {
    const int norb = psiMinv.rows();
    const int nw   = engines.size();
    prepare_inv_row_buffer_H2D.resize(sizeof(T*) * 7 * nw);
    resize_fill_constant_arrays(nw);

    const int lda_Binv = Binv_gpu.cols();
    Matrix<T*> ptr_buffer(reinterpret_cast<T**>(prepare_inv_row_buffer_H2D.data()), 7, nw);
    for (int iw = 0; iw < nw; iw++)
    {
      This_t& engine    = engines[iw];
      ptr_buffer[0][iw] = engine.psiMinv.device_data() + rowchanged * psiMinv.cols();
      ptr_buffer[1][iw] = engine.invRow.device_data();
      ptr_buffer[2][iw] = engine.U_gpu.data();
      ptr_buffer[3][iw] = engine.p_gpu.data();
      ptr_buffer[4][iw] = engine.Binv_gpu.data();
      ptr_buffer[5][iw] = engine.Binv_gpu.data() + delay_count * lda_Binv;
      ptr_buffer[6][iw] = engine.V_gpu.data();
    }

    cudaErrorCheck(cudaMemcpyAsync(prepare_inv_row_buffer_H2D.device_data(), prepare_inv_row_buffer_H2D.data(),
                                   prepare_inv_row_buffer_H2D.size(), cudaMemcpyHostToDevice, cuda_handles_->hstream),
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
    cudaErrorCheck(cuBLAS_MFs::copy_batched(cuda_handles_->hstream, norb, oldRow_mw_ptr, 1, invRow_mw_ptr, 1, nw),
                   "cuBLAS_MFs::copy_batched failed!");
    // multiply V (NxK) Binv(KxK) U(KxN) invRow right to the left
    //BLAS::gemv('T', norb, delay_count, cone, U_gpu.data(), norb, invRow.data(), 1, czero, p_gpu.data(), 1);
    //BLAS::gemv('N', delay_count, delay_count, -cone, Binv.data(), lda_Binv, p.data(), 1, czero, Binv[delay_count], 1);
    //BLAS::gemv('N', norb, delay_count, cone, V.data(), norb, Binv[delay_count], 1, cone, invRow.data(), 1);
    cudaErrorCheck(cuBLAS_MFs::gemv_batched(cuda_handles_->hstream, 'T', norb, delay_count, cone_vec.device_data(),
                                            U_mw_ptr, norb, invRow_mw_ptr, 1, czero_vec.device_data(), p_mw_ptr, 1, nw),
                   "cuBLAS_MFs::gemv_batched failed!");
    cudaErrorCheck(cuBLAS_MFs::gemv_batched(cuda_handles_->hstream, 'N', delay_count, delay_count,
                                            cminusone_vec.device_data(), Binv_mw_ptr, lda_Binv, p_mw_ptr, 1,
                                            czero_vec.device_data(), BinvRow_mw_ptr, 1, nw),
                   "cuBLAS_MFs::gemv_batched failed!");
    cudaErrorCheck(cuBLAS_MFs::gemv_batched(cuda_handles_->hstream, 'N', norb, delay_count, cone_vec.device_data(),
                                            V_mw_ptr, norb, BinvRow_mw_ptr, 1, cone_vec.device_data(), invRow_mw_ptr, 1,
                                            nw),
                   "cuBLAS_MFs::gemv_batched failed!");
    // mark row prepared
    invRow_id = rowchanged;
  }

  void mw_updateRow(const RefVector<This_t>& engines,
                    const int rowchanged,
                    const std::vector<T*>& psiM_g_list,
                    const std::vector<T*>& psiM_l_list,
                    const std::vector<bool>& isAccepted,
                    const T* phi_vgl_v_dev_ptr,
                    const size_t phi_vgl_stride,
                    const std::vector<T>& ratios)
  {
    guard_no_delay();

    const size_t n_accepted = psiM_g_list.size();
    if (n_accepted == 0)
      return;

    const int norb = psiMinv.rows();
    const int lda  = psiMinv.cols();
    resize_updateRow_scratch_arrays(norb, n_accepted);
    updateRow_buffer_H2D.resize((sizeof(T*) * 8 + sizeof(T)) * n_accepted);

    // to handle T** of Ainv, psi_v, temp, rcopy
    Matrix<T*> ptr_buffer(reinterpret_cast<T**>(updateRow_buffer_H2D.data()), 8, n_accepted);
    T* c_ratio_inv = reinterpret_cast<T*>(updateRow_buffer_H2D.data() + sizeof(T*) * 8 * n_accepted);
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
    resize_fill_constant_arrays(n_accepted);

    cudaErrorCheck(cudaMemcpyAsync(updateRow_buffer_H2D.device_data(), updateRow_buffer_H2D.data(),
                                   updateRow_buffer_H2D.size(), cudaMemcpyHostToDevice, cuda_handles_->hstream),
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
      cudaErrorCheck(cuBLAS_MFs::gemv_batched(cuda_handles_->hstream, 'T', norb, norb, cone_vec.device_data(),
                                              Ainv_mw_ptr, lda, phiV_mw_ptr, 1, czero_vec.device_data(), temp_mw_ptr, 1,
                                              n_accepted),
                     "cuBLAS_MFs::gemv_batched failed!");

      cudaErrorCheck(CUDA::copyAinvRow_saveGL_cuda(cuda_handles_->hstream, rowchanged, norb, Ainv_mw_ptr, lda,
                                                   temp_mw_ptr, rcopy_mw_ptr, dpsiM_mw_in, d2psiM_mw_in, dpsiM_mw_out,
                                                   d2psiM_mw_out, n_accepted),
                     "CUDA::copyAinvRow_saveGL_cuda failed!");


      cudaErrorCheck(cuBLAS_MFs::ger_batched(cuda_handles_->hstream, norb, norb, ratio_inv_mw, rcopy_mw_ptr, 1,
                                             temp_mw_ptr, 1, Ainv_mw_ptr, lda, n_accepted),
                     "cuBLAS_MFs::ger_batched failed!");
    }
  }

  Handles& getHandles() { return *cuda_handles_; }

public:
  /// default constructor
  MatrixDelayedUpdateCUDA() : invRow_id(-1), delay_count(0) {}

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
    psiMinv.resize(norb, getAlignedSize<T>(norb));
  }

  void createResource(ResourceCollection& collection)
  {
    auto resource_index = collection.addResource(std::make_unique<CUDALinearAlgebraHandles>());
    app_log() << "    CUDALinearAlgebraHandles resource created in MatrixDelayedUpdateCUDA. Index " << resource_index
              << std::endl;
    auto resource_index_det_eng = collection.addResource(std::make_unique<DiracMatrixComputeCUDA<T_FP>>());
    app_log() << "    DiracMatrixComputeCUDAed det_inverter_ created in MatrixDelayedUpdateCUDA. Index "
              << resource_index_det_eng << std::endl;
  }

  void acquireResource(ResourceCollection& collection)
  {
    auto res_ptr = dynamic_cast<CUDALinearAlgebraHandles*>(collection.lendResource().release());
    if (!res_ptr)
      throw std::runtime_error("MatrixDelayedUpdateCUDA::acquireResource dynamic_cast failed");
    cuda_handles_.reset(res_ptr);
    auto det_eng_ptr = dynamic_cast<DiracMatrixComputeCUDA<T_FP>*>(collection.lendResource().release());
    if (!det_eng_ptr)
      throw std::runtime_error(
          "MatrixDelayedUpdateCUDA::acquireResource dynamic_cast to DiracMatrixComputeCUDAed<T_FP>* failed");
    det_inverter_.reset(det_eng_ptr);
  }

  void releaseResource(ResourceCollection& collection)
  {
    collection.takebackResource(std::move(cuda_handles_));
    collection.takebackResource(std::move(det_inverter_));
  }

  inline OffloadPinnedValueMatrix_t& get_psiMinv() { return psiMinv; }

  inline T* getRow_psiMinv_offload(int row_id) { return psiMinv.device_data() + row_id * psiMinv.cols(); }

  /** make this class unit tests friendly without the need of setup resources.
   *  belongs in a friend class in test
   */
  inline void checkResourcesForTest()
  {    
    if (!cuda_handles_)
    {
      app_warning() << "MatrixDelayedUpdateCUDA local cuda_handles_ made : This message should not be seen in "
                       "production (performance bug) runs "
                       "but only unit tests (expected)."
                    << std::endl;
      cuda_handles_ = std::make_unique<CUDALinearAlgebraHandles>();
    }

    if (!det_inverter_)
    {
      app_warning() << "MatrixDelayedUpdateCUDA local det_inverter_ made : This message should not be seen in "
                       "production (performance bug) runs "
                       "but only unit tests (expected)."
                    << std::endl;

      det_inverter_ = std::make_unique<DiracMatrixComputeCUDA<T_FP>>();
    }
  }
    
  
  /** compute the inverse of the transpose of matrix logdetT, result is in psiMinv
   *
   *  This does not get called constantly so get real benchmark data that redirection to mw
   *  is a big deal before optimizing.
   * @param logdetT orbital value matrix
   * @param LogValue log(det(logdetT))
   */

  inline void invert_transpose(OffloadPinnedValueMatrix_t& log_det, OffloadPinnedLogValueVector_t& log_values)
  {
    checkResourcesForTest();
    guard_no_delay();
    det_inverter_->invert_transpose(*cuda_handles_, log_det, psiMinv, log_values);
    // RefVectorWithLeader<MatrixDelayedUpdateCUDA<T, T_FP>> engines(*this);
    // RefVector<OffloadPinnedValueMatrix_t> log_dets;
    // engines.push_back(*this);
    // log_dets.push_back(std::ref(log_det));
    // mw_invertTranspose(engines, log_dets, log_values);
  }

  inline void mw_invertTranspose(RefVectorWithLeader<MatrixDelayedUpdateCUDA<T, T_FP>>& engines,
                                 RefVector<OffloadPinnedValueMatrix_t>& logdetT_list,
                                 OffloadPinnedLogValueVector_t& log_values)
  {
    checkResourcesForTest();
    guard_no_delay();
    auto& engine_leader = engines.getLeader();

    RefVector<OffloadPinnedValueMatrix_t> a_inv_refs;
    a_inv_refs.reserve(engines.size());
    for (int iw = 0; iw < engines.size(); iw++)
    {
      a_inv_refs.emplace_back(engines[iw].psiMinv);
      T* a_inv_ptr = a_inv_refs.back().get().data();
      // This seems likely to be inefficient
      PRAGMA_OFFLOAD("omp target update to(a_inv_ptr[:a_inv_refs.back().get().size()])")
    }
    PRAGMA_OFFLOAD("omp taskwait")
      det_inverter_->mw_invertTranspose(*(engine_leader.cuda_handles_), logdetT_list, a_inv_refs, log_values);
  }

  // prepare invRow and compute the old gradients.
  template<typename GT>
  void mw_evalGrad(const RefVector<This_t>& engines,
                   const std::vector<const T*>& dpsiM_row_list,
                   const int rowchanged,
                   std::vector<GT>& grad_now)
  {
    if (!isSM1())
      mw_prepareInvRow(engines, rowchanged);

    const int nw = engines.size();
    evalGrad_buffer_H2D.resize(sizeof(T*) * 2 * nw);
    Matrix<const T*> ptr_buffer(reinterpret_cast<const T**>(evalGrad_buffer_H2D.data()), 2, nw);
    for (int iw = 0; iw < nw; iw++)
    {
      if (isSM1())
        ptr_buffer[0][iw] = engines[iw].get().psiMinv.device_data() + rowchanged * psiMinv.cols();
      else
        ptr_buffer[0][iw] = engines[iw].get().invRow.device_data();
      ptr_buffer[1][iw] = dpsiM_row_list[iw];
    }

    cudaErrorCheck(cudaMemcpyAsync(evalGrad_buffer_H2D.device_data(), evalGrad_buffer_H2D.data(),
                                   evalGrad_buffer_H2D.size(), cudaMemcpyHostToDevice, cuda_handles_->hstream),
                   "cudaMemcpyAsync evalGrad_buffer_H2D failed!");

    if (grads_value_v.rows() != nw || grads_value_v.cols() != GT::Size)
      grads_value_v.resize(nw, GT::Size);

    const T** invRow_ptr    = reinterpret_cast<const T**>(evalGrad_buffer_H2D.device_data());
    const T** dpsiM_row_ptr = reinterpret_cast<const T**>(evalGrad_buffer_H2D.device_data()) + nw;

    const int norb = psiMinv.rows();
    cudaErrorCheck(CUDA::calcGradients_cuda(cuda_handles_->hstream, norb, invRow_ptr, dpsiM_row_ptr,
                                            grads_value_v.device_data(), nw),
                   "CUDA::calcGradients_cuda failed!");
    cudaErrorCheck(cudaMemcpyAsync(grads_value_v.data(), grads_value_v.device_data(), grads_value_v.size() * sizeof(T),
                                   cudaMemcpyDeviceToHost, cuda_handles_->hstream),
                   "cudaMemcpyAsync grads_value_v failed!");
    waitStream();

    for (int iw = 0; iw < nw; iw++)
      grad_now[iw] = {grads_value_v[iw][0], grads_value_v[iw][1], grads_value_v[iw][2]};
  }

  template<typename VVT, typename RATIOT>
  void updateRow(int rowchanged, const VVT& phiV, RATIOT c_ratio_in)
  {
    guard_no_delay();
    auto& Ainv = psiMinv;
    // update the inverse matrix
    constexpr T cone(1), czero(0);
    const int norb = Ainv.rows();
    const int lda  = Ainv.cols();
    resize_updateRow_scratch_arrays(norb, 1);
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

  void mw_accept_rejectRow(const RefVector<This_t>& engines,
                           const int rowchanged,
                           const std::vector<T*>& psiM_g_list,
                           const std::vector<T*>& psiM_l_list,
                           const std::vector<bool>& isAccepted,
                           const T* phi_vgl_v_dev_ptr,
                           const size_t phi_vgl_stride,
                           const std::vector<T>& ratios)
  {
    // invRow consumed, mark invRow_id unset
    invRow_id = -1;

    if (isSM1())
    {
      mw_updateRow(engines, rowchanged, psiM_g_list, psiM_l_list, isAccepted, phi_vgl_v_dev_ptr, phi_vgl_stride,
                   ratios);
      return;
    }

    const int lda_Binv   = Binv_gpu.cols();
    const int norb       = psiMinv.rows();
    const int lda        = psiMinv.cols();
    const int nw         = engines.size();
    const int n_accepted = psiM_g_list.size();
    accept_rejectRow_buffer_H2D.resize((sizeof(T*) * 14 + sizeof(T)) * nw);
    resize_fill_constant_arrays(nw);

    Matrix<T*> ptr_buffer(reinterpret_cast<T**>(accept_rejectRow_buffer_H2D.data()), 14, nw);
    T* c_ratio_inv = reinterpret_cast<T*>(accept_rejectRow_buffer_H2D.data() + sizeof(T*) * 14 * nw);
    for (int iw = 0, count_accepted = 0, count_rejected = 0; iw < nw; iw++)
    {
      This_t& engine = engines[iw];
      if (isAccepted[iw])
      {
        ptr_buffer[0][count_accepted]  = engine.psiMinv.device_data() + lda * rowchanged;
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
        ptr_buffer[0][n_accepted + count_rejected] = engine.psiMinv.device_data() + lda * rowchanged;
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
                                   accept_rejectRow_buffer_H2D.size(), cudaMemcpyHostToDevice, cuda_handles_->hstream),
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
    cudaErrorCheck(cuBLAS_MFs::copy_batched(cuda_handles_->hstream, norb, invRow_mw_ptr, 1, V_row_mw_ptr, 1, nw),
                   "cuBLAS_MFs::copy_batched failed!");
    // handle accepted walkers
    // the new Binv is [[X Y] [Z sigma]]
    //BLAS::gemv('T', norb, delay_count + 1, cminusone, V.data(), norb, psiV.data(), 1, czero, p.data(), 1);
    cudaErrorCheck(cuBLAS_MFs::gemv_batched(cuda_handles_->hstream, 'T', norb, delay_count, cminusone_vec.device_data(),
                                            V_mw_ptr, norb, phiV_mw_ptr, 1, czero_vec.device_data(), p_mw_ptr, 1,
                                            n_accepted),
                   "cuBLAS_MFs::gemv_batched failed!");
    // Y
    //BLAS::gemv('T', delay_count, delay_count, sigma, Binv.data(), lda_Binv, p.data(), 1, czero, Binv.data() + delay_count,
    //           lda_Binv);
    cudaErrorCheck(cuBLAS_MFs::gemv_batched(cuda_handles_->hstream, 'T', delay_count, delay_count, ratio_inv_mw_ptr,
                                            Binv_mw_ptr, lda_Binv, p_mw_ptr, 1, czero_vec.device_data(), BinvCol_mw_ptr,
                                            lda_Binv, n_accepted),
                   "cuBLAS_MFs::gemv_batched failed!");
    // X
    //BLAS::ger(delay_count, delay_count, cone, Binv[delay_count], 1, Binv.data() + delay_count, lda_Binv,
    //          Binv.data(), lda_Binv);
    cudaErrorCheck(cuBLAS_MFs::ger_batched(cuda_handles_->hstream, delay_count, delay_count, cone_vec.device_data(),
                                           BinvRow_mw_ptr, 1, BinvCol_mw_ptr, lda_Binv, Binv_mw_ptr, lda_Binv,
                                           n_accepted),
                   "cuBLAS_MFs::ger_batched failed!");
    // sigma and Z
    cudaErrorCheck(CUDA::add_delay_list_save_sigma_VGL_batched(cuda_handles_->hstream, delay_list_mw_ptr, rowchanged,
                                                               delay_count, Binv_mw_ptr, lda_Binv, ratio_inv_mw_ptr,
                                                               phiV_mw_ptr, dpsiM_mw_in, d2psiM_mw_in, U_row_mw_ptr,
                                                               dpsiM_mw_out, d2psiM_mw_out, norb, n_accepted, nw),
                   "CUDA::add_delay_list_save_y_VGL_batched failed!");
    delay_count++;
    // update Ainv when maximal delay is reached
    if (delay_count == lda_Binv)
      mw_updateInvMat(engines);
  }

  /** update the full Ainv and reset delay_count
   * @param Ainv inverse matrix
   */
  void mw_updateInvMat(const RefVector<This_t>& engines)
  {
    if (delay_count == 0)
      return;
    // update the inverse matrix
    const int norb = psiMinv.rows();
    const int lda  = psiMinv.cols();
    const int nw   = engines.size();
    updateInv_buffer_H2D.resize(sizeof(T*) * 6 * nw);
    resize_fill_constant_arrays(nw);

    Matrix<T*> ptr_buffer(reinterpret_cast<T**>(updateInv_buffer_H2D.data()), 6, nw);
    for (int iw = 0; iw < nw; iw++)
    {
      This_t& engine    = engines[iw];
      ptr_buffer[0][iw] = engine.U_gpu.data();
      ptr_buffer[1][iw] = engine.psiMinv.device_data();
      ptr_buffer[2][iw] = engine.tempMat_gpu.data();
      ptr_buffer[3][iw] = reinterpret_cast<T*>(engine.delay_list_gpu.data());
      ptr_buffer[4][iw] = engine.V_gpu.data();
      ptr_buffer[5][iw] = engine.Binv_gpu.data();
    }

    cudaErrorCheck(cudaMemcpyAsync(updateInv_buffer_H2D.device_data(), updateInv_buffer_H2D.data(),
                                   updateInv_buffer_H2D.size(), cudaMemcpyHostToDevice, cuda_handles_->hstream),
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
      const int lda_Binv = Binv_gpu.cols();
      constexpr T cone(1), czero(0), cminusone(-1);
      cublasErrorCheck(cuBLAS::gemm_batched(cuda_handles_->h_cublas, CUBLAS_OP_T, CUBLAS_OP_N, delay_count, norb, norb,
                                            &cone, U_mw_ptr, norb, Ainv_mw_ptr, lda, &czero, tempMat_mw_ptr, lda_Binv,
                                            nw),
                       "cuBLAS::gemm_batched failed!");
      cudaErrorCheck(CUDA::applyW_batched(cuda_handles_->hstream, delay_list_mw_ptr, delay_count, tempMat_mw_ptr,
                                          lda_Binv, nw),
                     "CUDA::applyW_batched failed!");
      cublasErrorCheck(cuBLAS::gemm_batched(cuda_handles_->h_cublas, CUBLAS_OP_N, CUBLAS_OP_N, norb, delay_count,
                                            delay_count, &cone, V_mw_ptr, norb, Binv_mw_ptr, lda_Binv, &czero, U_mw_ptr,
                                            norb, nw),
                       "cuBLAS::gemm_batched failed!");
      cublasErrorCheck(cuBLAS::gemm_batched(cuda_handles_->h_cublas, CUBLAS_OP_N, CUBLAS_OP_N, norb, norb, delay_count,
                                            &cminusone, U_mw_ptr, norb, tempMat_mw_ptr, lda_Binv, &cone, Ainv_mw_ptr,
                                            lda, nw),
                       "cuBLAS::gemm_batched failed!");
    }
    delay_count = 0;
  }

  inline void print_Ainv(const RefVector<This_t>& engines)
  {
    for (This_t& engine : engines)
    {
      std::cout << "debug Ainv host  " << engine.psiMinv[0][0] << " " << engine.psiMinv[0][1] << " "
                << engine.psiMinv[1][0] << " " << engine.psiMinv[1][1] << std::endl;
      auto* temp_ptr = engine.psiMinv.data();
      PRAGMA_OFFLOAD("omp target update from(temp_ptr[:psiMinv.size()])")
      std::cout << "debug Ainv devi  " << engine.psiMinv[0][0] << " " << engine.psiMinv[0][1] << " "
                << engine.psiMinv[1][0] << " " << engine.psiMinv[1][1] << std::endl;
    }
  }

  /** return invRow host or device pointers based on on_host request
   * prepare invRow if not already.
   */
  std::vector<const T*> mw_getInvRow(const RefVector<This_t>& engines, const int row_id, bool on_host)
  {
    if (isSM1())
      waitStream();
    else if (invRow_id != row_id)
    {
      // this can be skipped if mw_evalGrad gets called already.
      mw_prepareInvRow(engines, row_id);
      waitStream();
    }

    const size_t nw = engines.size();
    std::vector<const T*> row_ptr_list;
    row_ptr_list.reserve(nw);
    if (on_host)
    {
      // copy values to host and return host pointer
      for (This_t& engine : engines)
        if (isSM1())
        {
          auto* ptr = engine.psiMinv.data();
          PRAGMA_OFFLOAD("omp target update from(ptr[row_id * psiMinv.cols():psiMinv.cols()])")
          row_ptr_list.push_back(ptr + row_id * psiMinv.cols());
        }
        else
        {
          auto* ptr = engine.invRow.data();
          PRAGMA_OFFLOAD("omp target update from(ptr[:invRow.size()])")
          row_ptr_list.push_back(ptr);
        }
    }
    else
    {
      // return device pointer
      for (This_t& engine : engines)
        if (isSM1())
          row_ptr_list.push_back(engine.psiMinv.device_data() + row_id * psiMinv.cols());
        else
          row_ptr_list.push_back(engine.invRow.device_data());
    }
    return row_ptr_list;
  }

  inline void mw_transferAinv_D2H(const RefVector<This_t>& engines)
  {
    guard_no_delay();

    for (This_t& engine : engines)
      cudaErrorCheck(cudaMemcpyAsync(engine.psiMinv.data(), engine.psiMinv.device_data(),
                                     engine.psiMinv.size() * sizeof(T), cudaMemcpyDeviceToHost, cuda_handles_->hstream),
                     "cudaMemcpyAsync Ainv failed!");
    waitStream();
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
