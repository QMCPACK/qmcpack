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

#ifndef QMCPLUSPLUS_MATRIX_DELAYED_UPDATE_CUDA_H
#define QMCPLUSPLUS_MATRIX_DELAYED_UPDATE_CUDA_H

#include "CPU/SIMD/aligned_allocator.hpp"
#include "Platforms/PinnedAllocator.h"
#include "OpenMP/OMPallocator.hpp"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "QMCWaveFunctions/Fermion/DiracMatrix.h"
#include "Platforms/OpenMP/ompBLAS.hpp"
#include <cuda_runtime_api.h>
#include "CUDA/cuBLAS.hpp"
#include "CUDA/cuBLAS_missing_functions.hpp"
#include "QMCWaveFunctions/detail/CUDA/matrix_update_helper.hpp"
#include "CUDA/CUDAallocator.hpp"

namespace qmcplusplus
{
/** implements dirac matrix delayed update using OpenMP offload and CUDA.
 * It is used as DET_ENGINE_TYPE in DiracDeterminantBatched.
 * @tparam T base precision for most computation
 * @tparam T_FP high precision for matrix inversion, T_FP >= T
 */
template<typename T, typename T_FP>
class MatrixDelayedUpdateCUDA
{
  using This_t = MatrixDelayedUpdateCUDA<T, T_FP>;

  template<typename DT>
  using OffloadAllocator = OMPallocator<DT, aligned_allocator<DT>>;
  template<typename DT>
  using OffloadPinnedAllocator     = OMPallocator<DT, PinnedAlignedAllocator<DT>>;
  using OffloadValueVector_t       = Vector<T, OffloadAllocator<T>>;
  using OffloadPinnedValueVector_t = Vector<T, OffloadPinnedAllocator<T>>;
  using OffloadPinnedValueMatrix_t = Matrix<T, OffloadPinnedAllocator<T>>;

  /// matrix inversion engine
  DiracMatrix<T_FP> detEng;
  /// inverse transpose of psiM(j,i) \f$= \psi_j({\bf r}_i)\f$
  OffloadPinnedValueMatrix_t psiMinv;
  /// device pointer of psiMinv data
  T* psiMinv_dev_ptr;
  /// scratch space for rank-1 update
  OffloadValueVector_t temp;
  /// device pointer of temp
  T* temp_dev_ptr;
  // row of up-to-date Ainv
  OffloadValueVector_t invRow;
  /** row id correspond to the up-to-date invRow. [0 norb), invRow is ready; -1, invRow is not valid.
   *  This id is set after calling getInvRow indicating invRow has been prepared for the invRow_id row
   *  ratioGrad checks if invRow_id is consistent. If not, invRow needs to be recomputed.
   *  acceptMove and completeUpdates mark invRow invalid by setting invRow_id to -1
   */
  int invRow_id;
  /// device pointer of invRow
  T* invRow_dev_ptr;
  // scratch space for keeping one row of Ainv
  OffloadValueVector_t rcopy;
  /// device pointer of rcopy
  T* rcopy_dev_ptr;
  // constant array value T(1)
  OffloadValueVector_t cone_vec;
  // device pointer of cone_vec
  T* cone_dev_ptr;
  // constant array value T(-1)
  OffloadValueVector_t cminusone_vec;
  // device pointer of cminusone_vec
  T* cminusone_dev_ptr;
  // constant array value T(0)
  OffloadValueVector_t czero_vec;
  // device pointer of czero_vec
  T* czero_dev_ptr;
  // multi walker of grads for transfer needs.
  OffloadPinnedValueMatrix_t grads_value_v;
  // device pointer of grads_value_v
  T* grads_value_dev_ptr;
  // mw_updateRow pointer buffer
  Vector<char, OffloadPinnedAllocator<char>> updateRow_buffer_H2D;
  // device pointer of updateRow_buffer_H2D
  char* updateRow_buffer_H2D_dev_ptr;
  // mw_prepareInvRow pointer buffer
  Vector<char, OffloadPinnedAllocator<char>> prepare_inv_row_buffer_H2D;
  // device pointer of prepare_inv_row_buffer_H2D
  char* prepare_inv_row_buffer_H2D_dev_ptr;
  // mw_accept_rejectRow pointer buffer
  Vector<char, OffloadPinnedAllocator<char>> accept_rejectRow_buffer_H2D;
  // device pointer of accept_rejectRow_buffer_H2D
  char* accept_rejectRow_buffer_H2D_dev_ptr;
  // mw_updateInv pointer buffer
  Vector<char, OffloadPinnedAllocator<char>> updateInv_buffer_H2D;
  // device pointer of updateInv_buffer_H2D
  char* updateInv_buffer_H2D_dev_ptr;

  // mw_evalGrad pointer buffer
  Vector<char, OffloadPinnedAllocator<char>> evalGrad_buffer_H2D;
  // device pointer of evalGrad_buffer_H2D
  char* evalGrad_buffer_H2D_dev_ptr;

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

  // CUDA specific variables
  cudaStream_t hstream;
  cublasHandle_t h_cublas;

  inline void waitStream() { cudaErrorCheck(cudaStreamSynchronize(hstream), "cudaStreamSynchronize failed!"); }
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
      cone_dev_ptr = getOffloadDevicePtr(cone_ptr);
      // cminusone
      cminusone_vec.resize(nw);
      std::fill_n(cminusone_vec.data(), nw, T(-1));
      T* cminusone_ptr = cminusone_vec.data();
      PRAGMA_OFFLOAD("omp target update to(cminusone_ptr[:nw])")
      cminusone_dev_ptr = getOffloadDevicePtr(cminusone_ptr);
      // czero
      czero_vec.resize(nw);
      std::fill_n(czero_vec.data(), nw, T(0));
      T* czero_ptr = czero_vec.data();
      PRAGMA_OFFLOAD("omp target update to(czero_ptr[:nw])")
      czero_dev_ptr = getOffloadDevicePtr(czero_ptr);
    }
  }

  void resize_prepareInvRow_scratch_arrays(size_t nw)
  {
    if (prepare_inv_row_buffer_H2D.size() < sizeof(T*) * 7 * nw)
    {
      prepare_inv_row_buffer_H2D.resize(sizeof(T*) * 7 * nw);
      prepare_inv_row_buffer_H2D_dev_ptr = getOffloadDevicePtr(prepare_inv_row_buffer_H2D.data());
    }
  }

  void resize_evalGrad_scratch_arrays(size_t nw)
  {
    if (evalGrad_buffer_H2D.size() < sizeof(T*) * 2 * nw)
    {
      evalGrad_buffer_H2D.resize(sizeof(T*) * 2 * nw);
      evalGrad_buffer_H2D_dev_ptr = getOffloadDevicePtr(evalGrad_buffer_H2D.data());
    }
  }

  void resize_updateRow_scratch_arrays(int norb, size_t nw)
  {
    size_t total_size = norb * nw;
    if (temp.size() < total_size)
    {
      temp.resize(total_size);
      temp_dev_ptr = getOffloadDevicePtr(temp.data());
      rcopy.resize(total_size);
      rcopy_dev_ptr = getOffloadDevicePtr(rcopy.data());
      updateRow_buffer_H2D.resize((sizeof(T*) * 8 + sizeof(T)) * nw);
      updateRow_buffer_H2D_dev_ptr = getOffloadDevicePtr(updateRow_buffer_H2D.data());
    }
  }

  void resize_accept_rejectRow_scratch_arrays(size_t nw)
  {
    const size_t total_size = (sizeof(T*) * 14 + sizeof(T)) * nw;
    if (accept_rejectRow_buffer_H2D.size() < total_size)
    {
      accept_rejectRow_buffer_H2D.resize(total_size);
      accept_rejectRow_buffer_H2D_dev_ptr = getOffloadDevicePtr(accept_rejectRow_buffer_H2D.data());
    }
  }

  void resize_updateInv_scratch_arrays(size_t nw)
  {
    const size_t total_size = sizeof(T*) * 6 * nw;
    if (updateInv_buffer_H2D.size() < total_size)
    {
      updateInv_buffer_H2D.resize(total_size);
      updateInv_buffer_H2D_dev_ptr = getOffloadDevicePtr(updateInv_buffer_H2D.data());
    }
  }

  void resizeGradsArray(size_t nw, size_t ndim)
  {
    if (grads_value_v.rows() != nw || grads_value_v.cols() != ndim)
    {
      grads_value_v.resize(nw, ndim);
      grads_value_dev_ptr = getOffloadDevicePtr(grads_value_v.data());
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
    resize_prepareInvRow_scratch_arrays(nw);
    resize_fill_constant_arrays(nw);

    const int lda_Binv = Binv_gpu.cols();
    Matrix<T*> ptr_buffer(reinterpret_cast<T**>(prepare_inv_row_buffer_H2D.data()), 7, nw);
    for (int iw = 0; iw < nw; iw++)
    {
      This_t& engine    = engines[iw];
      ptr_buffer[0][iw] = engine.psiMinv_dev_ptr + rowchanged * psiMinv.cols();
      ptr_buffer[1][iw] = engine.invRow_dev_ptr;
      ptr_buffer[2][iw] = engine.U_gpu.data();
      ptr_buffer[3][iw] = engine.p_gpu.data();
      ptr_buffer[4][iw] = engine.Binv_gpu.data();
      ptr_buffer[5][iw] = engine.Binv_gpu.data() + delay_count * lda_Binv;
      ptr_buffer[6][iw] = engine.V_gpu.data();
    }

    cudaErrorCheck(cudaMemcpyAsync(prepare_inv_row_buffer_H2D_dev_ptr, prepare_inv_row_buffer_H2D.data(),
                                   prepare_inv_row_buffer_H2D.size(), cudaMemcpyHostToDevice, hstream),
                   "cudaMemcpyAsync prepare_inv_row_buffer_H2D failed!");

    T** oldRow_mw_ptr  = reinterpret_cast<T**>(prepare_inv_row_buffer_H2D_dev_ptr);
    T** invRow_mw_ptr  = reinterpret_cast<T**>(prepare_inv_row_buffer_H2D_dev_ptr + sizeof(T*) * nw);
    T** U_mw_ptr       = reinterpret_cast<T**>(prepare_inv_row_buffer_H2D_dev_ptr + sizeof(T*) * nw * 2);
    T** p_mw_ptr       = reinterpret_cast<T**>(prepare_inv_row_buffer_H2D_dev_ptr + sizeof(T*) * nw * 3);
    T** Binv_mw_ptr    = reinterpret_cast<T**>(prepare_inv_row_buffer_H2D_dev_ptr + sizeof(T*) * nw * 4);
    T** BinvRow_mw_ptr = reinterpret_cast<T**>(prepare_inv_row_buffer_H2D_dev_ptr + sizeof(T*) * nw * 5);
    T** V_mw_ptr       = reinterpret_cast<T**>(prepare_inv_row_buffer_H2D_dev_ptr + sizeof(T*) * nw * 6);

    // save Ainv[rowchanged] to invRow
    //std::copy_n(Ainv[rowchanged], norb, invRow.data());
    cudaErrorCheck(cuBLAS_MFs::copy_batched(hstream, norb, oldRow_mw_ptr, 1, invRow_mw_ptr, 1, nw),
                   "cuBLAS_MFs::copy_batched failed!");
    // multiply V (NxK) Binv(KxK) U(KxN) invRow right to the left
    //BLAS::gemv('T', norb, delay_count, cone, U_gpu.data(), norb, invRow.data(), 1, czero, p_gpu.data(), 1);
    //BLAS::gemv('N', delay_count, delay_count, -cone, Binv.data(), lda_Binv, p.data(), 1, czero, Binv[delay_count], 1);
    //BLAS::gemv('N', norb, delay_count, cone, V.data(), norb, Binv[delay_count], 1, cone, invRow.data(), 1);
    cudaErrorCheck(cuBLAS_MFs::gemv_batched(hstream, 'T', norb, delay_count, cone_dev_ptr, U_mw_ptr, norb,
                                            invRow_mw_ptr, 1, czero_dev_ptr, p_mw_ptr, 1, nw),
                   "cuBLAS_MFs::gemv_batched failed!");
    cudaErrorCheck(cuBLAS_MFs::gemv_batched(hstream, 'N', delay_count, delay_count, cminusone_dev_ptr, Binv_mw_ptr,
                                            lda_Binv, p_mw_ptr, 1, czero_dev_ptr, BinvRow_mw_ptr, 1, nw),
                   "cuBLAS_MFs::gemv_batched failed!");
    cudaErrorCheck(cuBLAS_MFs::gemv_batched(hstream, 'N', norb, delay_count, cone_dev_ptr, V_mw_ptr, norb,
                                            BinvRow_mw_ptr, 1, cone_dev_ptr, invRow_mw_ptr, 1, nw),
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

    // to handle T** of Ainv, psi_v, temp, rcopy
    Matrix<T*> ptr_buffer(reinterpret_cast<T**>(updateRow_buffer_H2D.data()), 8, n_accepted);
    T* c_ratio_inv = reinterpret_cast<T*>(updateRow_buffer_H2D.data() + sizeof(T*) * 8 * n_accepted);
    for (int iw = 0, count = 0; iw < isAccepted.size(); iw++)
      if (isAccepted[iw])
      {
        ptr_buffer[0][count] = engines[iw].get().psiMinv_dev_ptr;
        ptr_buffer[1][count] = const_cast<T*>(phi_vgl_v_dev_ptr + norb * iw);
        ptr_buffer[2][count] = temp_dev_ptr + norb * count;
        ptr_buffer[3][count] = rcopy_dev_ptr + norb * count;
        ptr_buffer[4][count] = psiM_g_list[count];
        ptr_buffer[5][count] = psiM_l_list[count];
        ptr_buffer[6][count] = const_cast<T*>(phi_vgl_v_dev_ptr + phi_vgl_stride + norb * 3 * iw);
        ptr_buffer[7][count] = const_cast<T*>(phi_vgl_v_dev_ptr + phi_vgl_stride * 4 + norb * iw);

        c_ratio_inv[count] = T(-1) / ratios[iw];
        count++;
      }

    // update the inverse matrix
    resize_fill_constant_arrays(n_accepted);

    cudaErrorCheck(cudaMemcpyAsync(updateRow_buffer_H2D_dev_ptr, updateRow_buffer_H2D.data(),
                                   updateRow_buffer_H2D.size(), cudaMemcpyHostToDevice, hstream),
                   "cudaMemcpyAsync updateRow_buffer_H2D failed!");

    {
      T** Ainv_mw_ptr   = reinterpret_cast<T**>(updateRow_buffer_H2D_dev_ptr);
      T** phiV_mw_ptr   = reinterpret_cast<T**>(updateRow_buffer_H2D_dev_ptr + sizeof(T*) * n_accepted);
      T** temp_mw_ptr   = reinterpret_cast<T**>(updateRow_buffer_H2D_dev_ptr + sizeof(T*) * n_accepted * 2);
      T** rcopy_mw_ptr  = reinterpret_cast<T**>(updateRow_buffer_H2D_dev_ptr + sizeof(T*) * n_accepted * 3);
      T** dpsiM_mw_out  = reinterpret_cast<T**>(updateRow_buffer_H2D_dev_ptr + sizeof(T*) * n_accepted * 4);
      T** d2psiM_mw_out = reinterpret_cast<T**>(updateRow_buffer_H2D_dev_ptr + sizeof(T*) * n_accepted * 5);
      T** dpsiM_mw_in   = reinterpret_cast<T**>(updateRow_buffer_H2D_dev_ptr + sizeof(T*) * n_accepted * 6);
      T** d2psiM_mw_in  = reinterpret_cast<T**>(updateRow_buffer_H2D_dev_ptr + sizeof(T*) * n_accepted * 7);
      T* ratio_inv_mw   = reinterpret_cast<T*>(updateRow_buffer_H2D_dev_ptr + sizeof(T*) * n_accepted * 8);

      // invoke the Fahy's variant of Sherman-Morrison update.
      cudaErrorCheck(cuBLAS_MFs::gemv_batched(hstream, 'T', norb, norb, cone_dev_ptr, Ainv_mw_ptr, lda, phiV_mw_ptr, 1,
                                              czero_dev_ptr, temp_mw_ptr, 1, n_accepted),
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
  MatrixDelayedUpdateCUDA() : invRow_id(-1), delay_count(0)
  {
    cudaErrorCheck(cudaStreamCreate(&hstream), "cudaStreamCreate failed!");
    cublasErrorCheck(cublasCreate(&h_cublas), "cublasCreate failed!");
    cublasErrorCheck(cublasSetStream(h_cublas, hstream), "cublasSetStream failed!");
  }

  ~MatrixDelayedUpdateCUDA()
  {
    cublasErrorCheck(cublasDestroy(h_cublas), "cublasDestroy failed!");
    cudaErrorCheck(cudaStreamDestroy(hstream), "cudaStreamDestroy failed!");
  }

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
    invRow_dev_ptr = getOffloadDevicePtr(invRow.data());
    psiMinv.resize(norb, getAlignedSize<T>(norb));
    psiMinv_dev_ptr = getOffloadDevicePtr(psiMinv.data());
  }

  inline OffloadPinnedValueMatrix_t& get_psiMinv() { return psiMinv; }

  inline T* getRow_psiMinv_offload(int row_id) { return psiMinv_dev_ptr + row_id * psiMinv.cols(); }

  /** compute the inverse of the transpose of matrix A
   * @param logdetT orbital value matrix
   * @param Ainv inverse matrix
   */
  template<typename TREAL, typename OMPALLOC>
  inline void invert_transpose(const Matrix<T>& logdetT, Matrix<T, OMPALLOC>& Ainv, std::complex<TREAL>& LogValue)
  {
    guard_no_delay();
    Matrix<T> Ainv_host_view(Ainv.data(), Ainv.rows(), Ainv.cols());
    detEng.invert_transpose(logdetT, Ainv_host_view, LogValue);
    T* Ainv_ptr = Ainv.data();
    PRAGMA_OFFLOAD("omp target update to(Ainv_ptr[:Ainv.size()])")
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
    resize_evalGrad_scratch_arrays(nw);
    Matrix<const T*> ptr_buffer(reinterpret_cast<const T**>(evalGrad_buffer_H2D.data()), 2, nw);
    for (int iw = 0; iw < nw; iw++)
    {
      if (isSM1())
        ptr_buffer[0][iw] = engines[iw].get().psiMinv_dev_ptr + rowchanged * psiMinv.cols();
      else
        ptr_buffer[0][iw] = engines[iw].get().invRow_dev_ptr;
      ptr_buffer[1][iw] = dpsiM_row_list[iw];
    }

    cudaErrorCheck(cudaMemcpyAsync(evalGrad_buffer_H2D_dev_ptr, evalGrad_buffer_H2D.data(), evalGrad_buffer_H2D.size(),
                                   cudaMemcpyHostToDevice, hstream),
                   "cudaMemcpyAsync evalGrad_buffer_H2D failed!");

    resizeGradsArray(nw, GT::Size);

    const T** invRow_ptr    = reinterpret_cast<const T**>(evalGrad_buffer_H2D_dev_ptr);
    const T** dpsiM_row_ptr = reinterpret_cast<const T**>(evalGrad_buffer_H2D_dev_ptr) + nw;

    const int norb = psiMinv.rows();
    cudaErrorCheck(CUDA::calcGradients_cuda(hstream, norb, invRow_ptr, dpsiM_row_ptr, grads_value_dev_ptr, nw),
                   "CUDA::calcGradients_cuda failed!");
    cudaErrorCheck(cudaMemcpyAsync(grads_value_v.data(), grads_value_dev_ptr, grads_value_v.size() * sizeof(T),
                                   cudaMemcpyDeviceToHost, hstream),
                   "cudaMemcpyAsync grads_value_v failed!");
    waitStream();

    for (int iw = 0; iw < nw; iw++)
      grad_now[iw] = {grads_value_v[iw][0], grads_value_v[iw][1], grads_value_v[iw][2]};
  }

  template<typename VVT, typename RATIOT, typename OMPALLOC>
  void updateRow(Matrix<T, OMPALLOC>& Ainv, int rowchanged, const VVT& phiV, RATIOT c_ratio_in)
  {
    guard_no_delay();
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
    resize_accept_rejectRow_scratch_arrays(nw);
    resize_fill_constant_arrays(nw);

    Matrix<T*> ptr_buffer(reinterpret_cast<T**>(accept_rejectRow_buffer_H2D.data()), 14, nw);
    T* c_ratio_inv = reinterpret_cast<T*>(accept_rejectRow_buffer_H2D.data() + sizeof(T*) * 14 * nw);
    for (int iw = 0, count_accepted = 0, count_rejected = 0; iw < nw; iw++)
    {
      This_t& engine = engines[iw];
      if (isAccepted[iw])
      {
        ptr_buffer[0][count_accepted]  = engine.psiMinv_dev_ptr + lda * rowchanged;
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
        ptr_buffer[0][n_accepted + count_rejected] = engine.psiMinv_dev_ptr + lda * rowchanged;
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

    cudaErrorCheck(cudaMemcpyAsync(accept_rejectRow_buffer_H2D_dev_ptr, accept_rejectRow_buffer_H2D.data(),
                                   accept_rejectRow_buffer_H2D.size(), cudaMemcpyHostToDevice, hstream),
                   "cudaMemcpyAsync prepare_inv_row_buffer_H2D failed!");

    T** invRow_mw_ptr       = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D_dev_ptr);
    T** V_mw_ptr            = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D_dev_ptr + sizeof(T*) * nw);
    T** U_row_mw_ptr        = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D_dev_ptr + sizeof(T*) * nw * 2);
    T** p_mw_ptr            = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D_dev_ptr + sizeof(T*) * nw * 3);
    T** Binv_mw_ptr         = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D_dev_ptr + sizeof(T*) * nw * 4);
    T** BinvRow_mw_ptr      = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D_dev_ptr + sizeof(T*) * nw * 5);
    T** BinvCol_mw_ptr      = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D_dev_ptr + sizeof(T*) * nw * 6);
    int** delay_list_mw_ptr = reinterpret_cast<int**>(accept_rejectRow_buffer_H2D_dev_ptr + sizeof(T*) * nw * 7);
    T** V_row_mw_ptr        = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D_dev_ptr + sizeof(T*) * nw * 8);
    T** phiV_mw_ptr         = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D_dev_ptr + sizeof(T*) * nw * 9);
    T** dpsiM_mw_in         = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D_dev_ptr + sizeof(T*) * nw * 10);
    T** d2psiM_mw_in        = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D_dev_ptr + sizeof(T*) * nw * 11);
    T** dpsiM_mw_out        = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D_dev_ptr + sizeof(T*) * nw * 12);
    T** d2psiM_mw_out       = reinterpret_cast<T**>(accept_rejectRow_buffer_H2D_dev_ptr + sizeof(T*) * nw * 13);
    T* ratio_inv_mw_ptr     = reinterpret_cast<T*>(accept_rejectRow_buffer_H2D_dev_ptr + sizeof(T*) * nw * 14);

    //std::copy_n(Ainv[rowchanged], norb, V[delay_count]);
    cudaErrorCheck(cuBLAS_MFs::copy_batched(hstream, norb, invRow_mw_ptr, 1, V_row_mw_ptr, 1, nw),
                   "cuBLAS_MFs::copy_batched failed!");
    // handle accepted walkers
    // the new Binv is [[X Y] [Z y]]
    //BLAS::gemv('T', norb, delay_count + 1, cminusone, V.data(), norb, psiV.data(), 1, czero, p.data(), 1);
    cudaErrorCheck(cuBLAS_MFs::gemv_batched(hstream, 'T', norb, delay_count, cminusone_dev_ptr, V_mw_ptr, norb,
                                            phiV_mw_ptr, 1, czero_dev_ptr, p_mw_ptr, 1, n_accepted),
                   "cuBLAS_MFs::gemv_batched failed!");
    // Y
    //BLAS::gemv('T', delay_count, delay_count, y, Binv.data(), lda_Binv, p.data(), 1, czero, Binv.data() + delay_count,
    //           lda_Binv);
    cudaErrorCheck(cuBLAS_MFs::gemv_batched(hstream, 'T', delay_count, delay_count, ratio_inv_mw_ptr, Binv_mw_ptr,
                                            lda_Binv, p_mw_ptr, 1, czero_dev_ptr, BinvCol_mw_ptr, lda_Binv, n_accepted),
                   "cuBLAS_MFs::gemv_batched failed!");
    // X
    //BLAS::ger(delay_count, delay_count, cone, Binv[delay_count], 1, Binv.data() + delay_count, lda_Binv,
    //          Binv.data(), lda_Binv);
    cudaErrorCheck(cuBLAS_MFs::ger_batched(hstream, delay_count, delay_count, cone_dev_ptr, BinvRow_mw_ptr, 1,
                                           BinvCol_mw_ptr, lda_Binv, Binv_mw_ptr, lda_Binv, n_accepted),
                   "cuBLAS_MFs::ger_batched failed!");
    // y and Z
    cudaErrorCheck(CUDA::add_delay_list_save_y_VGL_batched(hstream, delay_list_mw_ptr, rowchanged, delay_count,
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
  void mw_updateInvMat(const RefVector<This_t>& engines)
  {
    if (delay_count == 0)
      return;
    // update the inverse matrix
    const int norb = psiMinv.rows();
    const int lda  = psiMinv.cols();
    const int nw   = engines.size();
    resize_updateInv_scratch_arrays(nw);
    resize_fill_constant_arrays(nw);

    Matrix<T*> ptr_buffer(reinterpret_cast<T**>(updateInv_buffer_H2D.data()), 6, nw);
    for (int iw = 0; iw < nw; iw++)
    {
      This_t& engine    = engines[iw];
      ptr_buffer[0][iw] = engine.U_gpu.data();
      ptr_buffer[1][iw] = engine.psiMinv_dev_ptr;
      ptr_buffer[2][iw] = engine.tempMat_gpu.data();
      ptr_buffer[3][iw] = reinterpret_cast<T*>(engine.delay_list_gpu.data());
      ptr_buffer[4][iw] = engine.V_gpu.data();
      ptr_buffer[5][iw] = engine.Binv_gpu.data();
    }

    cudaErrorCheck(cudaMemcpyAsync(updateInv_buffer_H2D_dev_ptr, updateInv_buffer_H2D.data(),
                                   updateInv_buffer_H2D.size(), cudaMemcpyHostToDevice, hstream),
                   "cudaMemcpyAsync updateInv_buffer_H2D failed!");

    T** U_mw_ptr            = reinterpret_cast<T**>(updateInv_buffer_H2D_dev_ptr);
    T** Ainv_mw_ptr         = reinterpret_cast<T**>(updateInv_buffer_H2D_dev_ptr + sizeof(T*) * nw);
    T** tempMat_mw_ptr      = reinterpret_cast<T**>(updateInv_buffer_H2D_dev_ptr + sizeof(T*) * nw * 2);
    int** delay_list_mw_ptr = reinterpret_cast<int**>(updateInv_buffer_H2D_dev_ptr + sizeof(T*) * nw * 3);
    T** V_mw_ptr            = reinterpret_cast<T**>(updateInv_buffer_H2D_dev_ptr + sizeof(T*) * nw * 4);
    T** Binv_mw_ptr         = reinterpret_cast<T**>(updateInv_buffer_H2D_dev_ptr + sizeof(T*) * nw * 5);

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
          row_ptr_list.push_back(engine.psiMinv_dev_ptr + row_id * psiMinv.cols());
        else
          row_ptr_list.push_back(engine.invRow_dev_ptr);
    }
    return row_ptr_list;
  }

  inline void mw_transferAinv_D2H(const RefVector<This_t>& engines)
  {
    guard_no_delay();

    for (This_t& engine : engines)
      cudaErrorCheck(cudaMemcpyAsync(engine.psiMinv.data(), engine.psiMinv_dev_ptr, engine.psiMinv.size() * sizeof(T),
                                     cudaMemcpyDeviceToHost, hstream),
                     "cudaMemcpyAsync Ainv failed!");
    waitStream();
  }
};
} // namespace qmcplusplus

#endif // QMCPLUSPLUS_MATRIX_DELAYED_UPDATE_CUDA_H
