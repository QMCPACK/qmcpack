//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_CUSOLVERINVERTOR_H
#define QMCPLUSPLUS_CUSOLVERINVERTOR_H

#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "CUDA/CUDAruntime.hpp"
#include "CUDA/CUDAallocator.hpp"
#include "CUDA/cusolver.hpp"
#include "QMCWaveFunctions/detail/CUDA/delayed_update_helper.h"

namespace qmcplusplus
{
/** implements matrix inversion via cuSolverDN
 * @tparam T_FP high precision for matrix inversion, T_FP >= T
 */
template<typename T_FP>
class cuSolverInverter
{
  /// scratch memory for cusolverDN
  Matrix<T_FP, CUDAAllocator<T_FP>> Mat1_gpu;
  /// scratch memory for cusolverDN
  Matrix<T_FP, CUDAAllocator<T_FP>> Mat2_gpu;
  /// pivot array + info
  Vector<int, CUDAHostAllocator<int>> ipiv;
  Vector<int, CUDAAllocator<int>> ipiv_gpu;
  /// diagonal terms of LU matrix
  Vector<T_FP, CUDAHostAllocator<T_FP>> LU_diag;
  Vector<T_FP, CUDAAllocator<T_FP>> LU_diag_gpu;
  /// workspace
  Vector<T_FP, CUDAAllocator<T_FP>> work_gpu;

  // CUDA specific variables
  cusolverDnHandle_t h_cusolver_;
  cudaStream_t hstream_;

  /** resize the internal storage
   * @param norb number of electrons/orbitals
   * @param delay, maximum delay 0<delay<=norb
   */
  inline void resize(int norb)
  {
    if (Mat1_gpu.rows() != norb)
    {
      Mat1_gpu.resize(norb, norb);
      // prepare cusolver auxiliary arrays
      ipiv.resize(norb + 1);
      ipiv_gpu.resize(norb + 1);
      LU_diag.resize(norb);
      LU_diag_gpu.resize(norb);
      int lwork;
      cusolverErrorCheck(cusolver::getrf_bufferSize(h_cusolver_, norb, norb, Mat1_gpu.data(), norb, &lwork),
                         "cusolver::getrf_bufferSize failed!");
      work_gpu.resize(lwork);
    }
  }

public:
  /// default constructor
  cuSolverInverter()
  {
    cudaErrorCheck(cudaStreamCreate(&hstream_), "cudaStreamCreate failed!");
    cusolverErrorCheck(cusolverDnCreate(&h_cusolver_), "cusolverCreate failed!");
    cusolverErrorCheck(cusolverDnSetStream(h_cusolver_, hstream_), "cusolverSetStream failed!");
  }

  ~cuSolverInverter()
  {
    cusolverErrorCheck(cusolverDnDestroy(h_cusolver_), "cusolverDestroy failed!");
    cudaErrorCheck(cudaStreamDestroy(hstream_), "cudaStreamDestroy failed!");
  }

  /** compute the inverse of the transpose of matrix A and its determinant value in log
   * when T_FP and TMAT are the same
   * @tparam TREAL real type
   */
  template<typename TMAT, typename TREAL, typename = std::enable_if_t<std::is_same<TMAT, T_FP>::value>>
  std::enable_if_t<std::is_same<TMAT, T_FP>::value> invert_transpose(const Matrix<TMAT>& logdetT,
                                                                     Matrix<TMAT>& Ainv,
                                                                     Matrix<TMAT, CUDAAllocator<TMAT>>& Ainv_gpu,
                                                                     std::complex<TREAL>& log_value)
  {
    const int norb = logdetT.rows();
    resize(norb);
    cudaErrorCheck(cudaMemcpyAsync(Mat1_gpu.data(), logdetT.data(), logdetT.size() * sizeof(TMAT),
                                   cudaMemcpyHostToDevice, hstream_),
                   "cudaMemcpyAsync failed!");
    cusolverErrorCheck(cusolver::getrf(h_cusolver_, norb, norb, Mat1_gpu.data(), norb, work_gpu.data(),
                                       ipiv_gpu.data() + 1, ipiv_gpu.data()),
                       "cusolver::getrf failed!");
    cudaErrorCheck(cudaMemcpyAsync(ipiv.data(), ipiv_gpu.data(), ipiv_gpu.size() * sizeof(int), cudaMemcpyDeviceToHost,
                                   hstream_),
                   "cudaMemcpyAsync failed!");
    extract_matrix_diagonal_cuda(norb, Mat1_gpu.data(), norb, LU_diag_gpu.data(), hstream_);
    cudaErrorCheck(cudaMemcpyAsync(LU_diag.data(), LU_diag_gpu.data(), LU_diag.size() * sizeof(T_FP),
                                   cudaMemcpyDeviceToHost, hstream_),
                   "cudaMemcpyAsync failed!");
    // check LU success
    cudaErrorCheck(cudaStreamSynchronize(hstream_), "cudaStreamSynchronize failed!");
    if (ipiv[0] != 0)
    {
      std::ostringstream err;
      err << "cusolver::getrf calculation failed with devInfo = " << ipiv[0] << std::endl;
      std::cerr << err.str();
      throw std::runtime_error(err.str());
    }
    make_identity_matrix_cuda(norb, Ainv_gpu.data(), norb, hstream_);
    cusolverErrorCheck(cusolver::getrs(h_cusolver_, CUBLAS_OP_T, norb, norb, Mat1_gpu.data(), norb, ipiv_gpu.data() + 1,
                                       Ainv_gpu.data(), norb, ipiv_gpu.data()),
                       "cusolver::getrs failed!");
    cudaErrorCheck(cudaMemcpyAsync(ipiv.data(), ipiv_gpu.data(), sizeof(int), cudaMemcpyDeviceToHost, hstream_),
                   "cudaMemcpyAsync failed!");
    computeLogDet(LU_diag.data(), norb, ipiv.data() + 1, log_value);
    cudaErrorCheck(cudaMemcpyAsync(Ainv.data(), Ainv_gpu.data(), Ainv.size() * sizeof(TMAT), cudaMemcpyDeviceToHost,
                                   hstream_),
                   "cudaMemcpyAsync failed!");
    // no need to wait because : For transfers from device memory to pageable host memory, the function will return only once the copy has completed.
    if (ipiv[0] != 0)
    {
      std::ostringstream err;
      err << "cusolver::getrs calculation failed with devInfo = " << ipiv[0] << std::endl;
      std::cerr << err.str();
      throw std::runtime_error(err.str());
    }
  }

  /** compute the inverse of the transpose of matrix A and its determinant value in log
   * when T_FP and TMAT are not the same
   * @tparam TREAL real type
   */
  template<typename TMAT, typename TREAL, typename = std::enable_if_t<!std::is_same<TMAT, T_FP>::value>>
  std::enable_if_t<!std::is_same<TMAT, T_FP>::value> invert_transpose(const Matrix<TMAT>& logdetT,
                                                                      Matrix<TMAT>& Ainv,
                                                                      Matrix<TMAT, CUDAAllocator<TMAT>>& Ainv_gpu,
                                                                      std::complex<TREAL>& log_value)
  {
    const int norb = logdetT.rows();
    resize(norb);
    Mat2_gpu.resize(norb, norb);
    cudaErrorCheck(cudaMemcpyAsync(Mat2_gpu.data(), logdetT.data(), logdetT.size() * sizeof(TMAT),
                                   cudaMemcpyHostToDevice, hstream_),
                   "cudaMemcpyAsync failed!");
    copy_matrix_cuda(norb, norb, (TMAT*)Mat2_gpu.data(), norb, Mat1_gpu.data(), norb, hstream_);
    cusolverErrorCheck(cusolver::getrf(h_cusolver_, norb, norb, Mat1_gpu.data(), norb, work_gpu.data(),
                                       ipiv_gpu.data() + 1, ipiv_gpu.data()),
                       "cusolver::getrf failed!");
    cudaErrorCheck(cudaMemcpyAsync(ipiv.data(), ipiv_gpu.data(), ipiv_gpu.size() * sizeof(int), cudaMemcpyDeviceToHost,
                                   hstream_),
                   "cudaMemcpyAsync failed!");
    extract_matrix_diagonal_cuda(norb, Mat1_gpu.data(), norb, LU_diag_gpu.data(), hstream_);
    cudaErrorCheck(cudaMemcpyAsync(LU_diag.data(), LU_diag_gpu.data(), LU_diag.size() * sizeof(T_FP),
                                   cudaMemcpyDeviceToHost, hstream_),
                   "cudaMemcpyAsync failed!");
    // check LU success
    cudaErrorCheck(cudaStreamSynchronize(hstream_), "cudaStreamSynchronize failed!");
    if (ipiv[0] != 0)
    {
      std::ostringstream err;
      err << "cusolver::getrf calculation failed with devInfo = " << ipiv[0] << std::endl;
      std::cerr << err.str();
      throw std::runtime_error(err.str());
    }
    make_identity_matrix_cuda(norb, Mat2_gpu.data(), norb, hstream_);
    cusolverErrorCheck(cusolver::getrs(h_cusolver_, CUBLAS_OP_T, norb, norb, Mat1_gpu.data(), norb, ipiv_gpu.data() + 1,
                                       Mat2_gpu.data(), norb, ipiv_gpu.data()),
                       "cusolver::getrs failed!");
    copy_matrix_cuda(norb, norb, Mat2_gpu.data(), norb, Ainv_gpu.data(), norb, hstream_);
    cudaErrorCheck(cudaMemcpyAsync(ipiv.data(), ipiv_gpu.data(), sizeof(int), cudaMemcpyDeviceToHost, hstream_),
                   "cudaMemcpyAsync failed!");
    computeLogDet(LU_diag.data(), norb, ipiv.data() + 1, log_value);
    cudaErrorCheck(cudaMemcpyAsync(Ainv.data(), Ainv_gpu.data(), Ainv.size() * sizeof(TMAT), cudaMemcpyDeviceToHost,
                                   hstream_),
                   "cudaMemcpyAsync failed!");
    // no need to wait because : For transfers from device memory to pageable host memory, the function will return only once the copy has completed.
    if (ipiv[0] != 0)
    {
      std::ostringstream err;
      err << "cusolver::getrs calculation failed with devInfo = " << ipiv[0] << std::endl;
      std::cerr << err.str();
      throw std::runtime_error(err.str());
    }
  }
};
} // namespace qmcplusplus

#endif // QMCPLUSPLUS_CUSOLVERINVERTOR_H
