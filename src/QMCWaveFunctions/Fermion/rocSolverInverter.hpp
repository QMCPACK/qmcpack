//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_ROCSOLVERINVERTER_H
#define QMCPLUSPLUS_ROCSOLVERINVERTER_H

#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"

#if !defined(QMC_CUDA2HIP)
#error rocSolverInverter.hpp expects QMC_CUDA2HIP to be defined
#endif
// This file assumes that QMC_CUDA2HIP is defined and that creates HIP versions of these functions (despite being labeled with "CUDA")
#include "CUDA/CUDAruntime.hpp"
#include "CUDA/CUDAallocator.hpp"
#include "ROCm/rocsolver.hpp"
#include "QMCWaveFunctions/detail/CUDA/delayed_update_helper.h"

namespace qmcplusplus
{
/** implements matrix inversion via rocSolver
 * @tparam T_FP high precision for matrix inversion, T_FP >= T
 */
template<typename T_FP>
class rocSolverInverter
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
  rocblas_handle h_rocsolver_;
  hipStream_t hstream_;

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

      // Memory for temporary storage for solver calls.
      // The rocSOLVER library handles this memory itself.
      //  If we need more control, there are API's to get the size and set the buffer memory
#if 0
      size_t memory_size;
      rocblas_start_device_memory_size_query(h_rocsolver_);
      rocsolverErrorCheck(rocsolver::dgetrf(h_rocsolver_, norb, norb, nullptr, norb, nullptr, nullptr);
      rocsolverErrorCheck(rocsolver::dgetri(h_rocsolver_, norb, norb, nullptr, norb, nullptr, nullptr);
      rocblas_stop_device_memory_size_query(h_rocsolver_, &memory_size);
#endif
    }
  }

public:
  /// default constructor
  rocSolverInverter()
  {
    cudaErrorCheck(hipStreamCreate(&hstream_), "hipStreamCreate failed!");
    rocsolverErrorCheck(rocblas_create_handle(&h_rocsolver_), "rocblas_create_handle failed!");
    rocsolverErrorCheck(rocblas_set_stream(h_rocsolver_, hstream_), "rocblas_set_stream failed!");
  }

  ~rocSolverInverter()
  {
    rocsolverErrorCheck(rocblas_destroy_handle(h_rocsolver_), "rocblas_destroy_handle failed!");
    cudaErrorCheck(hipStreamDestroy(hstream_), "hipStreamDestroy failed!");
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
    cudaErrorCheck(hipMemcpyAsync(Mat1_gpu.data(), logdetT.data(), logdetT.size() * sizeof(TMAT), hipMemcpyHostToDevice,
                                  hstream_),
                   "hipMemcpyAsync for logdetT to Mat1_gpu failed!");
    rocsolverErrorCheck(rocsolver::getrf(h_rocsolver_, norb, norb, Mat1_gpu.data(), norb, ipiv_gpu.data() + 1,
                                         ipiv_gpu.data()),
                        "rocsolver::getrf failed!");
    cudaErrorCheck(hipMemcpyAsync(ipiv.data(), ipiv_gpu.data(), ipiv_gpu.size() * sizeof(int), hipMemcpyDeviceToHost,
                                  hstream_),
                   "hipMemcpyAsync for ipiv failed!");
    extract_matrix_diagonal_cuda(norb, Mat1_gpu.data(), norb, LU_diag_gpu.data(), hstream_);
    cudaErrorCheck(hipMemcpyAsync(LU_diag.data(), LU_diag_gpu.data(), LU_diag.size() * sizeof(T_FP),
                                  hipMemcpyDeviceToHost, hstream_),
                   "hipMemcpyAsync for LU_diag failed!");
    // check LU success
    cudaErrorCheck(hipStreamSynchronize(hstream_), "hipStreamSynchronize failed!");
    if (ipiv[0] != 0)
    {
      std::ostringstream err;
      err << "rocsolver::getrf calculation failed with devInfo = " << ipiv[0] << std::endl;
      std::cerr << err.str();
      throw std::runtime_error(err.str());
    }
    make_identity_matrix_cuda(norb, Ainv_gpu.data(), norb, hstream_);
    rocsolverErrorCheck(rocsolver::getrs(h_rocsolver_, rocblas_operation_transpose, norb, norb, Mat1_gpu.data(), norb,
                                         ipiv_gpu.data() + 1, Ainv_gpu.data(), norb),
                        "rocsolver::getrs failed!");
    cudaErrorCheck(hipMemcpyAsync(ipiv.data(), ipiv_gpu.data(), sizeof(int), hipMemcpyDeviceToHost, hstream_),
                   "hipMemcpyAsync for ipiv failed!");
    computeLogDet(LU_diag.data(), norb, ipiv.data() + 1, log_value);
    cudaErrorCheck(hipMemcpyAsync(Ainv.data(), Ainv_gpu.data(), Ainv.size() * sizeof(TMAT), hipMemcpyDeviceToHost,
                                  hstream_),
                   "hipMemcpyAsync for Ainv failed!");
    // no need to wait because : For transfers from device memory to pageable host memory, the function will return only once the copy has completed.
    if (ipiv[0] != 0)
    {
      std::ostringstream err;
      err << "rocsolver::getrs calculation failed with devInfo = " << ipiv[0] << std::endl;
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
    cudaErrorCheck(hipMemcpyAsync(Mat2_gpu.data(), logdetT.data(), logdetT.size() * sizeof(TMAT), hipMemcpyHostToDevice,
                                  hstream_),
                   "hipMemcpyAsync failed!");
    copy_matrix_cuda(norb, norb, (TMAT*)Mat2_gpu.data(), norb, Mat1_gpu.data(), norb, hstream_);
    rocsolverErrorCheck(rocsolver::getrf(h_rocsolver_, norb, norb, Mat1_gpu.data(), norb, ipiv_gpu.data() + 1,
                                         ipiv_gpu.data()),
                        "rocsolver::getrf failed!");
    cudaErrorCheck(hipMemcpyAsync(ipiv.data(), ipiv_gpu.data(), ipiv_gpu.size() * sizeof(int), hipMemcpyDeviceToHost,
                                  hstream_),
                   "hipMemcpyAsync failed!");
    extract_matrix_diagonal_cuda(norb, Mat1_gpu.data(), norb, LU_diag_gpu.data(), hstream_);
    cudaErrorCheck(hipMemcpyAsync(LU_diag.data(), LU_diag_gpu.data(), LU_diag.size() * sizeof(T_FP),
                                  hipMemcpyDeviceToHost, hstream_),
                   "hipMemcpyAsync failed!");
    // check LU success
    cudaErrorCheck(hipStreamSynchronize(hstream_), "hipStreamSynchronize failed!");
    if (ipiv[0] != 0)
    {
      std::ostringstream err;
      err << "rocsolver::getrf calculation failed with devInfo = " << ipiv[0] << std::endl;
      std::cerr << err.str();
      throw std::runtime_error(err.str());
    }
    make_identity_matrix_cuda(norb, Mat2_gpu.data(), norb, hstream_);
    rocsolverErrorCheck(rocsolver::getrs(h_rocsolver_, rocblas_operation_transpose, norb, norb, Mat1_gpu.data(), norb,
                                         ipiv_gpu.data() + 1, Mat2_gpu.data(), norb),
                        "rocsolver::getrs failed!");
    copy_matrix_cuda(norb, norb, Mat2_gpu.data(), norb, Ainv_gpu.data(), norb, hstream_);
    cudaErrorCheck(hipMemcpyAsync(ipiv.data(), ipiv_gpu.data(), sizeof(int), hipMemcpyDeviceToHost, hstream_),
                   "hipMemcpyAsync failed!");
    computeLogDet(LU_diag.data(), norb, ipiv.data() + 1, log_value);
    cudaErrorCheck(hipMemcpyAsync(Ainv.data(), Ainv_gpu.data(), Ainv.size() * sizeof(TMAT), hipMemcpyDeviceToHost,
                                  hstream_),
                   "hipMemcpyAsync failed!");
    // no need to wait because : For transfers from device memory to pageable host memory, the function will return only once the copy has completed.
    if (ipiv[0] != 0)
    {
      std::ostringstream err;
      err << "rocsolver::getrs calculation failed with devInfo = " << ipiv[0] << std::endl;
      std::cerr << err.str();
      throw std::runtime_error(err.str());
    }
  }
};
} // namespace qmcplusplus

#endif // QMCPLUSPLUS_ROCSOLVERINVERTER_H
