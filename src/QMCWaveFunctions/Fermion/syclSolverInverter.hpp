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

#ifndef QMCPLUSPLUS_SYCL_MKL_SOLVERINVERTOR_H
#define QMCPLUSPLUS_SYCL_MKL_SOLVERINVERTOR_H

#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "SYCL/SYCLallocator.hpp"
#include "SYCL/syclBLAS.hpp"
#include "SYCL/syclSolver.hpp"
#include "QMCWaveFunctions/detail/SYCL/sycl_determinant_helper.hpp"

namespace qmcplusplus
{
/** implements matrix inversion via cuSolverDN
 * @tparam T_FP high precision for matrix inversion, T_FP >= T
 */
template<typename T_FP>
class syclSolverInverter
{
  /// scratch memory for cusolverDN
  Matrix<T_FP, SYCLAllocator<T_FP>> Mat1_gpu;
  /// pivot array + info
  Vector<std::int64_t, SYCLAllocator<std::int64_t>> ipiv;
  /// workspace
  Vector<T_FP, SYCLAllocator<T_FP>> workspace;
  std::int64_t getrf_ws = 0;
  std::int64_t getri_ws = 0;

  /** resize the internal storage
   * @param norb number of electrons/orbitals
   * @param delay, maximum delay 0<delay<=norb
   */
  inline void resize(int norb, sycl::queue& m_queue)
  {
    if (Mat1_gpu.rows() != norb)
    {
      Mat1_gpu.resize(norb, norb);
      ipiv.resize(norb);
      getrf_ws = syclSolver::getrf_scratchpad_size<T_FP>(m_queue, norb, norb, norb);
      getri_ws = syclSolver::getri_scratchpad_size<T_FP>(m_queue, norb, norb);
      workspace.resize(std::max(getrf_ws, getri_ws));
    }
  }

public:
  /** compute the inverse of the transpose of matrix A and its determinant value in log
   * when T_FP and TMAT are the same
   * @tparam TREAL real type
   */
  template<typename TMAT, typename TREAL, typename = std::enable_if_t<std::is_same<TMAT, T_FP>::value>>
  std::enable_if_t<std::is_same<TMAT, T_FP>::value> invert_transpose(const Matrix<TMAT>& logdetT,
                                                                     Matrix<TMAT>& Ainv,
                                                                     Matrix<TMAT, SYCLAllocator<TMAT>>& Ainv_gpu,
                                                                     std::complex<TREAL>& log_value,
                                                                     sycl::queue& m_queue)
  {
    const int norb = logdetT.rows();
    resize(norb, m_queue);

    auto c_event = m_queue.memcpy(Mat1_gpu.data(), logdetT.data(), logdetT.size() * sizeof(TMAT));
    auto t_event = syclBLAS::transpose(m_queue, Mat1_gpu.data(), norb, Mat1_gpu.cols(), Ainv_gpu.data(), norb,
                                       Ainv_gpu.cols(), {c_event});
    try
    {
      syclSolver::getrf(m_queue, norb, norb, Ainv_gpu.data(), norb, ipiv.data(), workspace.data(), getrf_ws, {t_event})
          .wait();
    }
    catch (sycl::exception const& ex)
    {
      std::ostringstream err;
      err << "\t\tCaught synchronous SYCL exception during getrf:\n"
          << ex.what() << "  status: " << ex.code() << std::endl;
      std::cerr << err.str();
      throw std::runtime_error(err.str());
    }

    log_value = computeLogDet_sycl<TREAL>(m_queue, norb, Ainv_gpu.cols(), Ainv_gpu.data(), ipiv.data());

    c_event = syclSolver::getri(m_queue, norb, Ainv_gpu.data(), norb, ipiv.data(), workspace.data(), getri_ws);

    m_queue.memcpy(Ainv.data(), Ainv_gpu.data(), Ainv.size() * sizeof(TMAT), {c_event}).wait();
  }

  /** compute the inverse of the transpose of matrix A and its determinant value in log
   * when T_FP and TMAT are not the same
   * @tparam TREAL real type
   */
  template<typename TMAT, typename TREAL, typename = std::enable_if_t<!std::is_same<TMAT, T_FP>::value>>
  std::enable_if_t<!std::is_same<TMAT, T_FP>::value> invert_transpose(const Matrix<TMAT>& logdetT,
                                                                      Matrix<TMAT>& Ainv,
                                                                      Matrix<TMAT, SYCLAllocator<TMAT>>& Ainv_gpu,
                                                                      std::complex<TREAL>& log_value,
                                                                      sycl::queue& m_queue)
  {
    const int norb = logdetT.rows();
    resize(norb, m_queue);
    //use Ainv_gpu for transpose
    auto c_event = m_queue.memcpy(Ainv_gpu.data(), logdetT.data(), logdetT.size() * sizeof(TMAT));
    //transpose
    auto t_event = syclBLAS::transpose(m_queue, Ainv_gpu.data(), norb, Ainv_gpu.cols(), Mat1_gpu.data(), norb,
                                       Mat1_gpu.cols(), {c_event});

    //getrf (LU) -> getri (inverse)
    try
    {
      syclSolver::getrf(m_queue, norb, norb, Mat1_gpu.data(), norb, ipiv.data(), workspace.data(), getrf_ws, {t_event})
          .wait();
    }
    catch (sycl::exception const& ex)
    {
      std::ostringstream err;
      err << "\t\tCaught synchronous SYCL exception during getrf:\n"
          << ex.what() << "  status: " << ex.code() << std::endl;
      std::cerr << err.str();
      throw std::runtime_error(err.str());
    }

    log_value = computeLogDet_sycl<TREAL>(m_queue, norb, Mat1_gpu.cols(), Mat1_gpu.data(), ipiv.data());

    c_event = syclSolver::getri(m_queue, norb, Mat1_gpu.data(), norb, ipiv.data(), workspace.data(), getri_ws);

    t_event = syclBLAS::copy_n(m_queue, Mat1_gpu.data(), Mat1_gpu.size(), Ainv_gpu.data(), {c_event});

    m_queue.memcpy(Ainv.data(), Ainv_gpu.data(), Ainv.size() * sizeof(TMAT), {t_event}).wait();
  }
};
} // namespace qmcplusplus

#endif // QMCPLUSPLUS_CUSOLVERINVERTOR_H
