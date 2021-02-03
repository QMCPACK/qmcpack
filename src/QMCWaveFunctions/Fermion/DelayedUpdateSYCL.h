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


#ifndef QMCPLUSPLUS_DELAYED_UPDATE_SYCL_H
#define QMCPLUSPLUS_DELAYED_UPDATE_SYCL_H

#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "SYCL/SYCLallocator.hpp"
#include "QMCWaveFunctions/detail/SYCL/delayed_update_helper.h"
#include "QMCWaveFunctions/Fermion/DiracMatrix.h"
#include "oneapi/mkl/lapack.hpp"
#include "oneapi/mkl/blas.hpp"
namespace qmcplusplus
{
/// helper class for the prefetched range of a vector
class Range
{
  // [first, last) rows of Ainv
  int first, last;

public:
  Range() : first(0), last(0){};
  void setRange(int first_in, int last_in)
  {
    first = first_in;
    last  = last_in;
  }
  inline void clear() { first = last = 0; };
  inline int getOffset(int index) const
  {
    if (!checkRange(index))
      throw std::runtime_error("index not in range \n");
    return index - first;
  }
  inline bool checkRange(int index) const { return (index >= first) && (index < last); };
};

/** implements delayed update on NVIDIA GPU using cuBLAS and cusolverDN
 * @tparam T base precision for most computation
 * @tparam T_FP high precision for matrix inversion, T_FP >= T
 */
template<typename T, typename T_FP>
class DelayedUpdateSYCL
{
  /// define real type
  using real_type = typename scalar_traits<T>::real_type;
  // Data staged during for delayed acceptRows
  Matrix<T, SYCLManagedAllocator<T>> U;
  Matrix<T, SYCLManagedAllocator<T>> Binv;
  Matrix<T> V;
  //Matrix<T> tempMat; // for debugging only
  Matrix<T, SYCLManagedAllocator<T>> temp_gpu;
  /// GPU copy of U, V, Binv, Ainv
  Matrix<T, SYCLManagedAllocator<T>> U_gpu;
  Matrix<T, SYCLManagedAllocator<T>> V_gpu;
  Matrix<T, SYCLManagedAllocator<T>> Binv_gpu;
  Matrix<T, SYCLManagedAllocator<T>> Ainv_gpu;
  // auxiliary arrays for B
  Vector<T> p;
  Vector<int, SYCLManagedAllocator<int>> delay_list;
  Vector<int, SYCLManagedAllocator<int>> delay_list_gpu;
  /// current number of delays, increase one for each acceptance, reset to 0 after updating Ainv
  int delay_count;
  /// scratch memory for cusolverDN
  Matrix<T_FP, SYCLManagedAllocator<T_FP>> Mat1_gpu;
  /** scratch memory for cusolverDN and shared storage with Ainv_gpu.
   * The full precision build has Mat2_gpu = Ainv_gpu
   * The mixed precision build has the first half of Mat2_gpu containing Ainv_gpu
   */
  Matrix<T_FP, SYCLManagedAllocator<T_FP>> Mat2_gpu;
  /// pivot array + info
  // TODO:ipiv_gpu was int, needed to change to long for oneapi mkl call
  //      maybe need to check any transfer between ipiv/ipiv_gpu to make sure it's still ok?
  Vector<int, SYCLManagedAllocator<int>> ipiv;
  Vector<long, SYCLManagedAllocator<long>> ipiv_gpu;
  /// workspace
  Vector<T_FP, SYCLManagedAllocator<T_FP>> work_gpu;
  /// diagonal terms of LU matrix
  Vector<T_FP, SYCLManagedAllocator<T_FP>> LU_diag;
  Vector<T_FP, SYCLManagedAllocator<T_FP>> LU_diag_gpu;

  // the range of prefetched_Ainv_rows
  Range prefetched_range;
  // Ainv prefetch buffer
  Matrix<T, SYCLManagedAllocator<T>> Ainv_buffer;

  // SYCL specific variables
  //cublasHandle_t h_cublas;
  //cusolverDnHandle_t h_cusolver;
  //cudaStream_t hstream;
  //SYCLManagedAllocator<T> myalloc;
  sycl::queue q = SYCLManagedAllocator<T_FP>::q;

  inline void waitStream() { q.wait(); }

public:
  /// default constructor
  DelayedUpdateSYCL() : delay_count(0) {}

  ~DelayedUpdateSYCL() {}

  /** resize the internal storage
   * @param norb number of electrons/orbitals
   * @param delay, maximum delay 0<delay<=norb
   */
  inline void resize(int norb, int delay)
  {
    //tempMat.resize(norb, delay);
    V.resize(delay, norb);
    U.resize(delay, norb);
    p.resize(delay);
    Binv.resize(delay, delay);
    // prefetch 8% more rows corresponding to roughly 96% acceptance ratio
    Ainv_buffer.resize(std::min(static_cast<int>(delay * 1.08), norb), norb);

    temp_gpu.resize(norb, delay);
    delay_list.resize(delay);
    U_gpu.resize(delay, norb);
    V_gpu.resize(delay, norb);
    Binv_gpu.resize(delay, delay);
    delay_list_gpu.resize(delay);
    Mat1_gpu.resize(norb, norb);
    Mat2_gpu.resize(norb, norb);
    LU_diag.resize(norb);
    LU_diag_gpu.resize(norb);
    // Ainv_gpu reuses all or part of the memory of Mat2_gpu
    Ainv_gpu.attachReference(reinterpret_cast<T*>(Mat2_gpu.data()), norb, norb);
    // prepare cusolver auxiliary arrays
    ipiv.resize(norb + 1);
    ipiv_gpu.resize(norb + 1);
    const int lwork = oneapi::mkl::lapack::getrf_scratchpad_size<T_FP>(q, norb, norb, norb);
    work_gpu.resize(lwork);
  }

  /** compute the inverse of the transpose of matrix A and its determinant value in log
   * when T_FP and T are the same
   * @tparam TREAL real type
   */
  template<typename TMAT, typename TREAL, typename = std::enable_if_t<std::is_same<TMAT, T>::value>>
  std::enable_if_t<std::is_same<TMAT, T_FP>::value> invert_transpose(const Matrix<TMAT>& logdetT,
                                                                     Matrix<TMAT>& Ainv,
                                                                     std::complex<TREAL>& LogValue)
  {
    // safe mechanism
    delay_count = 0;
    int norb    = logdetT.rows();
    q.memcpy(Mat1_gpu.data(), logdetT.data(), logdetT.size() * sizeof(T));
    q.wait();

    oneapi::mkl::lapack::getrf(q, norb, norb, Mat1_gpu.data(), norb, ipiv_gpu.data(), work_gpu.data(), work_gpu.size());
    q.wait();

    q.memcpy(ipiv.data(), ipiv_gpu.data(), ipiv_gpu.size() * sizeof(int));
    q.wait();

    extract_matrix_diagonal_sycl(norb, Mat1_gpu.data(), norb, LU_diag_gpu.data(), q);
    q.wait();

    q.memcpy(LU_diag.data(), LU_diag_gpu.data(), LU_diag.size() * sizeof(T_FP));
    q.wait();

    if (ipiv[0] != 0)
    {
      // TODO: this was copied from CUDA implementation, oneapi doesn't return info in the same way
      std::ostringstream err;
      err << "oneapi::mkl::lapack::getrf calculation failed with devInfo = " << ipiv[0] << std::endl;
      std::cerr << err.str();
      throw std::runtime_error(err.str());
    }
    else
    {
      std::ostringstream err;
      err << "oneapi::mkl::lapack::getrf calculation succeeded" << std::endl;
      std::cerr << err.str();
    }
    make_identity_matrix_sycl(norb, Mat2_gpu.data(), q);

    //Todo check if CUBLAS_OP_T == oneapi::mkl::transpose::trans
    //TODO: this is still incorrect (see mixed version below)
    oneapi::mkl::lapack::getrs(q, oneapi::mkl::transpose::trans, norb, norb, Mat1_gpu.data(), norb, ipiv_gpu.data(),
                               work_gpu.data(), work_gpu.size());
    q.wait();

    q.memcpy(ipiv.data(), ipiv_gpu.data(), sizeof(int));
    q.wait();

    computeLogDet(LU_diag.data(), norb, ipiv.data() + 1, LogValue);

    q.memcpy(Ainv.data(), Ainv_gpu.data(), Ainv.size() * sizeof(T));
    // no need to wait because : For transfers from device memory to pageable host memory, the function will return only once the copy has completed.
    //waitStream();
    q.wait();
    if (ipiv[0] != 0)
    {
      // TODO: fix this
      std::ostringstream err;
      err << "oneapi::mkl::lapack::getrs calculation failed with devInfo = " << ipiv[0] << std::endl;
      std::cerr << err.str();
      throw std::runtime_error(err.str());
    }
    else
    {
      std::ostringstream err;
      err << "oneapi::mkl::lapack::getrs calculation succeeded" << std::endl;
      std::cerr << err.str();
    }
  }

  /** compute the inverse of the transpose of matrix A and its determinant value in log
   * when T_FP and T are the same
   * @tparam TREAL real type
   */
  template<typename TMAT, typename TREAL, typename = std::enable_if_t<std::is_same<TMAT, T>::value>>
  std::enable_if_t<!std::is_same<TMAT, T_FP>::value> invert_transpose(const Matrix<TMAT>& logdetT,
                                                                      Matrix<TMAT>& Ainv,
                                                                      std::complex<TREAL>& LogValue)
  {
    // safe mechanism
    delay_count = 0;
    int norb    = logdetT.rows();

    q.memcpy(Mat1_gpu.data(), logdetT.data(), logdetT.size() * sizeof(T));
    q.wait();

    copy_matrix_sycl(norb, norb, (T*)Mat1_gpu.data(), norb, Mat2_gpu.data(), norb, q);
    q.wait();

    oneapi::mkl::lapack::getrf(q, norb, norb, Mat2_gpu.data(), norb, ipiv_gpu.data(), work_gpu.data(), work_gpu.size());
    q.wait();


    q.memcpy(ipiv.data(), ipiv_gpu.data(), ipiv_gpu.size() * sizeof(int));
    q.wait();
    extract_matrix_diagonal_sycl(norb, Mat2_gpu.data(), norb, LU_diag_gpu.data(), q);
    q.memcpy(LU_diag.data(), LU_diag_gpu.data(), LU_diag.size() * sizeof(T_FP));
    q.wait();
    // check LU success
    waitStream();
    if (ipiv[0] != 0)
    {
      //TODO: fix
      //std::ostringstream err;
      //err << "cusolver::getrf calculation failed with devInfo = " << ipiv[0] << std::endl;
      std::cout << "cusolver::getrf calculation failed with devInfo = " << ipiv[0] << std::endl;
      //std::cerr << err.str();
      //throw std::runtime_error(err.str());
    }
    else
    {
      //std::ostringstream err;
      std::cout << "cusolver::getrf calculation succeeded" << std::endl;
      //err << "cusolver::getrf calculation succeeded" << std::endl;
      //std::cerr << err.str();
    }
    make_identity_matrix_sycl(norb, Mat1_gpu.data(), norb, q);
    oneapi::mkl::lapack::getrs(q, oneapi::mkl::transpose::trans, norb, norb, Mat2_gpu.data(), norb, ipiv_gpu.data(),
                               Mat1_gpu.data(), norb, work_gpu.data(), work_gpu.size());
    q.wait();

    copy_matrix_sycl(norb, norb, Mat1_gpu.data(), norb, (T*)Mat2_gpu.data(), norb, q);
    q.memcpy(ipiv.data(), ipiv_gpu.data(), sizeof(int));
    computeLogDet(LU_diag.data(), norb, ipiv.data() + 1, LogValue);

    q.memcpy(Ainv.data(), Ainv_gpu.data(), Ainv.size() * sizeof(T));
    q.wait();

    if (ipiv[0] != 0)
    {
      //TODO: fix
      //std::ostringstream err;
      std::cout << "cusolver::getrs calculation failed with devInfo = " << ipiv[0] << std::endl;
      //err << "cusolver::getrs calculation failed with devInfo = " << ipiv[0] << std::endl;
      //std::cerr << err.str();
      //throw std::runtime_error(err.str());
    }
  }

  /** initialize internal objects when Ainv is refreshed
   * @param Ainv inverse matrix
   */
  inline void initializeInv(const Matrix<T>& Ainv)
  {
    q.memcpy(Ainv_gpu.data(), Ainv.data(), Ainv.size() * sizeof(T));
    q.wait();

    delay_count = 0;
    prefetched_range.clear();
  }

  /** compute the row of up-to-date Ainv
   * @param Ainv inverse matrix
   * @param rowchanged the row id corresponding to the proposed electron
   */
  template<typename VVT>
  inline void getInvRow(const Matrix<T>& Ainv, int rowchanged, VVT& invRow)
  {
    if (!prefetched_range.checkRange(rowchanged))
    {
      int last_row = std::min(rowchanged + Ainv_buffer.rows(), Ainv.rows());

      q.memcpy(Ainv_buffer.data(), Ainv_gpu[rowchanged], invRow.size() * (last_row - rowchanged) * sizeof(T));
      prefetched_range.setRange(rowchanged, last_row);
      waitStream();
    }
    // save AinvRow to new_AinvRow
    std::copy_n(Ainv_buffer[prefetched_range.getOffset(rowchanged)], invRow.size(), invRow.data());
    if (delay_count > 0)
    {
      const T cone(1);
      const T czero(0);
      const int norb     = Ainv.rows();
      const int lda_Binv = Binv.cols();
      // multiply V (NxK) Binv(KxK) U(KxN) AinvRow right to the left
      BLAS::gemv('T', norb, delay_count, cone, U.data(), norb, invRow.data(), 1, czero, p.data(), 1);
      BLAS::gemv('N', delay_count, delay_count, cone, Binv.data(), lda_Binv, p.data(), 1, czero, Binv[delay_count], 1);
      BLAS::gemv('N', norb, delay_count, -cone, V.data(), norb, Binv[delay_count], 1, cone, invRow.data(), 1);
    }
  }

  /** accept a move with the update delayed
   * @param Ainv inverse matrix
   * @param rowchanged the row id corresponding to the proposed electron
   * @param psiV new orbital values
   *
   * Before delay_count reaches the maximum delay, only Binv is updated with a recursive algorithm
   */
  template<typename VVT>
  inline void acceptRow(Matrix<T>& Ainv, int rowchanged, const VVT& psiV)
  {
    // update Binv from delay_count to delay_count+1
    const T cminusone(-1);
    const T czero(0);
    const int norb     = Ainv.rows();
    const int lda_Binv = Binv.cols();
    std::copy_n(Ainv_buffer[prefetched_range.getOffset(rowchanged)], norb, V[delay_count]);
    std::copy_n(psiV.data(), norb, U[delay_count]);
    delay_list[delay_count] = rowchanged;
    // the new Binv is [[X Y] [Z x]]
    BLAS::gemv('T', norb, delay_count + 1, cminusone, V.data(), norb, psiV.data(), 1, czero, p.data(), 1);
    // x
    T y = -p[delay_count];
    for (int i = 0; i < delay_count; i++)
      y += Binv[delay_count][i] * p[i];
    Binv[delay_count][delay_count] = y = T(1) / y;
    // Y
    BLAS::gemv('T', delay_count, delay_count, y, Binv.data(), lda_Binv, p.data(), 1, czero, Binv.data() + delay_count,
               lda_Binv);
    // X
    BLAS::ger(delay_count, delay_count, cminusone, Binv[delay_count], 1, Binv.data() + delay_count, lda_Binv,
              Binv.data(), lda_Binv);
    // Z
    for (int i = 0; i < delay_count; i++)
      Binv[delay_count][i] *= -y;
    delay_count++;
    // update Ainv when maximal delay is reached
    if (delay_count == lda_Binv)
      updateInvMat(Ainv, false);
  }

  /** update the full Ainv and reset delay_count
   * @param Ainv inverse matrix
   */
  inline void updateInvMat(Matrix<T>& Ainv, bool transfer_to_host = true)
  {
    // update the inverse matrix
    if (delay_count > 0)
    {
      const T cone(1);
      const T czero(0);
      const int norb     = Ainv.rows();
      const int lda_Binv = Binv.cols();
      const T cminusone(-1);

      q.memcpy(U_gpu.data(), U.data(), norb * delay_count * sizeof(T));
      q.wait();

      oneapi::mkl::blas::gemm(q, oneapi::mkl::transpose::trans, oneapi::mkl::transpose::nontrans, delay_count, norb,
                              norb, cone, U_gpu.data(), norb, Ainv_gpu.data(), norb, czero, temp_gpu.data(), lda_Binv);
      q.wait();

      q.memcpy(delay_list_gpu.data(), delay_list.data(), delay_count * sizeof(int));
      q.wait();
      applyW_stageV_sycl(delay_list_gpu.data(), delay_count, temp_gpu.data(), norb, temp_gpu.cols(), V_gpu.data(),
                         Ainv_gpu.data(), q);
      q.wait();
      q.memcpy(Binv_gpu.data(), Binv.data(), lda_Binv * delay_count * sizeof(T));
      q.wait();
      oneapi::mkl::blas::gemm(q, oneapi::mkl::transpose::trans, oneapi::mkl::transpose::nontrans, norb, delay_count,
                              delay_count, cone, V_gpu.data(), norb, Binv_gpu.data(), lda_Binv, czero, U_gpu.data(),
                              norb);
      q.wait();

      oneapi::mkl::blas::gemm(q, oneapi::mkl::transpose::trans, oneapi::mkl::transpose::nontrans, norb, norb,
                              delay_count, cminusone, U_gpu.data(), norb, temp_gpu.data(), lda_Binv, cone,
                              Ainv_gpu.data(), norb);
      q.wait();

      delay_count = 0;
      // Ainv is invalid, reset range
      prefetched_range.clear();
    }

    // transfer Ainv_gpu to Ainv and wait till completion
    if (transfer_to_host)
    {
      q.memcpy(Ainv.data(), Ainv_gpu.data(), Ainv.size() * sizeof(T));
      q.wait();
    }
  }
};
} // namespace qmcplusplus

#endif // QMCPLUSPLUS_DELAYED_UPDATE_SYCL_H
