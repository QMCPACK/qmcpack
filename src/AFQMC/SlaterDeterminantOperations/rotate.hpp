////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////


#ifndef AFQMC_ROTATE_HPP
#define AFQMC_ROTATE_HPP

#include <numeric>
#include "AFQMC/config.h"
#include "Utilities/FairDivide.h"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Matrix/csr_matrix.hpp"
#include "AFQMC/Numerics/ma_operations.hpp"
#include "AFQMC/Numerics/csr_blas.hpp"
#include "AFQMC/Matrix/matrix_emplace_wrapper.hpp"
#include "multi/array_ref.hpp"
#include "AFQMC/Utilities/afqmc_TTI.hpp"
#include "mpi.h"

namespace qmcplusplus
{
namespace afqmc
{
namespace sparse_rotate
{
/*
 *  Performs a (left) half rotation (and a possible transposition) of a Cholesky matrix.
 *  The rotated matrix is sotored in "tuple" form in the provided container using emplace_back.
 *  The generation of the matrix is distributed among the containers
 *  of all the cores in the node.
 *  Input:
 *    -alpha: Sparse hermitian of rotation matrix for spin up, e.g. transpose(conj(Aup)).
 *    -beta: Sparse hermitian of rotation matrix for spin down, e.g. transpose(conj(Aup)).
 *    -A: Input sparse cholesky matrix.
 *    -cutoff: Value below which elements of rotated matrix are ignored.
 *  Output:
 *    -B: Container (e.g. std::vector) with a segment of the non-zero terms of the rotated matrix.
 *        If transpose==false, the elements in the container are ordered by (a,k) values.
 *
 *  If transposed==true:
 *     B(n,ka-ka0) = sum_i^M alpha(a,i) * Spvn(ik,n)
 *     B(n,ka+N*M-ka0) = sum_i^M beta(a,i) * Spvn(ki,n)
 *  else:
 *     B(ka-ka0,n) = sum_i A(a,i) * Spvn(ik,n)
 *     B(ka+N*M-ka0,n) = sum_i^M beta(a,i) * Spvn(ik,n),
 *  where M/N is the number of rows/columns of alpha and beta.
 *  The number of rows of Spvn should be equal to M*M.
 *  If conjV=true, Spvn(ik,n) -> conj(Spvn(ki,n)) in the expressions above
 */
template<class Container, class task_group, class PsiT_Type>
void halfRotateCholeskyMatrix(WALKER_TYPES type,
                              task_group& TG,
                              int k0,
                              int kN,
                              Container& Q,
                              PsiT_Type* Alpha,
                              PsiT_Type* Beta,
                              SpVType_shm_csr_matrix const& CholMat,
                              bool transpose,
                              bool conjV           = false,
                              double cutoff        = 1e-6,
                              bool reserve_to_fit_ = true)
{
  int NAEA = Alpha->size(0);
  int NAEB = Alpha->size(0);
  int NMO  = Alpha->size(1);
  if (type == COLLINEAR)
    NAEB = Beta->size(0);
  int nvec   = std::get<1>(CholMat.sizes());
  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();

  assert(std::get<0>(CholMat.sizes()) == NMO * NMO);
  assert(kN > k0);
  if (type == CLOSED && kN > NMO)
    APP_ABORT(" Error: kN > NMO in halfRotateCholeskyMatrix. \n");

  // map from [0:2*NMO) to [0:NMO) in collinear case
  int k0_alpha = 0, k0_beta = 0, kN_alpha = 0, kN_beta = 0;
  if (k0 >= 0)
  {
    k0_alpha = std::min(k0, NMO);
    kN_alpha = std::min(kN, NMO);
    if (type == COLLINEAR)
    {
      kN_beta = std::max(kN, NMO) - NMO;
      k0_beta = std::max(k0, NMO) - NMO;
    }
  }

  int ak0, ak1;
  int Qdim = NAEA * (kN_alpha - k0_alpha) + NAEB * (kN_beta - k0_beta);
  if (transpose)
  {
    //if(not check_shape(Q,{nvec,Qdim}))
    if (not(Q.size(0) == nvec && Q.size(1) == Qdim))
      APP_ABORT(" Error: Container Q has incorrect dimensions in halfRotateCholeskyMatrix. \n");
  }
  else
  {
    //if(not check_shape(Q,{Qdim,nvec}))
    if (not(Q.size(0) == Qdim && Q.size(1) == nvec))
      APP_ABORT(" Error: Container Q has incorrect dimensions in halfRotateCholeskyMatrix. \n");
  }
  std::tie(ak0, ak1) = FairDivideBoundary(coreid, Qdim, ncores);

  if (type == NONCOLLINEAR)
    APP_ABORT(" GHF not yet implemented. \n");

  boost::multi::array<SPComplexType, 1> vec(iextensions<1u>{nvec});
  if (reserve_to_fit_)
  {
    std::vector<std::size_t> sz_per_row(Qdim);
    int cnt = 0;
    for (int k = k0_alpha; k < kN_alpha; k++)
    {
      for (int a = 0; a < NAEA; a++, cnt++)
      {
        if (cnt < ak0)
          continue;
        if (cnt >= ak1)
          break;
        std::fill_n(vec.origin(), vec.size(), SPComplexType(0, 0));
        auto Aa = (*Alpha)[a];
        for (int ip = 0; ip < Aa.num_non_zero_elements(); ++ip)
        {
          auto Aai = Aa.non_zero_values_data()[ip];
          auto i   = Aa.non_zero_indices2_data()[ip];
          if (conjV)
            csr::axpy('C', Aai, CholMat[k * NMO + i], vec);
          else
            csr::axpy('N', Aai, CholMat[i * NMO + k], vec);
        }
        if (transpose)
        {
          for (int n = 0; n < nvec; n++)
            if (std::abs(vec[n]) > cutoff)
              ++sz_per_row[n];
        }
        else
        {
          for (int n = 0; n < nvec; n++)
            if (std::abs(vec[n]) > cutoff)
              ++sz_per_row[cnt];
        }
      }
    }
    // reset "amount of work done" to full alpha piece
    cnt = NAEA * (kN_alpha - k0_alpha);
    if (type == COLLINEAR)
    {
      // reset "shift"
      for (int k = k0_beta; k < kN_beta; k++)
      {
        for (int a = 0; a < NAEB; a++, cnt++)
        {
          if (cnt < ak0)
            continue;
          if (cnt >= ak1)
            break;
          std::fill_n(vec.origin(), vec.size(), SPComplexType(0, 0));
          auto Aa = (*Beta)[a];
          for (int ip = 0; ip < Aa.num_non_zero_elements(); ++ip)
          {
            auto Aai = Aa.non_zero_values_data()[ip];
            auto i   = Aa.non_zero_indices2_data()[ip];
            if (conjV)
              csr::axpy('C', Aai, CholMat[k * NMO + i], vec);
            else
              csr::axpy('N', Aai, CholMat[i * NMO + k], vec);
          }
          if (transpose)
          {
            for (int n = 0; n < nvec; n++)
              if (std::abs(vec[n]) > cutoff)
                ++sz_per_row[n];
          }
          else
          {
            for (int n = 0; n < nvec; n++)
              if (std::abs(vec[n]) > cutoff)
                ++sz_per_row[cnt];
          }
        }
      }
    }
    TG.Node().all_reduce_in_place_n(sz_per_row.begin(), sz_per_row.size(), std::plus<>());
    reserve_to_fit(Q, sz_per_row);
  } // else
  //  Q.reserve( std::size_t(0.1 * nvec * NEL * NMO / ncores ) ); // by default, assume 10% sparsity

  int cnt = 0;
  for (int k = k0_alpha; k < kN_alpha; k++)
  {
    for (int a = 0; a < NAEA; a++, cnt++)
    {
      if (cnt < ak0)
        continue;
      if (cnt >= ak1)
        break;
      std::fill_n(vec.origin(), vec.size(), SPComplexType(0, 0));
      auto Aa = (*Alpha)[a];
      for (int ip = 0; ip < Aa.num_non_zero_elements(); ++ip)
      {
        auto Aai = Aa.non_zero_values_data()[ip];
        auto i   = Aa.non_zero_indices2_data()[ip];
        if (conjV)
          csr::axpy('C', Aai, CholMat[k * NMO + i], vec);
        else
          csr::axpy('N', Aai, CholMat[i * NMO + k], vec);
      }
      if (transpose)
      {
        for (int n = 0; n < nvec; n++)
          if (std::abs(vec[n]) > cutoff)
            emplace(Q, std::forward_as_tuple(n, cnt, vec[n]));
      }
      else
      {
        for (int n = 0; n < nvec; n++)
          if (std::abs(vec[n]) > cutoff)
            emplace(Q, std::forward_as_tuple(cnt, n, vec[n]));
      }
    }
  }
  // reset "amount of work done" to full alpha piece
  cnt = NAEA * (kN_alpha - k0_alpha);
  if (type == COLLINEAR)
  {
    // reset "shift"
    for (int k = k0_beta; k < kN_beta; k++)
    {
      for (int a = 0; a < NAEB; a++, cnt++)
      {
        if (cnt < ak0)
          continue;
        if (cnt >= ak1)
          break;
        std::fill_n(vec.origin(), vec.size(), SPComplexType(0, 0));
        auto Aa = (*Beta)[a];
        for (int ip = 0; ip < Aa.num_non_zero_elements(); ++ip)
        {
          auto Aai = Aa.non_zero_values_data()[ip];
          auto i   = Aa.non_zero_indices2_data()[ip];
          if (conjV)
            csr::axpy('C', Aai, CholMat[k * NMO + i], vec);
          else
            csr::axpy('N', Aai, CholMat[i * NMO + k], vec);
        }
        if (transpose)
        {
          for (int n = 0; n < nvec; n++)
            if (std::abs(vec[n]) > cutoff)
              emplace(Q, std::forward_as_tuple(n, cnt, vec[n]));
        }
        else
        {
          for (int n = 0; n < nvec; n++)
            if (std::abs(vec[n]) > cutoff)
              emplace(Q, std::forward_as_tuple(cnt, n, vec[n]));
        }
      }
    }
  }
}

/*
 * Calculates the rotated Cholesky matrix used in the calculation of the vias potential.
 *     v(n,ak) = sum_i A(a,i) * Spvn(ik,n)     A(a,i)=PsiT(i,a)*
 *     v(n,N*M+ak) = sum_i^M beta(a,i) * Spvn(ik,n),
 *  where M/N is the number of rows/columns of alpha and beta.
 *  The number of rows of Spvn should be equal to M*M.
 *  Since we only rotate the "local" cholesky matrix, the algorithm is only parallelized
 *  over a node. In principle, it is possible to spread this over equivalent nodes.
 */
template<typename T,
         class task_group,
         class PsiT_Type,
         typename = typename std::enable_if_t<std::is_same<T, std::complex<double>>::value or
                                              std::is_same<T, std::complex<float>>::value>>
SpCType_shm_csr_matrix halfRotateCholeskyMatrixForBias(WALKER_TYPES type,
                                                       task_group& TG,
                                                       PsiT_Type* Alpha,
                                                       PsiT_Type* Beta,
                                                       SpVType_shm_csr_matrix const& CholMat,
                                                       double cutoff = 1e-6)
{
  int NAEA = Alpha->size(0);
  int NAEB = Alpha->size(0);
  int NMO  = Alpha->size(1);
  if (type != CLOSED)
    NAEB = Beta->size(0);
  int nvec   = std::get<1>(CholMat.sizes());
  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();

  // to speed up, generate new communicator for eqv_nodes and split full work among all
  // cores in this comm. Then build from distributed container?

  assert(std::get<0>(CholMat.sizes()) == NMO * NMO);

  std::size_t Qdim = NAEA * NMO;
  if (type == COLLINEAR)
    Qdim += NAEB * NMO;
  if (type == NONCOLLINEAR)
    Qdim = 2 * NMO * (NAEA + NAEB);
  std::size_t ak0, ak1;
  std::tie(ak0, ak1) = FairDivideBoundary(std::size_t(coreid), Qdim, std::size_t(ncores));

  if (type == NONCOLLINEAR)
    APP_ABORT(" GHF not yet implemented. \n");

  boost::multi::array<SPComplexType, 1> vec(iextensions<1u>{nvec});
  std::vector<std::size_t> sz_per_row(nvec);
  std::size_t cnt = 0;
  for (int a = 0; a < NAEA; a++)
  {
    for (int k = 0; k < NMO; k++, cnt++)
    {
      if (cnt < ak0)
        continue;
      if (cnt >= ak1)
        break;
      std::fill_n(vec.origin(), vec.size(), SPComplexType(0, 0));
      auto Aa = (*Alpha)[a];
      for (int ip = 0; ip < Aa.num_non_zero_elements(); ++ip)
      {
        auto Aai = Aa.non_zero_values_data()[ip];
        auto i   = Aa.non_zero_indices2_data()[ip];
        csr::axpy('N', Aai, CholMat[i * NMO + k], vec);
      }
      for (std::size_t n = 0; n < nvec; n++)
        if (std::abs(vec[n]) > cutoff)
          ++(sz_per_row[n]);
    }
  }

  // reset "amount of work done" to full alpha piece
  cnt = NAEA * NMO;
  if (type == COLLINEAR)
  {
    // reset "shift"
    for (int a = 0; a < NAEB; a++)
    {
      for (int k = 0; k < NMO; k++, cnt++)
      {
        if (cnt < ak0)
          continue;
        if (cnt >= ak1)
          break;
        std::fill_n(vec.origin(), vec.size(), SPComplexType(0, 0));
        auto Aa = (*Beta)[a];
        for (int ip = 0; ip < Aa.num_non_zero_elements(); ++ip)
        {
          auto Aai = Aa.non_zero_values_data()[ip];
          auto i   = Aa.non_zero_indices2_data()[ip];
          csr::axpy('N', Aai, CholMat[i * NMO + k], vec);
        }
        for (std::size_t n = 0; n < nvec; n++)
          if (std::abs(vec[n]) > cutoff)
            ++sz_per_row[n];
      }
    }
  }
  TG.Node().all_reduce_in_place_n(sz_per_row.begin(), sz_per_row.size(), std::plus<>());

  using Alloc = shared_allocator<SPComplexType>;
  SpCType_shm_csr_matrix::base ucsr(tp_ul_ul{nvec, Qdim}, tp_ul_ul{0, 0}, sz_per_row, Alloc(TG.Node()));

  using mat_wrapper = csr::matrix_emplace_wrapper<SpCType_shm_csr_matrix::base>;
  mat_wrapper ucsr_wrapper(ucsr, TG.Node());

  cnt = 0;
  for (int a = 0; a < NAEA; a++)
  {
    for (int k = 0; k < NMO; k++, cnt++)
    {
      if (cnt < ak0)
        continue;
      if (cnt >= ak1)
        break;
      std::fill_n(vec.origin(), vec.size(), SPComplexType(0, 0));
      auto Aa = (*Alpha)[a];
      for (int ip = 0; ip < Aa.num_non_zero_elements(); ++ip)
      {
        auto Aai = Aa.non_zero_values_data()[ip];
        auto i   = Aa.non_zero_indices2_data()[ip];
        csr::axpy('N', Aai, CholMat[i * NMO + k], vec);
      }
      for (std::size_t n = 0; n < nvec; n++)
        if (std::abs(vec[n]) > cutoff)
          ucsr_wrapper.emplace(std::forward_as_tuple(n, cnt, vec[n]));
    }
  }
  // reset "amount of work done" to full alpha piece
  cnt = NAEA * NMO;
  if (type == COLLINEAR)
  {
    // reset "shift"
    for (int a = 0; a < NAEB; a++)
    {
      for (int k = 0; k < NMO; k++, cnt++)
      {
        if (cnt < ak0)
          continue;
        if (cnt >= ak1)
          break;
        std::fill_n(vec.origin(), vec.size(), SPComplexType(0, 0));
        auto Aa = (*Beta)[a];
        for (int ip = 0; ip < Aa.num_non_zero_elements(); ++ip)
        {
          auto Aai = Aa.non_zero_values_data()[ip];
          auto i   = Aa.non_zero_indices2_data()[ip];
          csr::axpy('N', Aai, CholMat[i * NMO + k], vec);
        }
        for (std::size_t n = 0; n < nvec; n++)
          if (std::abs(vec[n]) > cutoff)
            ucsr_wrapper.emplace(std::forward_as_tuple(n, cnt, vec[n]));
      }
    }
  }
  ucsr_wrapper.push_buffer();
  return SpCType_shm_csr_matrix(std::move(ucsr));
}

// This is needed to allow generic implementation of SparseTensorIO while catching errors at compile time.
template<typename T,
         class task_group,
         class PsiT_Type,
         typename = typename std::enable_if_t<not(std::is_same<T, std::complex<double>>::value or
                                                  std::is_same<T, std::complex<float>>::value)>,
         typename = void>
SpVType_shm_csr_matrix halfRotateCholeskyMatrixForBias(WALKER_TYPES type,
                                                       task_group& TG,
                                                       PsiT_Type* Alpha,
                                                       PsiT_Type* Beta,
                                                       SpVType_shm_csr_matrix const& CholMat,
                                                       double cutoff = 1e-6)
{
  print_stacktrace throw std::runtime_error(
      "Error: Incorrect template parameter in halfRotateCholeskyMatrixForBias. \n");
  using Alloc = shared_allocator<SPRealType>;
  SpVType_shm_csr_matrix csr(tp_ul_ul{1, 1}, tp_ul_ul{0, 0}, 1, Alloc(TG.Node()));
  return csr;
}

} // namespace sparse_rotate

namespace ma_rotate
{
template<class MultiArray2D, class task_group, class PsiT_Type>
void halfRotateCholeskyMatrix(WALKER_TYPES type,
                              task_group& TG,
                              int k0,
                              int kN,
                              MultiArray2D&& Q,
                              PsiT_Type* Alpha,
                              PsiT_Type* Beta,
                              SpVType_shm_csr_matrix const& CholMat,
                              bool transpose,
                              bool conjV    = false,
                              double cutoff = 1e-6)
{
  int NAEA = Alpha->size(0);
  int NAEB = 0;
  int NMO  = Alpha->size(1);
  if (type == COLLINEAR)
    NAEB = Beta->size(0);
  int nvec   = std::get<1>(CholMat.sizes());
  int ncores = TG.getTotalCores(), coreid = TG.getCoreID();

  assert(std::get<0>(CholMat.sizes()) == NMO * NMO);
  if (type == CLOSED && kN > NMO)
    APP_ABORT(" Error: kN > NMO in halfRotateCholeskyMatrix. \n");

  // map from [0:2*NMO) to [0:NMO) in collinear case
  int k0_alpha = 0, k0_beta = 0, kN_alpha = 0, kN_beta = 0;
  if (k0 >= 0)
  {
    k0_alpha = std::min(k0, NMO);
    kN_alpha = std::min(kN, NMO);
    if (type == COLLINEAR)
    {
      k0_beta = std::max(k0, NMO) - NMO;
      kN_beta = std::max(kN, NMO) - NMO;
    }
  }

  int ak0, ak1;
  int Qdim = NAEA * (kN_alpha - k0_alpha) + NAEB * (kN_beta - k0_beta);
  if (transpose)
  {
    assert(std::get<0>(Q.sizes()) == nvec);
    assert(std::get<1>(Q.sizes()) == Qdim);
  }
  else
  {
    assert(std::get<0>(Q.sizes()) == Qdim);
    assert(std::get<1>(Q.sizes()) == nvec);
  }
  std::tie(ak0, ak1) = FairDivideBoundary(coreid, Qdim, ncores);

  if (type == NONCOLLINEAR)
    APP_ABORT(" GHF not yet implemented. \n");

  int cnt = 0;
  for (int k = k0_alpha; k < kN_alpha; k++)
  {
    for (int a = 0; a < NAEA; a++, cnt++)
    {
      if (cnt < ak0)
        continue;
      if (cnt >= ak1)
        break;
      if (transpose)
      {
        auto vec = Q(Q.extension(0), cnt);
        for (auto& v : vec)
          v = SPComplexType(0, 0);
        auto Aa = (*Alpha)[a];
        for (int ip = 0; ip < Aa.num_non_zero_elements(); ++ip)
        {
          auto Aai = Aa.non_zero_values_data()[ip];
          auto i   = Aa.non_zero_indices2_data()[ip];
          if (conjV)
            csr::axpy('C', Aai, CholMat[k * NMO + i], vec);
          else
            csr::axpy('N', Aai, CholMat[i * NMO + k], vec);
        }
      }
      else
      {
        auto vec = Q[cnt];
        for (auto& v : vec)
          v = SPComplexType(0, 0);
        auto Aa = (*Alpha)[a];
        for (int ip = 0; ip < Aa.num_non_zero_elements(); ++ip)
        {
          auto Aai = Aa.non_zero_values_data()[ip];
          auto i   = Aa.non_zero_indices2_data()[ip];
          if (conjV)
            csr::axpy('C', Aai, CholMat[k * NMO + i], vec);
          else
            csr::axpy('N', Aai, CholMat[i * NMO + k], vec);
        }
      }
    }
  }
  // reset "amount of work done" to full alpha piece
  cnt = NAEA * (kN_alpha - k0_alpha);
  if (type == COLLINEAR)
  {
    // reset "shift"
    for (int k = k0_beta; k < kN_beta; k++)
    {
      for (int a = 0; a < NAEB; a++, cnt++)
      {
        if (cnt < ak0)
          continue;
        if (cnt >= ak1)
          break;
        if (transpose)
        {
          auto vec = Q(Q.extension(0), cnt);
          for (auto& v : vec)
            v = SPComplexType(0, 0);
          auto Aa = (*Beta)[a];
          for (int ip = 0; ip < Aa.num_non_zero_elements(); ++ip)
          {
            auto Aai = Aa.non_zero_values_data()[ip];
            auto i   = Aa.non_zero_indices2_data()[ip];
            if (conjV)
              csr::axpy('C', Aai, CholMat[k * NMO + i], vec);
            else
              csr::axpy('N', Aai, CholMat[i * NMO + k], vec);
          }
        }
        else
        {
          auto vec = Q[cnt];
          for (auto& v : vec)
            v = SPComplexType(0, 0);
          auto Aa = (*Beta)[a];
          for (int ip = 0; ip < Aa.num_non_zero_elements(); ++ip)
          {
            auto Aai = Aa.non_zero_values_data()[ip];
            auto i   = Aa.non_zero_indices2_data()[ip];
            if (conjV)
              csr::axpy('C', Aai, CholMat[k * NMO + i], vec);
            else
              csr::axpy('N', Aai, CholMat[i * NMO + k], vec);
          }
        }
      }
    }
  }
  TG.Node().barrier();
}

// design for compact arrays
/*
 * Constructs the following contraction of the Slater matrix of the trial wfn and the Cholesky Matrix (Likn).
 * Case:
 *   - Closed/Collinear:  L[a][n][k] = sum_i A[a][i] L[i][k][n]
 *       - In collinear case, two separate calls are made for each spin channel.
 *   - Non-collinear: L[a][n][sk] = sum_i A[a][si] L[i][k][n]   // [si] == [s][i] combined spinor index
 *       - In this case, to preserve matrix dimensions, [s][k] --> [sk] is kept as a single index.
 */
template<class MultiArray2DA, class MultiArray3DB, class MultiArray3DC, class MultiArray2D>
void getLank(MultiArray2DA&& Aai,
             MultiArray3DB&& Likn,
             MultiArray3DC&& Lank,
             MultiArray2D&& buff,
             bool noncollinear = false)
{
  int npol = noncollinear ? 2 : 1;
  int na   = std::get<0>(Aai.sizes());
  if (na == 0)
    return;
  int ni    = std::get<1>(Aai.sizes()) / npol;
  int nk    = std::get<1>(Likn.sizes());
  int nchol = std::get<2>(Likn.sizes());
  assert(std::get<0>(Likn.sizes()) == ni);
  assert(std::get<0>(Lank.sizes()) == na);
  assert(std::get<1>(Lank.sizes()) == nchol);
  assert(std::get<2>(Lank.sizes()) == nk * npol);
  assert(std::get<0>(buff.sizes()) >= npol * nk);
  assert(std::get<1>(buff.sizes()) >= nchol);
  if (noncollinear)
    assert(Aai.stride(0) == std::get<1>(Aai.sizes())); // make sure it is contiguous

  using elementA = typename std::decay<MultiArray2DA>::type::element;
  using element  = typename std::decay<MultiArray3DC>::type::element;
  boost::multi::array_ref<elementA, 2> Aas_i(to_address(Aai.origin()), {na * npol, ni});
  boost::multi::array_ref<element, 2> Li_kn(to_address(Likn.origin()), {ni, nk * nchol});
  boost::multi::array_ref<element, 2> Las_kn(to_address(Lank.origin()), {na * npol, nk * nchol});

  ma::product(Aas_i, Li_kn, Las_kn);
  for (int a = 0; a < na; a++)
  {
    boost::multi::array_ref<element, 2> Lskn(to_address(Lank[a].origin()), {npol * nk, nchol});
    boost::multi::array_ref<element, 2> Lnsk(to_address(Lank[a].origin()), {nchol, npol * nk});
    buff({0, npol * nk}, {0, nchol}) = Lskn;
    ma::transpose(buff({0, npol * nk}, {0, nchol}), Lnsk);
  }
}

/*
 * Constructs the following contraction of the Slater matrix of the trial wfn and the Cholesky Matrix (Likn).
 * Case:
 *   - Closed/Collinear:  L[a][n][k] = sum_i A[a][i] conj(L[k][i][n])
 *       - In collinear case, two separate calls are made for each spin channel.
 *   - Non-collinear: L[a][n][sk] = sum_i A[a][si] conj(L[k][i][n])   // [si] == [s][i] combined spinor index
 *       - In this case, to preserve matrix dimensions, [s][k] --> [sk] is kept as a single index.
 */
template<class MultiArray2DA, class MultiArray3DB, class MultiArray3DC, class MultiArray2D>
void getLank_from_Lkin(MultiArray2DA&& Aai,
                       MultiArray3DB&& Lkin,
                       MultiArray3DC&& Lank,
                       MultiArray2D&& buff,
                       bool noncollinear = false)
{
  int npol = noncollinear ? 2 : 1;
  int na   = std::get<0>(Aai.sizes());
  if (na == 0)
    return;
  int ni    = std::get<1>(Aai.sizes()) / npol;
  int nk    = std::get<0>(Lkin.sizes());
  int nchol = std::get<2>(Lkin.sizes());
  assert(std::get<1>(Lkin.sizes()) == ni);
  assert(std::get<0>(Lank.sizes()) == na);
  assert(std::get<1>(Lank.sizes()) == nchol);
  assert(std::get<2>(Lank.sizes()) == nk * npol);
  assert(buff.num_elements() >= na * npol * nchol);
  if (noncollinear)
    assert(Aai.stride(0) == std::get<1>(Aai.sizes())); // make sure it is contiguous

  using Type     = typename std::decay<MultiArray3DC>::type::element;
  using elementA = typename std::decay<MultiArray2DA>::type::element;
  boost::multi::array_ref<elementA, 2> Aas_i(to_address(Aai.origin()), {na * npol, ni});
  boost::multi::array_ref<Type, 2> bnas(to_address(buff.origin()), {nchol, na * npol});
  // Lank[a][n][k] = sum_i Aai[a][i] conj(Lkin[k][i][n])
  // Lank[as][n][k] = sum_i Aai[as][i] conj(Lkin[k][i][n])
  for (int k = 0; k < nk; k++)
  {
    ma::product(ma::H(Lkin[k]), ma::T(Aas_i), bnas);
    for (int a = 0; a < na; a++)
      for (int n = 0; n < nchol; n++)
        for (int p = 0; p < npol; p++)
          Lank[a][n][p * nk + k] = bnas[n][a * npol + p];
  }
}

} // namespace ma_rotate

namespace ma_rotate_padded
{
// designed for padded arrays
template<class MultiArray2DA, class MultiArray3DB, class MultiArray3DC>
void getLakn_Lank(MultiArray2DA&& Aai,
                  MultiArray3DB&& Likn,
                  MultiArray3DC&& Lakn,
                  MultiArray3DC&& Lank,
                  bool noncollinear = false)
{
  int npol = noncollinear ? 2 : 1;
  int na   = std::get<0>(Aai.sizes());
  if (na == 0)
    return;
  int ni = std::get<1>(Aai.sizes()) / npol;

  int nmo   = std::get<0>(Likn.sizes());
  int nchol = std::get<2>(Likn.sizes());
  assert(std::get<1>(Likn.sizes()) == nmo);

  assert(std::get<1>(Lakn.sizes()) == npol * nmo);
  assert(std::get<2>(Lakn.sizes()) == nchol);

  assert(std::get<0>(Lakn.sizes()) >= na);
  assert(std::get<0>(Lakn.sizes()) == std::get<0>(Lank.sizes()));
  assert(std::get<1>(Lank.sizes()) == nchol);
  assert(std::get<2>(Lank.sizes()) == npol * nmo);

  if (noncollinear)
    assert(Aai.stride(0) == std::get<1>(Aai.sizes())); // make sure it is contiguous

  using elmA = typename std::decay<MultiArray2DA>::type::element;
  using elmB = typename std::decay<MultiArray3DB>::type::element;
  using elmC = typename std::decay<MultiArray3DC>::type::element;

  boost::multi::array_ref<elmA, 2> Aas_i(to_address(Aai.origin()), {na * npol, ni});
  boost::multi::array_ref<elmB, 2, decltype(Likn.origin())> Li_kn(Likn.origin(), {ni, nmo * nchol});
  boost::multi::array_ref<elmC, 2, decltype(Lakn.origin())> Las_kn(Lakn.origin(), {na * npol, nmo * nchol});

  ma::product(Aas_i, Li_kn, Las_kn);
  for (int a = 0; a < na; a++)
    ma::transpose(Lakn[a], Lank[a]);
}

template<class MultiArray2DA, class MultiArray3DB, class MultiArray3DC, class MultiArray2D>
void getLakn_Lank_from_Lkin(MultiArray2DA&& Aai,
                            MultiArray3DB&& Lkin,
                            MultiArray3DC&& Lakn,
                            MultiArray3DC&& Lank,
                            MultiArray2D&& buff,
                            bool noncollinear = false)
{
  int npol = noncollinear ? 2 : 1;
  int na   = std::get<0>(Aai.sizes());
  if (na == 0)
    return;
  int ni = std::get<1>(Aai.sizes()) / npol;

  int nmo   = std::get<0>(Lkin.sizes());
  int nchol = std::get<2>(Lkin.sizes());
  assert(std::get<1>(Lkin.sizes()) == nmo);

  assert(std::get<1>(Lakn.sizes()) == npol * nmo);
  assert(std::get<2>(Lakn.sizes()) == nchol);

  assert(std::get<0>(Lakn.sizes()) >= na);
  assert(std::get<0>(Lakn.sizes()) == std::get<0>(Lank.sizes()));
  assert(std::get<1>(Lank.sizes()) == nchol);
  assert(std::get<2>(Lank.sizes()) == npol * nmo);

  if (noncollinear)
    assert(Aai.stride(0) == std::get<1>(Aai.sizes())); // make sure it is contiguous

  assert(buff.num_elements() >= na * npol * nchol);

  using ptr2 = typename std::decay<MultiArray2D>::type::element_ptr;
  using elm2 = typename std::decay<MultiArray2D>::type::element;
  using elmA = typename std::decay<MultiArray2DA>::type::element;

  boost::multi::array_ref<elmA, 2> Aas_i(to_address(Aai.origin()), {na * npol, ni});
  boost::multi::array_ref<elm2, 2, ptr2> bnas(buff.origin(), {nchol, na * npol});
  // Lakn[a][sk][n] = sum_i Aai[as][i] conj(Lkin[k][i][n])
  for (int k = 0; k < nmo; k++)
  {
    ma::product(ma::H(Lkin[k].sliced(0, ni)), ma::T(Aas_i), bnas);
    for (int a = 0; a < na; a++)
      for (int p = 0; p < npol; p++)
        Lakn[a][p * nmo + k] = bnas({0, nchol}, a * npol + p);
  }
  for (int a = 0; a < na; a++)
    ma::transpose(Lakn[a], Lank[a]);
}


} // namespace ma_rotate_padded


} // namespace afqmc

} // namespace qmcplusplus

#endif
