////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
// Alfredo Correa, correaa@llnl.gov
//    Lawrence Livermore National Laboratory
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////


#ifndef AFQMC_APPLY_EXPM_HPP
#define AFQMC_APPLY_EXPM_HPP

#include "AFQMC/Numerics/ma_operations.hpp"
#include "Utilities/FairDivide.h"

namespace qmcplusplus
{
namespace afqmc
{
namespace SlaterDeterminantOperations
{
namespace base
{
/*
 * Calculate S = exp(im*V)*S using a Taylor expansion of exp(V)
 */
template<class MatA, class MatB, class MatC>
inline void apply_expM(const MatA& V, MatB&& S, MatC& T1, MatC& T2, int order = 6, char TA = 'N')
{
  assert(std::get<0>(V.sizes()) == std::get<1>(V.sizes()));
  assert(std::get<1>(V.sizes()) == std::get<0>(S.sizes()));
  assert(std::get<0>(S.sizes()) == std::get<0>(T1.sizes()));
  assert(std::get<1>(S.sizes()) == std::get<1>(T1.sizes()));
  assert(std::get<0>(S.sizes()) == std::get<0>(T2.sizes()));
  assert(std::get<1>(S.sizes()) == std::get<1>(T2.sizes()));

  using ma::H;
  using ma::T;
  using ComplexType = typename std::decay<MatB>::type::element;
  ComplexType zero(0.);
  auto pT1(std::addressof(T1));
  auto pT2(std::addressof(T2));

  ComplexType im(0.0, 1.0);
  if (TA == 'H' || TA == 'h')
    im = ComplexType(0.0, -1.0);

  // getting around issue in multi, fix later
  //T1 = S;
  T1.sliced(0, std::get<0>(T1.sizes())) = S;
  for (int n = 1; n <= order; n++)
  {
    ComplexType fact = im * static_cast<ComplexType>(1.0 / static_cast<double>(n));
    if (TA == 'H' || TA == 'h')
      ma::product(fact, ma::H(V), *pT1, zero, *pT2);
    else if (TA == 'T' || TA == 't')
      ma::product(fact, ma::T(V), *pT1, zero, *pT2);
    else
      ma::product(fact, V, *pT1, zero, *pT2);
    ma::add(ComplexType(1.0), *pT2, ComplexType(1.0), S, S);
    std::swap(pT1, pT2);
  }
}

} // namespace base

namespace shm
{
/*
 * Calculate S = exp(im*V)*S using a Taylor expansion of exp(V)
 * V, S, T1, T2 are expected to be in shared memory.  
 */
template<class MatA, class MatB, class MatC, class communicator>
inline void apply_expM(const MatA& V, MatB&& S, MatC& T1, MatC& T2, communicator& comm, int order = 6, char TA = 'N')
{
  assert(std::get<0>(V.sizes()) == std::get<0>(S.sizes()));
  assert(std::get<1>(V.sizes()) == std::get<0>(S.sizes()));
  assert(std::get<0>(S.sizes()) == std::get<0>(T1.sizes()));
  assert(std::get<1>(S.sizes()) == std::get<1>(T1.sizes()));
  assert(std::get<0>(S.sizes()) == std::get<0>(T2.sizes()));
  assert(std::get<1>(S.sizes()) == std::get<1>(T2.sizes()));

  using ComplexType = typename std::decay<MatB>::type::element;

  const ComplexType zero(0.);
  ComplexType im(0.0, 1.0);
  if (TA == 'H' || TA == 'h')
    im = ComplexType(0.0, -1.0);

  auto pT1(std::addressof(T1));
  auto pT2(std::addressof(T2));

  int M0, Mn;
  std::tie(M0, Mn) = FairDivideBoundary(comm.rank(), int(S.size()), comm.size());

  assert(M0 <= Mn);
  assert(M0 >= 0);

  T1.sliced(M0, Mn) = S.sliced(M0, Mn);
  comm.barrier();
  for (int n = 1; n <= order; n++)
  {
    const ComplexType fact = im * static_cast<ComplexType>(1.0 / static_cast<double>(n));
    if (TA == 'H' || TA == 'h')
      ma::product(fact, ma::H(V(V.extension(0), {M0, Mn})), *pT1, zero, (*pT2).sliced(M0, Mn));
    else if (TA == 'T' || TA == 't')
      ma::product(fact, ma::T(V(V.extension(0), {M0, Mn})), *pT1, zero, (*pT2).sliced(M0, Mn));
    else
      ma::product(fact, V.sliced(M0, Mn), *pT1, zero, (*pT2).sliced(M0, Mn));
    // overload += ???
    for (int i = M0; i < Mn; i++)
      for (int j = 0, je = std::get<1>(S.sizes()); j < je; j++)
        S[i][j] += (*pT2)[i][j];
    comm.barrier();
    std::swap(pT1, pT2);
  }
}


} // namespace shm

namespace batched
{
/*
 * Calculate S = exp(im*V)*S using a Taylor expansion of exp(V)
 */
template<class MatA, class MatB, class MatC>
inline void apply_expM(const MatA& V, MatB&& S, MatC& T1, MatC& T2, int order = 6, char TA = 'N')
{
  static_assert(std::decay<MatA>::type::dimensionality == 3, " batched::apply_expM::dimenionality == 3");
  static_assert(std::decay<MatB>::type::dimensionality == 3, " batched::apply_expM::dimenionality == 3");
  static_assert(std::decay<MatC>::type::dimensionality == 3, " batched::apply_expM::dimenionality == 3");
  assert(std::get<0>(V.sizes()) == std::get<0>(S.sizes()));
  assert(std::get<0>(V.sizes()) == std::get<0>(T1.sizes()));
  assert(std::get<0>(V.sizes()) == std::get<0>(T2.sizes()));
  assert(std::get<1>(V.sizes()) == std::get<2>(V.sizes()));
  assert(std::get<2>(V.sizes()) == std::get<1>(S.sizes()));
  assert(std::get<1>(S.sizes()) == std::get<1>(T1.sizes()));
  assert(std::get<2>(S.sizes()) == std::get<2>(T1.sizes()));
  assert(std::get<1>(S.sizes()) == std::get<1>(T2.sizes()));
  assert(std::get<2>(S.sizes()) == std::get<2>(T2.sizes()));
  // for now limit to continuous
  assert(S.stride(0) == std::get<1>(S.sizes()) * std::get<2>(S.sizes()));
  assert(T1.stride(0) == std::get<1>(T1.sizes()) * std::get<2>(T1.sizes()));
  assert(T2.stride(0) == std::get<1>(T2.sizes()) * std::get<2>(T2.sizes()));
  assert(S.stride(1) == std::get<2>(S.sizes()));
  assert(T1.stride(1) == std::get<2>(T1.sizes()));
  assert(T2.stride(1) == std::get<2>(T2.sizes()));
  assert(S.stride(2) == 1);
  assert(T1.stride(2) == 1);
  assert(T2.stride(2) == 1);

  using ComplexType = typename std::decay<MatB>::type::element;
  ComplexType zero(0.);
  ComplexType im(0.0, 1.0);
  if (TA == 'H' || TA == 'h')
    im = ComplexType(0.0, -1.0);
  auto pT1(std::addressof(T1));
  auto pT2(std::addressof(T2));

  // getting around issue in multi, fix later
  //T1 = S;
  using std::copy_n;
  copy_n(S.origin(), S.num_elements(), T1.origin());
  for (int n = 1; n <= order; n++)
  {
    ComplexType fact = im * static_cast<ComplexType>(1.0 / static_cast<double>(n));
    if (TA == 'H' || TA == 'h')
      ma::productStridedBatched(fact, ma::H(V), *pT1, zero, *pT2);
    else if (TA == 'T' || TA == 't')
      ma::productStridedBatched(fact, ma::T(V), *pT1, zero, *pT2);
    else
      ma::productStridedBatched(fact, V, *pT1, zero, *pT2);
    //ma::add(ComplexType(1.0),*pT2,ComplexType(1.0),S,S);
    using ma::axpy;
    axpy(S.num_elements(), ComplexType(1.0), (*pT2).origin(), 1, S.origin(), 1);
    std::swap(pT1, pT2);
  }
}

/*
 * Calculate S = exp(im*V)*S using a Taylor expansion of exp(V)
 * Version for non_collinear calculations, where there are 2 S matrices per V in the batch.
 */
template<class MatA, class MatB, class MatC>
inline void apply_expM_noncollinear(const MatA& V, MatB&& S, MatC& T1, MatC& T2, int order = 6, char TA = 'N')
{
  static_assert(std::decay<MatA>::type::dimensionality == 3, " batched::apply_expM::dimenionality == 3");
  static_assert(std::decay<MatB>::type::dimensionality == 3, " batched::apply_expM::dimenionality == 3");
  static_assert(std::decay<MatC>::type::dimensionality == 3, " batched::apply_expM::dimenionality == 3");
  assert(std::get<0>(V.sizes()) * 2 == std::get<0>(S.sizes()));
  assert(std::get<0>(V.sizes()) * 2 == std::get<0>(T1.sizes()));
  assert(std::get<0>(V.sizes()) * 2 == std::get<0>(T2.sizes()));
  assert(std::get<1>(V.sizes()) == std::get<2>(V.sizes()));
  assert(std::get<2>(V.sizes()) == std::get<1>(S.sizes()));
  assert(std::get<1>(S.sizes()) == std::get<1>(T1.sizes()));
  assert(std::get<2>(S.sizes()) == std::get<2>(T1.sizes()));
  assert(std::get<1>(S.sizes()) == std::get<1>(T2.sizes()));
  assert(std::get<2>(S.sizes()) == std::get<2>(T2.sizes()));
  // for now limit to continuous
  assert(S.stride(0) == std::get<1>(S.sizes()) * std::get<2>(S.sizes()));
  assert(T1.stride(0) == std::get<1>(T1.sizes()) * std::get<2>(T1.sizes()));
  assert(T2.stride(0) == std::get<1>(T2.sizes()) * std::get<2>(T2.sizes()));
  assert(S.stride(1) == std::get<2>(S.sizes()));
  assert(T1.stride(1) == std::get<2>(T1.sizes()));
  assert(T2.stride(1) == std::get<2>(T2.sizes()));
  assert(S.stride(2) == 1);
  assert(T1.stride(2) == 1);
  assert(T2.stride(2) == 1);

  using ComplexType = typename std::decay<MatB>::type::element;
  ComplexType zero(0.);
  ComplexType im(0.0, 1.0);
  if (TA == 'H' || TA == 'h')
    im = ComplexType(0.0, -1.0);
  auto pT1(std::addressof(T1));
  auto pT2(std::addressof(T2));

  using pointerA = typename std::decay<MatA>::type::element_const_ptr;
  using pointerC = typename std::decay<MatC>::type::element_ptr;

  int nbatch = S.size();
  int ldv    = V.stride(1);
  int M      = std::get<2>(T2.sizes());
  int N      = std::get<1>(T2.sizes());
  int K      = std::get<1>(T1.sizes());

  std::vector<pointerA> Vi;
  std::vector<pointerC> T1i;
  std::vector<pointerC> T2i;
  Vi.reserve(2 * V.size());
  T1i.reserve(T1.size());
  T2i.reserve(T2.size());
  for (int i = 0; i < V.size(); i++)
  {
    Vi.emplace_back(ma::pointer_dispatch(V[i].origin()));
    Vi.emplace_back(ma::pointer_dispatch(V[i].origin()));
  }
  for (int i = 0; i < T1.size(); i++)
    T1i.emplace_back(ma::pointer_dispatch(T1[i].origin()));
  for (int i = 0; i < T2.size(); i++)
    T2i.emplace_back(ma::pointer_dispatch(T2[i].origin()));

  auto pT1i(std::addressof(T1i));
  auto pT2i(std::addressof(T2i));

  // getting around issue in multi, fix later
  //T1 = S;
  using std::copy_n;
  copy_n(S.origin(), S.num_elements(), T1.origin());
  for (int n = 1; n <= order; n++)
  {
    ComplexType fact = im * static_cast<ComplexType>(1.0 / static_cast<double>(n));
    using ma::gemmBatched;
    // careful with fortran ordering
    gemmBatched('N', TA, M, N, K, fact, pT1i->data(), (*pT1).stride(1), Vi.data(), ldv, zero, pT2i->data(),
                (*pT2).stride(1), nbatch);
    using ma::axpy;
    axpy(S.num_elements(), ComplexType(1.0), (*pT2).origin(), 1, S.origin(), 1);
    std::swap(pT1, pT2);
    std::swap(pT1i, pT2i);
  }
}

} // namespace batched

} // namespace SlaterDeterminantOperations

} // namespace afqmc

} // namespace qmcplusplus

#endif
