//////////////////////////////////////////////////////////////////////////////
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

#ifndef QMCPLUSPLUS_AFQMC_CREATEHAMILTONIAN_HELPER_H
#define QMCPLUSPLUS_AFQMC_CREATEHAMILTONIAN_HELPER_H

#include <cstdlib>
#include <algorithm>
#include <complex>
#include <iostream>
#include <vector>
#include <numeric>
#if defined(USE_MPI)
#include <mpi.h>
#endif

#include "Configuration.h"
#include "AFQMC/config.h"

#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Hamiltonians/Hamiltonian_Utilities.hpp"

namespace qmcplusplus
{
namespace afqmc
{
template<typename T>
inline bool isComplex(const T& a)
{
  return std::is_same<T, std::complex<RealType>>::value;
}

template<typename T>
inline ComplexType toComplex(const T& a)
{
  return a;
}

template<>
inline ComplexType toComplex(const RealType& a)
{
  return ComplexType(a, 0.0);
}

inline void count_allowed_terms(int NMO,
                                std::vector<std::size_t>& sz,
                                std::vector<s4D<ValueType>>& vs4D,
                                std::map<IndexType, std::pair<bool, IndexType>>& occ_a)
{
  for (auto& v : vs4D)
    if (occ_a[std::get<0>(v)].first && occ_a[std::get<1>(v)].first)
      ++sz[occ_a[std::get<0>(v)].second * NMO + std::get<2>(v)];
}

inline void count_allowed_terms(std::vector<std::size_t>& sz,
                                std::vector<s4D<ValueType>>& vs4D,
                                std::map<IndexType, std::pair<bool, IndexType>>& occ_a,
                                std::map<IndexType, std::pair<bool, IndexType>>& occ_b)
{}

template<class Container>
inline void add_allowed_terms(int NMO,
                              std::vector<s4D<ValueType>>& vs4D,
                              std::map<IndexType, std::pair<bool, IndexType>>& occ_a,
                              Container& V)
{
  for (auto& vi : vs4D)
    if (occ_a[std::get<0>(vi)].first && occ_a[std::get<1>(vi)].first)
      V.emplace(occ_a[std::get<0>(vi)].second * NMO + std::get<2>(vi),
                occ_a[std::get<1>(vi)].second * NMO + std::get<3>(vi), static_cast<SPValueType>(std::get<4>(vi)));
  //V.emplace( std::get<0>(vi)*NMO + std::get<2>(vi),
  //           std::get<1>(vi)*NMO + std::get<3>(vi),
  //           static_cast<SPValueType>(std::get<4>(vi)));
}

template<class Container>
inline void add_allowed_terms(int NMO,
                              std::vector<s4D<ValueType>>& vs4D,
                              std::map<IndexType, std::pair<bool, IndexType>>& occ_a,
                              std::map<IndexType, std::pair<bool, IndexType>>& occ_b,
                              Container& V,
                              bool needs_locks = false,
                              bool GHF         = false)
{}

inline void push_ijkl(int NMO,
                      OrbitalType i,
                      OrbitalType j,
                      OrbitalType k,
                      OrbitalType l,
                      ValueType V,
                      std::vector<s4D<ValueType>>& v,
                      bool GHF = false)
{
  long ik, jl;
  if (GHF)
  {
    ik = i * 2 * NMO + k;
    jl = j * 2 * NMO + l;
  }
  else
  {
    ik = (i < NMO) ? (i * NMO + k) : (NMO * NMO + (i - NMO) * NMO + k - NMO);
    jl = (j < NMO) ? (j * NMO + l) : (NMO * NMO + (j - NMO) * NMO + l - NMO);
  }
  if (ik == jl)
    v.push_back(std::make_tuple(i, j, k, l, V));
  else if (ik < jl)
    v.push_back(std::make_tuple(i, j, k, l, ValueType(2.0) * V));
  else
    v.push_back(std::make_tuple(j, i, l, k, ValueType(2.0) * V));
}

// For a given quartet ijkl and the 3 (for real) different permutations
// in the reduced list, add all possible contributions to the Hamiltonian
// taking into account cutoff and (ik)<->(jl) symmetry
// J1 = <ij|kl>
// J2 = <ij|lk>
// J3 = <ik|jl>
//  routine currently assumes that ijkl is the smallest of the 3 non-symmetric terms
//  For complex, there are 3 extra non-symmetric terms:
//  J1a = <il|kj>
//  J2a = <ik|lj>
//  J3a = <il|jk>
//    In this case, make sure you use the fact that: <ij|kl> = ma::conj( <kl|ij> )
inline void find_all_contributions_to_hamiltonian_closed_shell(int NMO,
                                                               OrbitalType i,
                                                               OrbitalType j,
                                                               OrbitalType k,
                                                               OrbitalType l,
                                                               ValueType J1,
                                                               ValueType J2,
                                                               ValueType J3,
                                                               ValueType J1a,
                                                               ValueType J2a,
                                                               ValueType J3a,
                                                               double cut,
                                                               std::vector<s4D<ValueType>>& v)
{
  // simple algorithm for now
  // 1. add all contributions blindly
  // 2. check for repeated and remove
  // 3. apply symmetry elimination

  using ma::conj;

  v.reserve(24);

#ifndef QMC_COMPLEX
  J1a = J1;
  J2a = J2;
  J3a = J3;
#endif

  // <ij||kl> -> (ijkl)
  // Symmetries:
  //   (ijkl) = (jilk) = (klij)* = (lkji)* = -(ijlk) = -(jikl) = -(lkij)* = -(klji)*

  ValueType J1J2 = 4.0 * J1 - 2.0 * J2;
  if (std::abs(J1J2) > cut)
  {
    push_ijkl(NMO, i, j, k, l, J1J2, v);
    push_ijkl(NMO, k, l, i, j, ma::conj(J1J2), v);
  }
  ValueType J1J3 = 4.0 * J1a - 2.0 * J3a;
  if (std::abs(J1J3) > cut)
  {
    push_ijkl(NMO, i, l, k, j, J1J3, v);
    push_ijkl(NMO, j, k, l, i, ma::conj(J1J3), v);
  }

  ValueType J2J1 = 4.0 * J2 - 2.0 * J1;
  if (std::abs(J2J1) > cut)
  {
    push_ijkl(NMO, i, j, l, k, J2J1, v);
    push_ijkl(NMO, k, l, j, i, ma::conj(J2J1), v);
  }
  ValueType J2J3 = 4.0 * J2a - 2.0 * J3;
  if (std::abs(J2J3) > cut)
  {
    push_ijkl(NMO, i, k, l, j, J2J3, v);
    push_ijkl(NMO, j, l, k, i, ma::conj(J2J3), v);
  }

  ValueType J3J1 = 4.0 * J3a - 2.0 * J1a;
  if (std::abs(J3J1) > cut)
  {
    push_ijkl(NMO, i, l, j, k, J3J1, v);
    push_ijkl(NMO, k, j, l, i, ma::conj(J3J1), v);
  }
  ValueType J3J2 = 4.0 * J3 - 2.0 * J2a;
  if (std::abs(J3J2) > cut)
  {
    push_ijkl(NMO, i, k, j, l, J3J2, v);
    push_ijkl(NMO, l, j, k, i, ma::conj(J3J2), v);
  }

  // order to remove consecutive repreated
  std::sort(v.begin(), v.end(), [](const s4D<ValueType>& lhs, const s4D<ValueType>& rhs) {
    return std::forward_as_tuple(std::get<0>(lhs), std::get<1>(lhs), std::get<2>(lhs), std::get<3>(lhs)) <
        std::forward_as_tuple(std::get<0>(rhs), std::get<1>(rhs), std::get<2>(rhs), std::get<3>(rhs));
  });
  // remove consecutive repeated
  std::vector<s4D<ValueType>>::iterator it;
  it = std::unique(v.begin(), v.end());
  // resize array, since unique does not remove elements
  // just reorders them
  v.resize(std::distance(v.begin(), it));
}

// For a given quartet ijkl and the 3 (for real) different permutations
// in the reduced list, add all possible contributions to the Hamiltonian
// taking into account cutoff and (ik)<->(jl) symmetry
// J1 = <ij|kl>
// J2 = <ij|lk>
// J3 = <ik|jl>
//  For complex, there are 3 extra non-symmetric terms:
//  J1a = <il|kj>
//  J2a = <ik|lj>
//  J3a = <il|jk>
inline void find_all_contributions_to_hamiltonian_collinear(int NMO,
                                                            OrbitalType i,
                                                            OrbitalType j,
                                                            OrbitalType k,
                                                            OrbitalType l,
                                                            ValueType J1,
                                                            ValueType J2,
                                                            ValueType J3,
                                                            ValueType J1a,
                                                            ValueType J2a,
                                                            ValueType J3a,
                                                            double cut,
                                                            std::vector<s4D<ValueType>>& v)
{
  // simple algorithm for now
  // 1. add all contributions blindly
  // 2. check for repeated and remove
  // 3. apply symmetry elimination
  //

  v.reserve(48);
  v.clear();

  int i2 = i + NMO;
  int j2 = j + NMO;
  int k2 = k + NMO;
  int l2 = l + NMO;
#ifndef QMC_COMPLEX
  J1a = J1;
  J2a = J2;
  J3a = J3;
#else
  APP_ABORT("Not yet working for std::complex \n\n\n");
#endif
  // 2bar terms
  ValueType J1J2 = J1 - J2;
  if (std::abs(J1J2) > cut)
  {
    push_ijkl(NMO, i, j, k, l, J1J2, v);
    push_ijkl(NMO, k, l, i, j, J1J2, v);
    push_ijkl(NMO, i2, j2, k2, l2, J1J2, v);
    push_ijkl(NMO, k2, l2, i2, j2, J1J2, v);
  }
  ValueType J1J3 = J1 - J3;
  if (std::abs(J1J3) > cut)
  {
    push_ijkl(NMO, i, l, k, j, J1J3, v);
    push_ijkl(NMO, j, k, l, i, J1J3, v);
    push_ijkl(NMO, i2, l2, k2, j2, J1J3, v);
    push_ijkl(NMO, j2, k2, l2, i2, J1J3, v);
  }

  ValueType J2J1 = J2 - J1;
  if (std::abs(J2J1) > cut)
  {
    push_ijkl(NMO, i, j, l, k, J2J1, v);
    push_ijkl(NMO, k, l, j, i, J2J1, v);
    push_ijkl(NMO, i2, j2, l2, k2, J2J1, v);
    push_ijkl(NMO, k2, l2, j2, i2, J2J1, v);
  }
  ValueType J2J3 = J2 - J3;
  if (std::abs(J2J3) > cut)
  {
    push_ijkl(NMO, i, k, l, j, J2J3, v);
    push_ijkl(NMO, j, l, k, i, J2J3, v);
    push_ijkl(NMO, i2, k2, l2, j2, J2J3, v);
    push_ijkl(NMO, j2, l2, k2, i2, J2J3, v);
  }

  ValueType J3J1 = J3 - J1;
  if (std::abs(J3J1) > cut)
  {
    push_ijkl(NMO, i, l, j, k, J3J1, v);
    push_ijkl(NMO, k, j, l, i, J3J1, v);
    push_ijkl(NMO, i2, l2, j2, k2, J3J1, v);
    push_ijkl(NMO, k2, j2, l2, i2, J3J1, v);
  }
  ValueType J3J2 = J3 - J2;
  if (std::abs(J3J2) > cut)
  {
    push_ijkl(NMO, i, k, j, l, J3J2, v);
    push_ijkl(NMO, l, j, k, i, J3J2, v);
    push_ijkl(NMO, i2, k2, j2, l2, J3J2, v);
    push_ijkl(NMO, l2, j2, k2, i2, J3J2, v);
  }

  // 1bar terms
  if (std::abs(J1) > cut)
  { // && !aa_only) {
    push_ijkl(NMO, i, j2, k, l2, J1, v);
    push_ijkl(NMO, k, l2, i, j2, J1, v);
    push_ijkl(NMO, i, l2, k, j2, J1, v);
    push_ijkl(NMO, j, k2, l, i2, J1, v);
    push_ijkl(NMO, j, i2, l, k2, J1, v);
    push_ijkl(NMO, l, k2, j, i2, J1, v);
    push_ijkl(NMO, l, i2, j, k2, J1, v);
    push_ijkl(NMO, k, j2, i, l2, J1, v);
  }
  if (std::abs(J2) > cut)
  { // && !aa_only) {
    push_ijkl(NMO, i, j2, l, k2, J2, v);
    push_ijkl(NMO, k, l2, j, i2, J2, v);
    push_ijkl(NMO, i, k2, l, j2, J2, v);
    push_ijkl(NMO, j, l2, k, i2, J2, v);
    push_ijkl(NMO, j, i2, k, l2, J2, v);
    push_ijkl(NMO, l, k2, i, j2, J2, v);
    push_ijkl(NMO, k, i2, j, l2, J2, v);
    push_ijkl(NMO, l, j2, i, k2, J2, v);
  }
  if (std::abs(J3) > cut)
  { // && !aa_only) {
    push_ijkl(NMO, i, l2, j, k2, J3, v);
    push_ijkl(NMO, k, j2, l, i2, J3, v);
    push_ijkl(NMO, i, k2, j, l2, J3, v);
    push_ijkl(NMO, l, j2, k, i2, J3, v);
    push_ijkl(NMO, l, i2, k, j2, J3, v);
    push_ijkl(NMO, j, k2, i, l2, J3, v);
    push_ijkl(NMO, k, i2, l, j2, J3, v);
    push_ijkl(NMO, j, l2, i, k2, J3, v);
  }


  // order to remove consecutive repreated
  std::sort(v.begin(), v.end(), [](const s4D<ValueType>& lhs, const s4D<ValueType>& rhs) {
    return std::forward_as_tuple(std::get<0>(lhs), std::get<1>(lhs), std::get<2>(lhs), std::get<3>(lhs)) <
        std::forward_as_tuple(std::get<0>(rhs), std::get<1>(rhs), std::get<2>(rhs), std::get<3>(rhs));
  });
  // remove consecutive repeated
  std::vector<s4D<ValueType>>::iterator it;
  it = std::unique(v.begin(), v.end());
  // resize array, since unique does not remove elements
  // just reorders them
  v.resize(std::distance(v.begin(), it));
}

// For a given quartet ijkl and the 3 (for real) different permutations
// in the reduced list, add all possible contributions to the Hamiltonian
// taking into account cutoff and (ik)<->(jl) symmetry
// J1 = <ij|kl>
// J2 = <ij|lk>
// J3 = <ik|jl>
//  For complex, there are 3 extra non-symmetric terms:
//  J1a = <il|kj>
//  J2a = <ik|lj>
//  J3a = <il|jk>
inline void find_all_contributions_to_hamiltonian_ghf(int NMO,
                                                      OrbitalType i,
                                                      OrbitalType j,
                                                      OrbitalType k,
                                                      OrbitalType l,
                                                      ValueType J1,
                                                      ValueType J2,
                                                      ValueType J3,
                                                      ValueType J1a,
                                                      ValueType J2a,
                                                      ValueType J3a,
                                                      double cut,
                                                      std::vector<s4D<ValueType>>& v)
{
  // simple algorithm for now
  // 1. add all contributions blindly
  // 2. check for repeated and remove
  // 3. apply symmetry elimination
  //
  // GHF: all terms from _spinRestricted routine plus: M_(i-alpha k-beta)_(j-beta l-alpha) = <ij|lk>

  v.reserve(48);
  v.clear();

  int i2 = i + NMO;
  int j2 = j + NMO;
  int k2 = k + NMO;
  int l2 = l + NMO;
#ifndef QMC_COMPLEX
  J1a = J1;
  J2a = J2;
  J3a = J3;
#else
  APP_ABORT("Not yet working for std::complex \n\n\n");
#endif
  // 2bar terms
  ValueType J1J2 = J1 - J2;
  if (std::abs(J1J2) > cut)
  {
    push_ijkl(NMO, i, j, k, l, J1J2, v, true);
    push_ijkl(NMO, k, l, i, j, J1J2, v, true);
    push_ijkl(NMO, i2, j2, k2, l2, J1J2, v, true);
    push_ijkl(NMO, k2, l2, i2, j2, J1J2, v, true);
  }
  ValueType J1J3 = J1 - J3;
  if (std::abs(J1J3) > cut)
  {
    push_ijkl(NMO, i, l, k, j, J1J3, v, true);
    push_ijkl(NMO, j, k, l, i, J1J3, v, true);
    push_ijkl(NMO, i2, l2, k2, j2, J1J3, v, true);
    push_ijkl(NMO, j2, k2, l2, i2, J1J3, v, true);
  }

  ValueType J2J1 = J2 - J1;
  if (std::abs(J2J1) > cut)
  {
    push_ijkl(NMO, i, j, l, k, J2J1, v, true);
    push_ijkl(NMO, k, l, j, i, J2J1, v, true);
    push_ijkl(NMO, i2, j2, l2, k2, J2J1, v, true);
    push_ijkl(NMO, k2, l2, j2, i2, J2J1, v, true);
  }
  ValueType J2J3 = J2 - J3;
  if (std::abs(J2J3) > cut)
  {
    push_ijkl(NMO, i, k, l, j, J2J3, v, true);
    push_ijkl(NMO, j, l, k, i, J2J3, v, true);
    push_ijkl(NMO, i2, k2, l2, j2, J2J3, v, true);
    push_ijkl(NMO, j2, l2, k2, i2, J2J3, v, true);
  }

  ValueType J3J1 = J3 - J1;
  if (std::abs(J3J1) > cut)
  {
    push_ijkl(NMO, i, l, j, k, J3J1, v, true);
    push_ijkl(NMO, k, j, l, i, J3J1, v, true);
    push_ijkl(NMO, i2, l2, j2, k2, J3J1, v, true);
    push_ijkl(NMO, k2, j2, l2, i2, J3J1, v, true);
  }
  ValueType J3J2 = J3 - J2;
  if (std::abs(J3J2) > cut)
  {
    push_ijkl(NMO, i, k, j, l, J3J2, v, true);
    push_ijkl(NMO, l, j, k, i, J3J2, v, true);
    push_ijkl(NMO, i2, k2, j2, l2, J3J2, v, true);
    push_ijkl(NMO, l2, j2, k2, i2, J3J2, v, true);
  }

  // 1bar terms: <alpha,beta|alpha,beta>
  if (std::abs(J1) > cut)
  {
    push_ijkl(NMO, i, j2, k, l2, J1, v, true);
    push_ijkl(NMO, k, l2, i, j2, J1, v, true);
    push_ijkl(NMO, i, l2, k, j2, J1, v, true);
    push_ijkl(NMO, j, k2, l, i2, J1, v, true);
    push_ijkl(NMO, j, i2, l, k2, J1, v, true);
    push_ijkl(NMO, l, k2, j, i2, J1, v, true);
    push_ijkl(NMO, l, i2, j, k2, J1, v, true);
    push_ijkl(NMO, k, j2, i, l2, J1, v, true);
  }
  if (std::abs(J2) > cut)
  {
    push_ijkl(NMO, i, j2, l, k2, J2, v, true);
    push_ijkl(NMO, k, l2, j, i2, J2, v, true);
    push_ijkl(NMO, i, k2, l, j2, J2, v, true);
    push_ijkl(NMO, j, l2, k, i2, J2, v, true);
    push_ijkl(NMO, j, i2, k, l2, J2, v, true);
    push_ijkl(NMO, l, k2, i, j2, J2, v, true);
    push_ijkl(NMO, k, i2, j, l2, J2, v, true);
    push_ijkl(NMO, l, j2, i, k2, J2, v, true);
  }
  if (std::abs(J3) > cut)
  {
    push_ijkl(NMO, i, l2, j, k2, J3, v, true);
    push_ijkl(NMO, k, j2, l, i2, J3, v, true);
    push_ijkl(NMO, i, k2, j, l2, J3, v, true);
    push_ijkl(NMO, l, j2, k, i2, J3, v, true);
    push_ijkl(NMO, l, i2, k, j2, J3, v, true);
    push_ijkl(NMO, j, k2, i, l2, J3, v, true);
    push_ijkl(NMO, k, i2, l, j2, J3, v, true);
    push_ijkl(NMO, j, l2, i, k2, J3, v, true);
  }

  // GHF extra exchage terms: <alpha,beta|beta,alpha>
  // M_iakb_jbla = <iajb|lakb>
  // J1 = <ij|kl>
  // J2 = <ij|lk>
  // J3 = <ik|jl>
  J1 *= RealType(-1);
  J2 *= RealType(-1);
  J3 *= RealType(-1);
  if (std::abs(J1) > cut)
  {
    push_ijkl(NMO, i, j2, l2, k, J1, v, true);
    push_ijkl(NMO, k, j2, l2, i, J1, v, true);
    push_ijkl(NMO, i, l2, j2, k, J1, v, true);
    push_ijkl(NMO, k, l2, j2, i, J1, v, true);
    push_ijkl(NMO, j, i2, k2, l, J1, v, true);
    push_ijkl(NMO, l, i2, k2, j, J1, v, true);
    push_ijkl(NMO, j, k2, i2, l, J1, v, true);
    push_ijkl(NMO, l, k2, i2, j, J1, v, true);
  }
  if (std::abs(J2) > cut)
  {
    push_ijkl(NMO, i, j2, k2, l, J2, v, true);
    push_ijkl(NMO, l, j2, k2, i, J2, v, true);
    push_ijkl(NMO, i, k2, j2, l, J2, v, true);
    push_ijkl(NMO, l, k2, j2, i, J2, v, true);
    push_ijkl(NMO, j, i2, l2, k, J2, v, true);
    push_ijkl(NMO, k, i2, l2, j, J2, v, true);
    push_ijkl(NMO, j, l2, i2, k, J2, v, true);
    push_ijkl(NMO, k, l2, i2, j, J2, v, true);
  }
  if (std::abs(J3) > cut)
  {
    push_ijkl(NMO, i, k2, l2, j, J3, v, true);
    push_ijkl(NMO, j, k2, l2, i, J3, v, true);
    push_ijkl(NMO, i, l2, k2, j, J3, v, true);
    push_ijkl(NMO, j, l2, k2, i, J3, v, true);
    push_ijkl(NMO, k, i2, j2, l, J3, v, true);
    push_ijkl(NMO, l, i2, j2, k, J3, v, true);
    push_ijkl(NMO, k, j2, i2, l, J3, v, true);
    push_ijkl(NMO, l, j2, i2, k, J3, v, true);
  }


  // order to remove consecutive repreated
  std::sort(v.begin(), v.end(), [](const s4D<ValueType>& lhs, const s4D<ValueType>& rhs) {
    return std::forward_as_tuple(std::get<0>(lhs), std::get<1>(lhs), std::get<2>(lhs), std::get<3>(lhs)) <
        std::forward_as_tuple(std::get<0>(rhs), std::get<1>(rhs), std::get<2>(rhs), std::get<3>(rhs));
  });
  // remove consecutive repeated
  std::vector<s4D<ValueType>>::iterator it;
  it = std::unique(v.begin(), v.end());
  // resize array, since unique does not remove elements
  // just reorders them
  v.resize(std::distance(v.begin(), it));
}

// For a given quartet ijkl and the 3 (for real) different permutations
// in the reduced list, add all possible contributions to the Hamiltonian
// taking into account cutoff and (ik)<->(jl) symmetry
// J1 = <ij|kl>
// J2 = <ij|lk>
// J3 = <ik|jl>
//  For complex, there are 3 extra non-symmetric terms:
//  J1a = <il|kj>
//  J2a = <ik|lj>
//  J3a = <il|jk>
inline void find_all_contributions_to_hamiltonian_general(int NMO,
                                                          OrbitalType i,
                                                          OrbitalType j,
                                                          OrbitalType k,
                                                          OrbitalType l,
                                                          ValueType J1,
                                                          ValueType J2,
                                                          ValueType J3,
                                                          ValueType J1a,
                                                          ValueType J2a,
                                                          ValueType J3a,
                                                          double cut,
                                                          std::vector<s4D<ValueType>>& v)
{
  APP_ABORT("Finsigh implementation. \n\n\n");
}

// looks for all equivalent terms associated with Vijkl
// Applies std::complex ma::conjugate when needed
// Quite inefficient, so don't use outside initialization
inline void find_equivalent_OneBar_for_integral_list(s4D<ValueType> Vijkl, std::vector<s4D<ValueType>>& v)
{
  v.reserve(24);
  v.clear();
  IndexType i, j, k, l;
  ValueType V;
  std::tie(i, j, k, l, V) = Vijkl;
  if (isComplex(std::get<4>(Vijkl)))
  {
    // only alpha/alpha sector is stored
    // so there should bever be any beta index here, should I check???
    // <ij|kl> = <ji|lk> = <kl|ij>* = <lk|ji>*
    s3D<IndexType> ijkl = std::make_tuple(i, j, k, l);
    s3D<IndexType> jilk = std::make_tuple(j, i, l, k);
    s3D<IndexType> klij = std::make_tuple(k, l, i, j);
    s3D<IndexType> lkji = std::make_tuple(l, k, j, i);

    v.push_back(Vijkl);
    if (jilk != ijkl)
      v.push_back(std::make_tuple(j, i, l, k, V));
    if (klij != ijkl && klij != jilk)
      v.push_back(std::make_tuple(k, l, i, j, ma::conj(V)));
    if (lkji != ijkl && lkji != jilk && lkji != klij)
      v.push_back(std::make_tuple(l, k, j, i, ma::conj(V)));
    // order, not needed but good measure
    std::sort(v.begin(), v.end(), [](const s4D<ValueType>& lhs, const s4D<ValueType>& rhs) {
      return std::forward_as_tuple(std::get<0>(lhs), std::get<1>(lhs), std::get<2>(lhs), std::get<3>(lhs)) <
          std::forward_as_tuple(std::get<0>(rhs), std::get<1>(rhs), std::get<2>(rhs), std::get<3>(rhs));
    });
  }
  else
  {
    // only alpha/alpha sector is stored
    // so there should bever be any beta index here, should I check???
    //  For one-bar integrals: <ij|kl> = <kj|il> = <il|kj> = <kl|ij>
    //                                 = <ji|lk> = <li|jk> = <jk|li> = <lk|ji>

    v.push_back(std::make_tuple(i, j, k, l, V));
    v.push_back(std::make_tuple(k, j, i, l, V));
    v.push_back(std::make_tuple(i, l, k, j, V));
    v.push_back(std::make_tuple(k, l, i, j, V));
    v.push_back(std::make_tuple(j, i, l, k, V));
    v.push_back(std::make_tuple(l, i, j, k, V));
    v.push_back(std::make_tuple(j, k, l, i, V));
    v.push_back(std::make_tuple(l, k, j, i, V));
    // order to remove consecutive repreated


    std::sort(v.begin(), v.end(), [](const s4D<ValueType>& lhs, const s4D<ValueType>& rhs) {
      return std::forward_as_tuple(std::get<0>(lhs), std::get<1>(lhs), std::get<2>(lhs), std::get<3>(lhs)) <
          std::forward_as_tuple(std::get<0>(rhs), std::get<1>(rhs), std::get<2>(rhs), std::get<3>(rhs));
    });
    //std::sort (v.begin(), v.end());
    // remove consecutive repeated
    std::vector<s4D<ValueType>>::iterator it;
    it = std::unique(v.begin(), v.end());
    // resize array, since unique does not remove elements
    // just reorders them
    v.resize(std::distance(v.begin(), it));
  }
}

// for symmetric terms (e.g. ik/jl - jl/ik, we keep 1 and multiply V by 2.
inline void find_equivalent_OneBar_for_hamiltonian_generation(s4D<ValueType> Vijkl, std::vector<s4D<ValueType>>& v)
{
  v.reserve(24);
  v.clear();
  IndexType i, j, k, l;
  ValueType V;
  std::tie(i, j, k, l, V) = Vijkl;
  if (isComplex(std::get<4>(Vijkl)))
  {
    // <ij|kl> = <ji|lk> = <kl|ij>* = <lk|ji>*
    s3D<IndexType> ijkl = std::make_tuple(i, j, k, l);
    s3D<IndexType> jilk = std::make_tuple(j, i, l, k);
    s3D<IndexType> klij = std::make_tuple(k, l, i, j);
    s3D<IndexType> lkji = std::make_tuple(l, k, j, i);


    if (jilk != ijkl)
    {
      v.push_back(std::make_tuple(i, j, k, l, static_cast<RealType>(2.0) * V));
    }
    else
      v.push_back(Vijkl);

    if (klij != ijkl && klij != jilk)
    {
      // need to add klij as ma::conj(V). How do I add it?
      if (klij != lkji)
      { // add once with factor of 2
        v.push_back(std::make_tuple(k, l, i, j, ma::conj(static_cast<RealType>(2.0) * V)));
        if (ijkl == lkji)
        {
          std::cerr << " Error in find_equivalent_OneBar_for_hamiltonian_generation: Not sure how you got here. (ijkl "
                       "== lkji): "
                    << i << " " << j << " " << k << " " << l << "  " << V << std::endl;
          APP_ABORT("Error in find_equivalent_OneBar_for_hamiltonian_generation: Not sure how you got here. (ijkl == "
                    "lkji) \n");
        }
      }
      else
        v.push_back(std::make_tuple(k, l, i, j, ma::conj(V)));
    }
    else
    {
      // just checking
      if (klij == jilk && toComplex(V).imag() > 1e-8)
      {
        std::cerr << " Error in find_equivalent_OneBar_for_hamiltonian_generation: Not sure how you got here. (klij == "
                     "jilk): "
                  << i << " " << j << " " << k << " " << l << "  " << V << std::endl;
        APP_ABORT(
            "Error in find_equivalent_OneBar_for_hamiltonian_generation: Not sure how you got here. (klij == jilk) \n");
      }
    }
  }
  else
  {
    v.push_back(std::make_tuple(i, j, k, l, V));
    v.push_back(std::make_tuple(k, j, i, l, V));
    v.push_back(std::make_tuple(i, l, k, j, V));
    v.push_back(std::make_tuple(k, l, i, j, V));
    v.push_back(std::make_tuple(j, i, l, k, V));
    v.push_back(std::make_tuple(l, i, j, k, V));
    v.push_back(std::make_tuple(j, k, l, i, V));
    v.push_back(std::make_tuple(l, k, j, i, V));
    // order to remove consecutive repreated
    std::sort(v.begin(), v.end(), [](const s4D<ValueType>& lhs, const s4D<ValueType>& rhs) {
      return std::forward_as_tuple(std::get<0>(lhs), std::get<1>(lhs), std::get<2>(lhs), std::get<3>(lhs)) <
          std::forward_as_tuple(std::get<0>(rhs), std::get<1>(rhs), std::get<2>(rhs), std::get<3>(rhs));
    });
    // remove consecutive repeated
    std::vector<s4D<ValueType>>::iterator it;
    it = std::unique(v.begin(), v.end());
    // resize array, since unique does not remove elements
    // just reorders them
    v.resize(std::distance(v.begin(), it));

    //return;

    // look for symmetric terms
    it = v.begin();
    std::vector<s4D<ValueType>>::iterator it2;
    do
    {
      it2 = it + 1;
      while (it2 != v.end())
      {
        // <i1j1|k1l1> -> i1k1/j1l1, so look for
        //   i1==j2, k1==l2, j1==i2, l1==k2
        if (std::get<0>(*it) == std::get<1>(*it2) && std::get<2>(*it) == std::get<3>(*it2) &&
            std::get<1>(*it) == std::get<0>(*it2) && std::get<3>(*it) == std::get<2>(*it2))
        {
          it2 = v.erase(it2); // since it2 > it, we can safely erase
          std::get<4>(*it) *= static_cast<RealType>(2.0);
        }
        else
        {
          ++it2;
        }
      }
      it++;
    } while (it != v.end());
  }
}

inline void find_equivalent_TwoBar_for_integral_list(s4D<ValueType> Vijkl, std::vector<s4D<ValueType>>& v)
{
  v.reserve(24);
  v.clear();
  IndexType i, j, k, l;
  ValueType V;
  std::tie(i, j, k, l, V) = Vijkl;

  // <ij||kl> -> (ijkl)
  // Symmetries:
  // (ijkl) = (jilk) = (klij)* = (lkji)* = -(ijlk) = -(jikl) = -(lkij)* = -(klji)*
  // This method doesn't work here, because I modify V.
  //v.push_back(std::make_tuple(i,j,k,l,V));
  //v.push_back(std::make_tuple(j,i,l,k,V));
  //v.push_back(std::make_tuple(k,l,i,j,conj(V)));
  //v.push_back(std::make_tuple(l,k,j,i,conj(V)));
  //v.push_back(std::make_tuple(i,j,l,k,static_cast<RealType>(-1.0)*V));
  //v.push_back(std::make_tuple(j,i,k,l,static_cast<RealType>(-1.0)*V));
  //v.push_back(std::make_tuple(l,k,i,j,conj(static_cast<RealType>(-1.0)*V)));
  //v.push_back(std::make_tuple(k,l,j,i,conj(static_cast<RealType>(-1.0)*V)));

  auto ijkl = std::make_tuple(i, j, k, l);
  auto jilk = std::make_tuple(j, i, l, k);
  auto klij = std::make_tuple(k, l, i, j);
  auto lkji = std::make_tuple(l, k, j, i);
  auto ijlk = std::make_tuple(i, j, l, k);
  auto jikl = std::make_tuple(j, i, k, l);
  auto lkij = std::make_tuple(l, k, i, j);
  auto klji = std::make_tuple(k, l, j, i);

  // doing it by hand
  v.push_back(std::make_tuple(i, j, k, l, V));
  // slow and inefficient, but EASY to write!!!
  if (ijkl != jilk)
    v.push_back(std::make_tuple(j, i, l, k, V));
  if (klij != ijkl && klij != jilk)
    v.push_back(std::make_tuple(k, l, i, j, ma::conj(V)));
  if (lkji != ijkl && lkji != jilk && lkji != klij)
    v.push_back(std::make_tuple(l, k, j, i, ma::conj(V)));
  if (ijlk != lkji && ijlk != ijkl && ijlk != jilk && ijlk != klij)
    v.push_back(std::make_tuple(i, j, l, k, static_cast<RealType>(-1.0) * V));
  if (jikl != ijlk && jikl != lkji && jikl != ijkl && jikl != jilk && jikl != klij)
    v.push_back(std::make_tuple(j, i, k, l, static_cast<RealType>(-1.0) * V));
  if (lkij != jikl && lkij != ijlk && lkij != lkji && lkij != ijkl && lkij != jilk && lkij != klij)
    v.push_back(std::make_tuple(l, k, i, j, ma::conj(static_cast<RealType>(-1.0) * V)));
  if (klji != lkij && klji != jikl && klji != ijlk && klji != lkji && klji != ijkl && klji != jilk && klji != klij)
    v.push_back(std::make_tuple(k, l, j, i, ma::conj(static_cast<RealType>(-1.0) * V)));


  // just in case
  std::sort(v.begin(), v.end(), [](const s4D<ValueType>& lhs, const s4D<ValueType>& rhs) {
    return std::forward_as_tuple(std::get<0>(lhs), std::get<1>(lhs), std::get<2>(lhs), std::get<3>(lhs)) <
        std::forward_as_tuple(std::get<0>(rhs), std::get<1>(rhs), std::get<2>(rhs), std::get<3>(rhs));
  });
}

// to do:
//  3. create list of 2-bar integrals based on sparse storage of smallest element per symmetry set
inline void find_equivalent_TwoBar_for_hamiltonian_generation(s4D<ValueType> Vijkl, std::vector<s4D<ValueType>>& v)
{
  v.reserve(24);
  v.clear();
  IndexType i, j, k, l;
  ValueType V;
  std::tie(i, j, k, l, V) = Vijkl;

  // <ij||kl> -> (ijkl)
  // Symmetries:
  // (ijkl) = (jilk) = (klij)* = (lkji)* = -(ijlk) = -(jikl) = -(lkij)* = -(klji)*
  // This method doesn't work here, because I modify V.
  //v.push_back(std::make_tuple(i,j,k,l,V));
  //v.push_back(std::make_tuple(j,i,l,k,V));
  //v.push_back(std::make_tuple(k,l,i,j,conj(V)));
  //v.push_back(std::make_tuple(l,k,j,i,conj(V)));
  //v.push_back(std::make_tuple(i,j,l,k,static_cast<RealType>(-1.0)*V));
  //v.push_back(std::make_tuple(j,i,k,l,static_cast<RealType>(-1.0)*V));
  //v.push_back(std::make_tuple(l,k,i,j,conj(static_cast<RealType>(-1.0)*V)));
  //v.push_back(std::make_tuple(k,l,j,i,conj(static_cast<RealType>(-1.0)*V)));

  auto ijkl = std::make_tuple(i, j, k, l);
  auto jilk = std::make_tuple(j, i, l, k);
  auto klij = std::make_tuple(k, l, i, j);
  auto lkji = std::make_tuple(l, k, j, i);
  auto ijlk = std::make_tuple(i, j, l, k);
  auto jikl = std::make_tuple(j, i, k, l);
  auto lkij = std::make_tuple(l, k, i, j);
  auto klji = std::make_tuple(k, l, j, i);

  // slow and inefficient, but EASY to write!!!
  if (ijkl != jilk)
    v.push_back(std::make_tuple(i, j, k, l, static_cast<RealType>(2.0) * V));
  else
    v.push_back(std::make_tuple(i, j, k, l, V));

  bool t = false;
  if (klij != ijkl && klij != jilk)
  {
    t = true;
    v.push_back(std::make_tuple(k, l, i, j, ma::conj(V)));
  }
  if (lkji != ijkl && lkji != jilk && lkji != klij)
  {
    if (t)
      std::get<4>(v.back()) *= static_cast<RealType>(2.0);
    else
      v.push_back(std::make_tuple(l, k, j, i, ma::conj(V)));
  }

  t = false;
  if (ijlk != lkji && ijlk != ijkl && ijlk != jilk && ijlk != klij)
  {
    t = true;
    v.push_back(std::make_tuple(i, j, l, k, static_cast<RealType>(-1.0) * V));
  }
  if (jikl != ijlk && jikl != lkji && jikl != ijkl && jikl != jilk && jikl != klij)
  {
    if (t)
      std::get<4>(v.back()) *= static_cast<RealType>(2.0);
    else
      v.push_back(std::make_tuple(j, i, k, l, static_cast<RealType>(-1.0) * V));
  }

  t = false;
  if (lkij != jikl && lkij != ijlk && lkij != lkji && lkij != ijkl && lkij != jilk && lkij != klij)
  {
    t = true;
    v.push_back(std::make_tuple(l, k, i, j, ma::conj(static_cast<RealType>(-1.0) * V)));
  }
  if (klji != lkij && klji != jikl && klji != ijlk && klji != lkji && klji != ijkl && klji != jilk && klji != klij)
  {
    if (t)
      std::get<4>(v.back()) *= static_cast<RealType>(2.0);
    else
      v.push_back(std::make_tuple(k, l, j, i, ma::conj(static_cast<RealType>(-1.0) * V)));
  }

  std::sort(v.begin(), v.end(), [](const s4D<ValueType>& lhs, const s4D<ValueType>& rhs) {
    return std::forward_as_tuple(std::get<0>(lhs), std::get<1>(lhs), std::get<2>(lhs), std::get<3>(lhs)) <
        std::forward_as_tuple(std::get<0>(rhs), std::get<1>(rhs), std::get<2>(rhs), std::get<3>(rhs));
  });
}

// looks for the indexes of all equivalent terms in V2 and returns
// the smaller one (based on s4D ordering).
inline s4D<ValueType> find_smaller_equivalent_OneBar_for_integral_list(s4D<ValueType> ijkl)
{
  s4D<ValueType> v = ijkl;
  find_smallest_permutation(v);
  return v;
  /*
    std::vector<s4D<ValueType> > v;
    find_equivalent_OneBar_for_integral_list(ijkl,v);
    std::sort (v.begin(), v.end(),
      [](const s4D<ValueType>& lhs, const s4D<ValueType>& rhs)
      {
        return std::forward_as_tuple(std::get<0>(lhs),std::get<1>(lhs),std::get<2>(lhs),std::get<3>(lhs)) < std::forward_as_tuple(std::get<0>(rhs),std::get<1>(rhs),std::get<2>(rhs),std::get<3>(rhs));
      }
    );
    return v[0];
*/
}

// looks for the indexes of all equivalent terms in V2_2bar and returns
// the smaller one (based on s4D ordering).
inline s4D<ValueType> find_smaller_equivalent_TwoBar_for_integral_list(s4D<ValueType> ijkl)
{
  std::vector<s4D<ValueType>> v;
  find_equivalent_TwoBar_for_integral_list(ijkl, v);
  std::sort(v.begin(), v.end(), [](const s4D<ValueType>& lhs, const s4D<ValueType>& rhs) {
    return std::forward_as_tuple(std::get<0>(lhs), std::get<1>(lhs), std::get<2>(lhs), std::get<3>(lhs)) <
        std::forward_as_tuple(std::get<0>(rhs), std::get<1>(rhs), std::get<2>(rhs), std::get<3>(rhs));
  });
  return v[0];
}

} // namespace afqmc
} // namespace qmcplusplus
#endif
