//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file FakeFunctor.h
 * @brief FakeFunctor for testing purpose
 */
#ifndef QMCPLUSPLUS_FAKEFUNCTOR_H
#define QMCPLUSPLUS_FAKEFUNCTOR_H

#include "OptimizableFunctorBase.h"
#include <cmath>
#include "OhmmsPETE/TinyVector.h"
#include "OMPTarget/OffloadAlignedAllocators.hpp"

namespace qmcplusplus
{
template<class T>
struct FakeFunctor : public OptimizableFunctorBase
{
  ///default constructor
  FakeFunctor(const std::string& my_name) : OptimizableFunctorBase(my_name) {}

  void setCusp(real_type cusp) override {}

  OptimizableFunctorBase* makeClone() const override { return new FakeFunctor(*this); }

  void reset() override {}

  inline real_type evaluate(real_type r) const { return 0.0; }


  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2) const
  {
    dudr   = 0.0;
    d2udr2 = 0.0;
    return 0.0;
  }

  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2, real_type& d3udr3) const
  {
    dudr   = 0.0;
    d2udr2 = 0.0;
    d3udr3 = 0.0;
    return 0.0;
  }


  inline real_type evaluateV(const int iat,
                             const int iStart,
                             const int iEnd,
                             const T* restrict _distArray,
                             T* restrict distArrayCompressed) const
  {
    return 0.0;
  }

  static void mw_evaluateV(const int num_groups,
                           const FakeFunctor* const functors[],
                           const int n_src,
                           const int* grp_ids,
                           const int num_pairs,
                           const int* ref_at,
                           const T* mw_dist,
                           const int dist_stride,
                           T* mw_vals,
                           Vector<char, OffloadPinnedAllocator<char>>& transfer_buffer)
  {}

  inline void evaluateVGL(const int iat,
                          const int iStart,
                          const int iEnd,
                          const T* distArray,
                          T* restrict valArray,
                          T* restrict gradArray,
                          T* restrict laplArray,
                          T* restrict distArrayCompressed,
                          int* restrict distIndices) const
  {}

  static void mw_evaluateVGL(const int iat,
                             const int num_groups,
                             const FakeFunctor* const functors[],
                             const int n_src,
                             const int* grp_ids,
                             const int nw,
                             T* mw_vgl, // [nw][DIM+2]
                             const int n_padded,
                             const T* mw_dist, // [nw][DIM+1][n_padded]
                             T* mw_cur_allu,   // [nw][3][n_padded]
                             Vector<char, OffloadPinnedAllocator<char>>& transfer_buffer)
  {}

  inline real_type f(real_type r) override { return evaluate(r); }

  inline real_type df(real_type r) override { return 0.0; }

  static void mw_updateVGL(const int iat,
                           const std::vector<bool>& isAccepted,
                           const int num_groups,
                           const FakeFunctor* const functors[],
                           const int n_src,
                           const int* grp_ids,
                           const int nw,
                           T* mw_vgl, // [nw][DIM+2]
                           const int n_padded,
                           const T* mw_dist, // [nw][DIM+1][n_padded]
                           T* mw_allUat,     // [nw][DIM+2][n_padded]
                           T* mw_cur_allu,   // [nw][3][n_padded]
                           Vector<char, OffloadPinnedAllocator<char>>& transfer_buffer)
  {}

  inline bool evaluateDerivatives(real_type r, std::vector<TinyVector<real_type, 3>>& derivs) override
  {
    derivs = derivs_;
    return true;
  }


  /// compute derivatives with respect to variational parameters
  inline bool evaluateDerivatives(real_type r, std::vector<real_type>& derivs) override { return true; }

  bool put(xmlNodePtr cur) override { return true; }


  void checkInVariablesExclusive(opt_variables_type& active) override {}

  void checkOutVariables(const opt_variables_type& active) override {}

  void resetParametersExclusive(const opt_variables_type& active) override {}


  std::vector<TinyVector<real_type, 3>> derivs_;
};


} // namespace qmcplusplus
#endif
