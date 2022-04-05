//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/*
 DO NOT MAKE PERMANENT EDITS IN THIS FILE
 This file is generated from src/QMCWavefunctions/Jastrow/codegen/user_jastrow.py and UserFunctor.h.in

 To make changes, edit UserFunctor.h.in and rerun user_jastrow.py.
*/


/** @file UserFunctor.h
 * @brief User-defined functor
 */
#ifndef QMCPLUSPLUS_USERFUNCTOR_H
#define QMCPLUSPLUS_USERFUNCTOR_H
#include "OptimizableFunctorBase.h"
#include "OhmmsData/AttributeSet.h"
#include <cmath>
// #include <vector>
#include "OhmmsPETE/TinyVector.h"
#include "OMPTarget/OffloadAlignedAllocators.hpp"


namespace qmcplusplus
{
/** Implements the function
 *  f = A*r/(B*r + 1) - A/B

 *
 */
template<class T>
struct UserFunctor : public OptimizableFunctorBase
{
  /// Is optimizable
  bool Opt_A;
  /// Value
  real_type A;
  /// Id in XML input
  std::string ID_A;


  /// Is optimizable
  bool Opt_B;
  /// Value
  real_type B;
  /// Id in XML input
  std::string ID_B;


  ///default constructor
  UserFunctor() { reset(); }

  // void setCusp(real_type cusp)

  void setCusp(real_type cusp) override
  {
    A     = cusp;
    Opt_A = false;
    reset();
  }


  OptimizableFunctorBase* makeClone() const override { return new UserFunctor(*this); }

  void reset() override {}


  // inline real_type evaluate(real_type r) const

  inline real_type evaluate(real_type r) const { return A * r / (B * r + 1) - A / B; }


  // const inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2) const

  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2) const
  {
    dudr   = -A * B * r / ((B * r + 1) * (B * r + 1)) + A / (B * r + 1);
    d2udr2 = 2 * A * B * (B * r / (B * r + 1) - 1) / ((B * r + 1) * (B * r + 1));
    return A * r / (B * r + 1) - A / B;
  }


  // inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2, real_type& d3udr3) const

  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2, real_type& d3udr3) const
  {
    dudr   = -A * B * r / ((B * r + 1) * (B * r + 1)) + A / (B * r + 1);
    d2udr2 = 2 * A * B * (B * r / (B * r + 1) - 1) / ((B * r + 1) * (B * r + 1));
    d3udr3 = 6 * A * ((B) * (B)) * (-B * r / (B * r + 1) + 1) / std::pow(B * r + 1, 3);
    return A * r / (B * r + 1) - A / B;
  }


  inline real_type evaluateV(const int iat,
                             const int iStart,
                             const int iEnd,
                             const T* restrict _distArray,
                             T* restrict distArrayCompressed) const
  {
    // specialized evaluation loop?
    real_type sum(0);
    for (int idx = iStart; idx < iEnd; idx++)
      if (idx != iat)
        sum += evaluate(_distArray[idx]);
    return sum;
  }

  /** evaluate sum of the pair potentials FIXME
   * @return \f$\sum u(r_j)\f$ for r_j < cutoff_radius
   */
  static void mw_evaluateV(const int num_groups,
                           const UserFunctor* const functors[],
                           const int n_src,
                           const int* grp_ids,
                           const int num_pairs,
                           const int* ref_at,
                           const T* mw_dist,
                           const int dist_stride,
                           T* mw_vals,
                           Vector<char, OffloadPinnedAllocator<char>>& transfer_buffer)
  {
    for (int ip = 0; ip < num_pairs; ip++)
    {
      mw_vals[ip] = 0;
      for (int j = 0; j < n_src; j++)
      {
        const int ig = grp_ids[j];
        auto& functor(*functors[ig]);
        if (j != ref_at[ip])
          mw_vals[ip] += functor.evaluate(mw_dist[ip * dist_stride + j]);
      }
    }
  }

  inline void evaluateVGL(const int iat,
                          const int iStart,
                          const int iEnd,
                          const T* distArray,
                          T* restrict valArray,
                          T* restrict gradArray,
                          T* restrict laplArray,
                          T* restrict distArrayCompressed,
                          int* restrict distIndices) const
  {
    // specialized evaluation loop?
    for (int idx = iStart; idx < iEnd; idx++)
    {
      valArray[idx] = evaluate(distArray[idx], gradArray[idx], laplArray[idx]);
      gradArray[idx] /= distArray[idx];
    }
    if (iat >= iStart && iat < iEnd)
      valArray[iat] = gradArray[iat] = laplArray[iat] = T(0);
  }

  static void mw_evaluateVGL(const int iat,
                             const int num_groups,
                             const UserFunctor* const functors[],
                             const int n_src,
                             const int* grp_ids,
                             const int nw,
                             T* mw_vgl, // [nw][DIM+2]
                             const int n_padded,
                             const T* mw_dist, // [nw][DIM+1][n_padded]
                             T* mw_cur_allu,   // [nw][3][n_padded]
                             Vector<char, OffloadPinnedAllocator<char>>& transfer_buffer)
  {
    throw std::runtime_error("UserFunctor mw_evaluateVGL not implemented!");
  }

  inline real_type f(real_type r) override { return evaluate(r); }

  inline real_type df(real_type r) override
  {
    real_type dudr, d2udr2;
    real_type res = evaluate(r, dudr, d2udr2);
    return dudr;
  }

  static void mw_updateVGL(const int iat,
                           const std::vector<bool>& isAccepted,
                           const int num_groups,
                           const UserFunctor* const functors[],
                           const int n_src,
                           const int* grp_ids,
                           const int nw,
                           T* mw_vgl, // [nw][DIM+2]
                           const int n_padded,
                           const T* mw_dist, // [nw][DIM+1][n_padded]
                           T* mw_allUat,     // [nw][DIM+2][n_padded]
                           T* mw_cur_allu,   // [nw][3][n_padded]
                           Vector<char, OffloadPinnedAllocator<char>>& transfer_buffer)
  {
    throw std::runtime_error("UserFunctor mw_updateVGL not implemented!");
  }

  // inline bool evaluateDerivatives(real_type r, std::vector<TinyVector<real_type, 3>>& derivs)

  inline bool evaluateDerivatives(real_type r, std::vector<TinyVector<real_type, 3>>& derivs) override
  {
    int i = 0;

    if (Opt_A)
    {
      derivs[i][0] = r / (B * r + 1) - 1 / B;
      derivs[i][1] = -B * r / ((B * r + 1) * (B * r + 1)) + 1.0 / (B * r + 1);
      derivs[i][2] = 2 * B * (B * r / (B * r + 1) - 1) / ((B * r + 1) * (B * r + 1));
      ++i;
    }

    if (Opt_B)
    {
      derivs[i][0] = -A * ((r) * (r)) / ((B * r + 1) * (B * r + 1)) + A / ((B) * (B));
      derivs[i][1] = 2 * A * B * ((r) * (r)) / std::pow(B * r + 1, 3) - 2 * A * r / ((B * r + 1) * (B * r + 1));
      derivs[i][2] = -4 * A * B * r * (B * r / (B * r + 1) - 1) / std::pow(B * r + 1, 3) +
          2 * A * B * (-B * ((r) * (r)) / ((B * r + 1) * (B * r + 1)) + r / (B * r + 1)) / ((B * r + 1) * (B * r + 1)) +
          2 * A * (B * r / (B * r + 1) - 1) / ((B * r + 1) * (B * r + 1));
      ++i;
    }

    return true;
  }


  /// compute derivatives with respect to variational parameters
  inline bool evaluateDerivatives(real_type r, std::vector<real_type>& derivs) override
  {
    int i = 0;

    if (Opt_A)
    {
      derivs[i] = r / (B * r + 1) - 1 / B;
      ++i;
    }

    if (Opt_B)
    {
      derivs[i] = -A * ((r) * (r)) / ((B * r + 1) * (B * r + 1)) + A / ((B) * (B));
      ++i;
    }


    return true;
  }


  //  bool put(xmlNodePtr cur)

  bool put(xmlNodePtr cur) override
  {
    cur = cur->xmlChildrenNode;
    while (cur != NULL)
    {
      std::string cname((const char*)(cur->name));
      if (cname == "var") //only accept var
      {
        std::string id_in;
        std::string p_name;
        OhmmsAttributeSet rAttrib;
        rAttrib.add(id_in, "id");
        rAttrib.add(p_name, "name");
        rAttrib.put(cur);

        if (p_name == "A")
        {
          ID_A = id_in;
          putContent(A, cur);
          Opt_A = true;
        }

        if (p_name == "B")
        {
          ID_B = id_in;
          putContent(B, cur);
          Opt_B = true;
        }
      }
      cur = cur->next;
    }
    reset();
    myVars.clear();

    if (Opt_A)
      myVars.insert(ID_A, A, Opt_A, optimize::OTHER_P);

    if (Opt_B)
      myVars.insert(ID_B, B, Opt_B, optimize::OTHER_P);

    return true;
  }


  void checkInVariables(opt_variables_type& active) override
  {
    active.insertFrom(myVars);
    //myVars.print(std::cout);
  }

  void checkOutVariables(const opt_variables_type& active) override
  {
    myVars.getIndex(active);
    //myVars.print(std::cout);
  }

  //void resetParameters(const opt_variables_type& active)

  void resetParameters(const opt_variables_type& active) override
  {
    if (myVars.size())
    {
      int ia = myVars.where(0);
      if (ia > -1)
      {
        int i = 0;

        if (Opt_A)
          A = std::real(myVars[i++] = active[ia++]);

        if (Opt_B)
          B = std::real(myVars[i++] = active[ia++]);
      }
      reset();
    }
  }
};


} // namespace qmcplusplus
#endif
