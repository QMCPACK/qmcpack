//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/**@file NonLocalTOperator.h
 * @brief Declaration of NonLocalTOperator
 *
 * NonLocalTOperator has the off-diagonal transition probability matrix.
 */
#ifndef QMCPLUSPLUS_NONLOCALTRANSITIONOPERATOR_H
#define QMCPLUSPLUS_NONLOCALTRANSITIONOPERATOR_H

#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace qmcplusplus
{

struct NonLocalTOperator
{

  typedef NonLocalData::RealType RealType;
  typedef NonLocalData::PosType PosType;

  RealType Tau;
  RealType Alpha;
  RealType Gamma;
  RealType plusFactor;
  RealType minusFactor;

  std::vector<NonLocalData> Txy;

  NonLocalTOperator();
  inline int size() const
  {
    return Txy.size();
  }

  inline int id(int ibar) const
  {
    return Txy[ibar].PID;
  }

  inline PosType delta(int ibar) const
  {
    return Txy[ibar].Delta;
  }

  /** initialize the parameters */
  bool put(xmlNodePtr cur);

  /** reserve Txy for memory optimization */
  void reserve(int n);

  /** reset Txy for a new set of non-local moves
   *
   * Txy[0] is always 1 corresponding to the diagonal(no) move
   */
  void reset();

  /** select the move for a given probability
   * @param prob value [0,1)
   * @return the move index k for \f$\sum_i^K T/\sum_i^N < prob\f$
   */
  int selectMove(RealType prob);
  int selectMove(RealType prob, std::vector<NonLocalData> &txy);
};

}
#endif


