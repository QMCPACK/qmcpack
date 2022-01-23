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

#include "NonLocalData.h"

namespace qmcplusplus
{
/// Tmove options
enum
{
  TMOVE_OFF = 0, // no Tmove
  TMOVE_V0,      // M. Casula, PRB 74, 161102(R) (2006)
  TMOVE_V1,      // version 1, M. Casula et al., JCP 132, 154113 (2010)
  TMOVE_V3,      // an approximation to version 1 but much faster.
};

struct NonLocalTOperator
{
  using RealType = NonLocalData::RealType;
  using PosType  = NonLocalData::PosType;

  NonLocalTOperator();

  /** replacement for put because wouldn't it be cool to know what the classes configuration actually
   *  is.
   */
  int thingsThatShouldBeInMyConstructor(const std::string& non_local_move_option,
                                        const double tau,
                                        const double alpha,
                                        const double gamma);
  /** initialize the parameters */
  int put(xmlNodePtr cur);

  /** select the move for a given probability
   * @param prob value [0,1)
   * @param txy a given Txy collection
   * @return pointer to NonLocalData
   */
  const NonLocalData* selectMove(RealType prob, const std::vector<NonLocalData>& txy);

  /** select the move for a given probability using internal txy_by_elec_
   * @param prob value [0,1)
   * @param iel reference electron
   * @return pointer to NonLocalData
   */
  inline const NonLocalData* selectMove(RealType prob, int iel) { return selectMove(prob, txy_by_elec_[iel]); }

  /** sort all the Txy elements by electron */
  void groupByElectron(size_t num_elec, const std::vector<NonLocalData>& txy);

private:
  RealType tau_;
  RealType alpha_;
  RealType gamma_;
  /// factor applied on >0 weight
  RealType plusFactor;
  /// factor applied on <=0 weight
  RealType minusFactor;
  // for selecting a move
  std::vector<RealType> txy_scan_;
  // txy grouped by electron id
  std::vector<std::vector<NonLocalData>> txy_by_elec_;
};

} // namespace qmcplusplus
#endif
