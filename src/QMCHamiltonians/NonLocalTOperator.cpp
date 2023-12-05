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


/**@file NonLocalTOperator.cpp
 *@brief Definition of NonLocalTOperator
 */
#include "NonLocalTOperator.h"
#include "OhmmsData/ParameterSet.h"

namespace qmcplusplus
{
NonLocalTOperator::NonLocalTOperator() : tau_(0.01), alpha_(0.0), gamma_(0.0) {}

/** process options related to TMoves
 * @return Tmove version
 *  Turns out this wants the NodePtr to the entire driver block.
 */
int NonLocalTOperator::put(xmlNodePtr cur)
{
  std::string use_tmove = "no";
  ParameterSet m_param;
  m_param.add(tau_, "timeStep");
  m_param.add(tau_, "timestep");
  m_param.add(tau_, "Tau");
  m_param.add(tau_, "tau");
  m_param.add(alpha_, "alpha");
  m_param.add(gamma_, "gamma");
  m_param.add(use_tmove, "nonlocalmove");
  m_param.add(use_tmove, "nonlocalmoves");
  bool success = m_param.put(cur);
  plusFactor   = tau_ * gamma_;
  minusFactor  = -tau_ * (1.0 - alpha_ * (1.0 + gamma_));
  int v_tmove  = TMOVE_OFF;
  std::ostringstream o;
  if (use_tmove == "no")
  {
    v_tmove = TMOVE_OFF;
    o << "  Using Locality Approximation";
  }
  else if (use_tmove == "yes" || use_tmove == "v0")
  {
    v_tmove = TMOVE_V0;
    o << "  Using Non-local T-moves v0, M. Casula, PRB 74, 161102(R) (2006)";
  }
  else if (use_tmove == "v1")
  {
    v_tmove = TMOVE_V1;
    o << "  Using Non-local T-moves v1, M. Casula et al., JCP 132, 154113 (2010)";
  }
  else if (use_tmove == "v3")
  {
    v_tmove = TMOVE_V3;
    o << "  Using Non-local T-moves v3, an approximation to v1";
  }
  else
    throw std::runtime_error("NonLocalTOperator::put unknown nonlocalmove option " + use_tmove);

#pragma omp master
  app_log() << o.str() << std::endl;

  return v_tmove;
}

int NonLocalTOperator::thingsThatShouldBeInMyConstructor(const std::string& non_local_move_option,
                                                         const double tau,
                                                         const double alpha,
                                                         const double gamma)
{
  tau_        = tau;
  alpha_      = alpha;
  gamma_      = gamma;
  plusFactor  = tau_ * gamma_;
  minusFactor = -tau_ * (1.0 - alpha_ * (1.0 + gamma_));
  int v_tmove = TMOVE_OFF;

  if (non_local_move_option == "no")
    v_tmove = TMOVE_OFF;
  else if (non_local_move_option == "yes" || non_local_move_option == "v0")
    v_tmove = TMOVE_V0;
  else if (non_local_move_option == "v1")
    v_tmove = TMOVE_V1;
  else if (non_local_move_option == "v3")
    v_tmove = TMOVE_V3;
  else
    throw std::runtime_error("NonLocalTOperator::put unknown nonlocalmove option " + non_local_move_option);
  return v_tmove;
}

const NonLocalData* NonLocalTOperator::selectMove(RealType prob, const std::vector<NonLocalData>& txy)
{
  // txy_scan_[0] = 1.0, txy_scan_[i>0] = txy_scan_[i-1] + txy[i-1].Weight (modified)
  txy_scan_.resize(txy.size());
  RealType wgt_t = 1.0;
  for (int i = 0; i < txy.size(); i++)
  {
    txy_scan_[i] = wgt_t;
    if (txy[i].Weight > 0)
      wgt_t += txy[i].Weight * plusFactor;
    else
      wgt_t += txy[i].Weight * minusFactor;
  }

  const RealType target = prob * wgt_t;
  // find ibar which satisify txy_scan_[ibar-1] <= target < txy_scan_[ibar]
  int ibar = 0;
  while (ibar < txy_scan_.size() && txy_scan_[ibar] <= target)
    ibar++;

  return ibar > 0 ? &(txy[ibar - 1]) : nullptr;
}

void NonLocalTOperator::groupByElectron(size_t num_elec, const std::vector<NonLocalData>& txy)
{
  txy_by_elec_.resize(num_elec);
  for (int i = 0; i < num_elec; i++)
    txy_by_elec_[i].clear();

  for (int i = 0; i < txy.size(); i++)
  {
    assert(txy[i].PID >= 0 && txy[i].PID < num_elec);
    txy_by_elec_[txy[i].PID].push_back(txy[i]);
  }
}

} // namespace qmcplusplus
