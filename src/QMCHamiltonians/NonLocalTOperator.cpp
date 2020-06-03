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
#include "QMCHamiltonians/NonLocalTOperator.h"
#include "OhmmsData/ParameterSet.h"

namespace qmcplusplus
{
NonLocalTOperator::NonLocalTOperator(size_t N) : Tau(0.01), Alpha(0.0), Gamma(0.0), Nelec(N) {}

/** process options related to TMoves
 * @return Tmove version
 *  Turns out this wants the NodePtr to the entire driver block.
 */
int NonLocalTOperator::put(xmlNodePtr cur)
{
  std::string use_tmove = "no";
  ParameterSet m_param;
  m_param.add(Tau, "timeStep", "double");
  m_param.add(Tau, "timestep", "double");
  m_param.add(Tau, "Tau", "double");
  m_param.add(Tau, "tau", "double");
  m_param.add(Alpha, "alpha", "double");
  m_param.add(Gamma, "gamma", "double");
  m_param.add(use_tmove, "nonlocalmove", "double");
  m_param.add(use_tmove, "nonlocalmoves", "double");
  bool success = m_param.put(cur);
  plusFactor   = Tau * Gamma;
  minusFactor  = -Tau * (1.0 - Alpha * (1.0 + Gamma));
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
  {
    APP_ABORT("NonLocalTOperator::put unknown nonlocalmove option " + use_tmove);
  }
#pragma omp master
  app_log() << o.str() << std::endl;
  return v_tmove;
}

int NonLocalTOperator::thingsThatShouldBeInMyConstructor(const std::string& non_local_move_option,
                                                         const double tau,
                                                         const double alpha,
                                                         const double gamma)
{
  Tau         = tau;
  Alpha       = alpha;
  Gamma       = gamma;
  plusFactor  = Tau * Gamma;
  minusFactor = -Tau * (1.0 - Alpha * (1.0 + Gamma));
  int v_tmove = TMOVE_OFF;
  std::ostringstream o;

  if (non_local_move_option == "no")
  {
    v_tmove = TMOVE_OFF;
    o << "  Using Locality Approximation";
  }
  else if (non_local_move_option == "yes" || non_local_move_option == "v0")
  {
    v_tmove = TMOVE_V0;
    o << "  Using Non-local T-moves v0, M. Casula, PRB 74, 161102(R) (2006)";
  }
  else if (non_local_move_option == "v1")
  {
    v_tmove = TMOVE_V1;
    o << "  Using Non-local T-moves v1, M. Casula et al., JCP 132, 154113 (2010)";
  }
  else if (non_local_move_option == "v3")
  {
    v_tmove = TMOVE_V3;
    o << "  Using Non-local T-moves v3, an approximation to v1";
  }
  else
  {
    APP_ABORT("NonLocalTOperator::put unknown nonlocalmove option " + non_local_move_option);
  }
  app_log() << o.str() << std::endl;
  return v_tmove;
}
void NonLocalTOperator::reset() { Txy.clear(); }

const NonLocalData* NonLocalTOperator::selectMove(RealType prob, std::vector<NonLocalData>& txy) const
{
  RealType wgt_t = 1.0;
  for (int i = 0; i < txy.size(); i++)
  {
    if (txy[i].Weight > 0)
    {
      wgt_t += txy[i].Weight *= plusFactor;
    }
    else
    {
      wgt_t += txy[i].Weight *= minusFactor;
    }
  }
  prob *= wgt_t;
  RealType wsum = 1.0;
  int ibar      = 0;
  while (wsum < prob)
  {
    wsum += txy[ibar].Weight;
    ibar++;
  }
  return ibar > 0 ? &(txy[ibar - 1]) : nullptr;
}

void NonLocalTOperator::group_by_elec()
{
  Txy_by_elec.resize(Nelec);
  for (int i = 0; i < Nelec; i++)
  {
    Txy_by_elec[i].clear();
  }

  for (int i = 0; i < Txy.size(); i++)
  {
    Txy_by_elec[Txy[i].PID].push_back(Txy[i]);
  }
}

} // namespace qmcplusplus
