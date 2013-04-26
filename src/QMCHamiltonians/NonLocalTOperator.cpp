//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file NonLocalTOperator.cpp
 *@brief Definition of NonLocalTOperator
 */
#include "QMCHamiltonians/NonLocalTOperator.h"
#include "OhmmsData/ParameterSet.h"

namespace qmcplusplus
{

NonLocalTOperator::NonLocalTOperator():Tau(0.01),Alpha(0.0),Gamma(0.0)
{
}

/** process options related to TMoves
 * @return true, if TMove is used.
 */
bool NonLocalTOperator::put(xmlNodePtr cur)
{
  string use_tmove="no";
  ParameterSet m_param;
  m_param.add(Tau,"timeStep","double");
  m_param.add(Tau,"timestep","double");
  m_param.add(Tau,"Tau","double");
  m_param.add(Tau,"tau","double");
  m_param.add(Alpha,"alpha","double");
  m_param.add(Gamma,"gamma","double");
  m_param.add(use_tmove,"nonlocalmove","double");
  m_param.add(use_tmove,"nonlocalmoves","double");
  bool success = m_param.put(cur);
  plusFactor=Tau*Gamma;
  minusFactor=-Tau*(1.0-Alpha*(1.0+Gamma));
  return use_tmove=="yes";
}

void NonLocalTOperator::reset()
{
  Txy.erase(Txy.begin(),Txy.end());
  Txy.push_back(NonLocalData());
}

void NonLocalTOperator::reserve(int n)
{
  Txy.reserve(n);
  Txy.push_back(NonLocalData());
}

int NonLocalTOperator::selectMove(RealType prob)
{
  RealType wgt_t=1.0;
  for(int i=1; i<Txy.size(); i++)
  {
    if(Txy[i].Weight>0)
    {
      wgt_t += Txy[i].Weight *=plusFactor;
    }
    else
    {
      wgt_t += Txy[i].Weight *=minusFactor;
    }
  }
  prob *= wgt_t;
  RealType wsum=Txy[0].Weight;
  int ibar=0;;
  while(wsum<prob)
  {
    ibar++;
    wsum += Txy[ibar].Weight;
  }
  return ibar;
}

int NonLocalTOperator::selectMove(RealType prob,
                                  vector<NonLocalData> &txy)
{
  RealType wgt_t=1.0;
  for(int i=1; i<txy.size(); i++)
  {
    if(txy[i].Weight>0)
    {
      wgt_t += txy[i].Weight *=plusFactor;
    }
    else
    {
      wgt_t += txy[i].Weight *=minusFactor;
    }
  }
  prob *= wgt_t;
  RealType wsum=txy[0].Weight;
  int ibar=0;;
  while(wsum<prob)
  {
    ibar++;
    wsum += txy[ibar].Weight;
  }
  return ibar;
}


}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

