//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/**@file HamiltonianPool.cpp
 * @brief Implements HamiltonianPool operators.
 */
#include "QMCApp/HamiltonianPool.h"
#include "QMCApp/WaveFunctionPool.h"
#include "QMCApp/ParticleSetPool.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/OpenMP.h"
#include "Utilities/ProgressReportEngine.h"

namespace qmcplusplus
{

HamiltonianPool::HamiltonianPool(Communicate* c, const char* aname)
  : MPIObjectBase(c), curH(0), ptclPool(0), psiPool(0), curDoc(0)
{
  ClassName="HamiltonianPool";
  myName=aname;
}

bool HamiltonianPool::put(xmlNodePtr cur)
{
  ReportEngine PRE("HamiltonianPool","put");
  std::string id("h0"), target("e"),role("extra");
  OhmmsAttributeSet hAttrib;
  hAttrib.add(id,"id");
  hAttrib.add(id,"name");
  hAttrib.add(role,"role");
  hAttrib.add(target,"target");
  hAttrib.put(cur);
  ParticleSet* qp=ptclPool->getParticleSet(target);
  if(qp == 0)
  {
    //never a good thing
    PRE.error("No target particle "+ target+ " exists.");
    return false;
  }
  bool set2Primary=false;
  //first Hamiltonian is set to the primary Hamiltonian
  if(myPool.empty() || role == "primary" )
    set2Primary=true;
  HamiltonianFactory *curH=0;
  PoolType::iterator hit(myPool.find(id));
  if(hit == myPool.end())
  {
    curH= new HamiltonianFactory(qp, ptclPool->getPool(), psiPool->getPool(),myComm);
    curH->setName(id);
    myPool[id]=curH;
  }
  else
    curH=(*hit).second;
  bool success= curH->put(cur);
  if(set2Primary)
    primaryH=curH->targetH;
  return success;
}

bool HamiltonianPool::get(std::ostream& os) const
{
  PoolType::const_iterator it(myPool.begin()), it_end(myPool.end());
  while(it != it_end)
  {
    os << "\n  Hamiltonian " << (*it).first << std::endl;;
    (*it).second->targetH->get(os);
    ++it;
  }
  return true;
}

void HamiltonianPool::reset()
{
}
}
