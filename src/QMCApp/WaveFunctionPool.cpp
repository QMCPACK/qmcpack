//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/**@file WaveFunctionPool.cpp
 * @brief Implements WaveFunctionPool operators.
 */
#include "QMCApp/WaveFunctionPool.h"
#include "QMCApp/ParticleSetPool.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{


WaveFunctionPool::WaveFunctionPool(Communicate* c, const char* aname)
  : MPIObjectBase(c)
{
  ClassName="WaveFunctionPool";
  myName=aname;
}

WaveFunctionPool::~WaveFunctionPool()
{
  DEBUG_MEMORY("WaveFunctionPool::~WaveFunctionPool");
  PoolType::iterator it(myPool.begin());
  while(it != myPool.end())
  {
    delete (*it).second;
    ++it;
  }
}

bool WaveFunctionPool::put(xmlNodePtr cur)
{
  std::string id("psi0"), target("e"), role("extra");
  OhmmsAttributeSet pAttrib;
  pAttrib.add(id,"id");
  pAttrib.add(id,"name");
  pAttrib.add(target,"target");
  pAttrib.add(target,"ref");
  pAttrib.add(role,"role");
  pAttrib.put(cur);
  ParticleSet *qp = ptclPool->getParticleSet(target);
 
  {//check ESHDF should be used to initialize both target and associated ionic system
    xmlNodePtr tcur=cur->children;
    while(tcur != NULL)
    { //check <determinantset/> or <sposet_builder/> to extract the ionic and electronic structure
      std::string cname((const char*)tcur->name);
      if(cname == OrbitalBuilderBase::detset_tag || cname =="sposet_builder")
      { 
        qp=ptclPool->createESParticleSet(tcur,target,qp);
      }
      tcur=tcur->next;
    }
  }
  if(qp==0)
  {
    APP_ABORT("WaveFunctionPool::put Target ParticleSet is not found.");
  }
  std::map<std::string,WaveFunctionFactory*>::iterator pit(myPool.find(id));
  WaveFunctionFactory* psiFactory=0;
  bool isPrimary=true;
  if(pit == myPool.end())
  {
    psiFactory=new WaveFunctionFactory(qp,ptclPool->getPool(),myComm);
    psiFactory->setName(id);
    isPrimary = (myPool.empty() || role == "primary");
    myPool[id]=psiFactory;
    app_summary() << " Wavefunction setup: " << std::endl;
    app_summary() << " ------------------- " << std::endl;
    app_summary() << "  Name: " << psiFactory->getName() << std::endl;

  }
  else
  {
    psiFactory=(*pit).second;
  }
  bool success = psiFactory->put(cur);
  if(success && isPrimary)
  {
    primaryPsi=psiFactory->targetPsi;
  }
  return success;
}

void  WaveFunctionPool::addFactory(WaveFunctionFactory* psifac)
{
  PoolType::iterator oit(myPool.find(psifac->getName()));
  if(oit == myPool.end())
  {
    app_log() << "  Adding " << psifac->getName() << " WaveFunctionFactory to the pool" << std::endl;
    myPool[psifac->getName()]=psifac;
  }
  else
  {
    app_warning() << "  " << psifac->getName() << " exists. Ignore addition" << std::endl;
  }
}

xmlNodePtr WaveFunctionPool::getWaveFunctionNode(const std::string& id)
{
  if(myPool.empty())
    return NULL;
  std::map<std::string,WaveFunctionFactory*>::iterator it(myPool.find(id));
  if(it == myPool.end())
  {
    return (*myPool.begin()).second->myNode;
  }
  else
  {
    return (*it).second->myNode;
  }
}
}
