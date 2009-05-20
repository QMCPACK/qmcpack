//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file WaveFunctionPool.cpp
 * @brief Implements WaveFunctionPool operators.
 */
#include "QMCApp/WaveFunctionPool.h"
#include "QMCApp/ParticleSetPool.h"
using namespace std;
#include "OhmmsData/AttributeSet.h"
#include "Utilities/OhmmsInfo.h"

namespace qmcplusplus {
  

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

  bool WaveFunctionPool::put(xmlNodePtr cur) {

    string id("psi0"), target("e"), role("extra");
    OhmmsAttributeSet pAttrib;
    pAttrib.add(id,"id"); pAttrib.add(id,"name"); 
    pAttrib.add(target,"target"); pAttrib.add(target,"ref"); 
    pAttrib.add(role,"role");
    pAttrib.put(cur);

    ParticleSet *qp = ptclPool->getParticleSet(target);
    if(qp == 0) {
      ERRORMSG("Wavefunction cannot be created because of missing particle set " << target)
      return false;
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
    } else {
      psiFactory=(*pit).second;
    }

    bool success = psiFactory->put(cur);
    if(success && isPrimary) {
      primaryPsi=psiFactory->targetPsi;
    }

    return success;
  }

  void  WaveFunctionPool::addFactory(WaveFunctionFactory* psifac) {
    PoolType::iterator oit(myPool.find(psifac->getName()));
    if(oit == myPool.end()) {
      LOGMSG("  Adding " << psifac->getName() << " WaveFunctionFactory to the pool")
      myPool[psifac->getName()]=psifac;
    } else {
      WARNMSG("  " << psifac->getName() << " exists. Ignore addition")
    }
  }

  xmlNodePtr WaveFunctionPool::getWaveFunctionNode(const string& id) {
    if(myPool.empty()) return NULL;

    map<string,WaveFunctionFactory*>::iterator it(myPool.find(id));
    if(it == myPool.end()) {
      return (*myPool.begin()).second->myNode;
    } else {
      return (*it).second->myNode;
    }
  }

  bool WaveFunctionPool::put(std::istream& is) {
    return true;
  }

  bool WaveFunctionPool::get(std::ostream& os) const {
    return true;
  }

  void WaveFunctionPool::reset() {
 
  }

}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
